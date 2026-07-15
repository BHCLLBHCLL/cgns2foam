"""Scan CGNS zones for fluid/solid regions and coupling interface pairs.

Couplings are inferred from FaceCenter BC geometry: same-named (or geometrically
coincident) boundary faces shared by two zones form an interface.

When a sidecar ``<cgns>.json`` is supplied (required for ``--cht`` /
``--cht-direct``), zone types and OpenFOAM region names come from that file.
Interface methods are then forced by region types:

* fluid–fluid → ``cyclicAMI``
* fluid–solid / solid–solid → ``mappedWall``
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any

import numpy as np

from .reader import CGNSCase, CGNSZone
from .regions_config import RegionsConfig
from .topology import (
    _ZoneTopo,
    _build_zone_topology,
    _face_coord_key,
    _sanitize_patch_name,
    _trim_cross_zone_bc_overlaps,
    _zone_is_solid,
)


class CouplingKind(str, Enum):
    FLUID_FLUID = "fluid_fluid"
    FLUID_SOLID = "fluid_solid"
    SOLID_SOLID = "solid_solid"


class CouplingMethod(str, Enum):
    """Suggested OpenFOAM mesh-prep / BC method for the pair."""

    CYCLIC_AMI = "cyclicAMI"
    MAPPED_WALL = "mappedWall"
    STITCH = "stitch"


@dataclass
class RegionInfo:
    zone_name: str
    foam_name: str
    region_type: str  # "fluid" | "solid"
    n_cells: int
    n_vertices: int
    bc_names: list[str] = field(default_factory=list)


@dataclass
class CouplingPair:
    kind: CouplingKind
    method: CouplingMethod
    master_zone: str
    slave_zone: str
    master_bc: str
    slave_bc: str
    n_faces: int
    note: str = ""
    master_region: str = ""
    slave_region: str = ""

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["kind"] = self.kind.value
        d["method"] = self.method.value
        return d


@dataclass
class CouplingReport:
    regions: list[RegionInfo]
    couplings: list[CouplingPair]
    fluid_regions: list[str]
    solid_regions: list[str]
    source: str = ""
    regions_json: str = ""

    def to_dict(self) -> dict[str, Any]:
        by_kind: dict[str, list[dict]] = {
            CouplingKind.FLUID_FLUID.value: [],
            CouplingKind.FLUID_SOLID.value: [],
            CouplingKind.SOLID_SOLID.value: [],
        }
        for c in self.couplings:
            by_kind[c.kind.value].append(c.to_dict())
        return {
            "source": self.source,
            "regions_json": self.regions_json,
            "fluid_regions": self.fluid_regions,
            "solid_regions": self.solid_regions,
            "regions": [asdict(r) for r in self.regions],
            "couplings_by_kind": by_kind,
            "coupling_pairs": [c.to_dict() for c in self.couplings],
            "n_couplings": len(self.couplings),
        }


def classify_zone_type(
    zone_name: str,
    *,
    solid_patterns: list[str] | None = None,
    fluid_patterns: list[str] | None = None,
    regions_config: RegionsConfig | None = None,
) -> str:
    """Return ``\"solid\"`` or ``\"fluid\"`` for a CGNS zone name."""
    if regions_config is not None:
        hit = regions_config.region_type_for(zone_name)
        if hit is not None:
            return hit
    name = zone_name.lower()
    for pat in solid_patterns or []:
        if re.search(pat, zone_name, re.I) or re.search(pat, name, re.I):
            return "solid"
    for pat in fluid_patterns or []:
        if re.search(pat, zone_name, re.I) or re.search(pat, name, re.I):
            return "fluid"
    return "solid" if _zone_is_solid(zone_name) else "fluid"


def foam_name_for_zone(
    zone_name: str,
    *,
    regions_config: RegionsConfig | None = None,
) -> str:
    if regions_config is not None:
        hit = regions_config.foam_name_for(zone_name)
        if hit is not None:
            return hit
    return _sanitize_patch_name(zone_name)


def classify_interface_method(
    type_a: str,
    type_b: str,
    *,
    same_foam_region: bool = False,
) -> tuple[CouplingKind, CouplingMethod]:
    """Map region types to coupling kind / OpenFOAM interface method.

    Rule (CHT):

    * fluid–fluid → ``cyclicAMI`` (including AMI pairs inside a merged fluid
      region such as air_domain ↔ rotation*)
    * fluid–solid / solid–solid → ``mappedWall`` when the sides belong to
      different OpenFOAM regions; same-region solid interfaces → ``stitch``
    """
    if type_a == "fluid" and type_b == "fluid":
        return CouplingKind.FLUID_FLUID, CouplingMethod.CYCLIC_AMI
    if type_a == "solid" and type_b == "solid":
        kind = CouplingKind.SOLID_SOLID
        if same_foam_region:
            return kind, CouplingMethod.STITCH
        return kind, CouplingMethod.MAPPED_WALL
    return CouplingKind.FLUID_SOLID, CouplingMethod.MAPPED_WALL


def _classify_pair(
    type_a: str,
    type_b: str,
    bc_a: str,
    bc_b: str,
    zone_a: str,
    zone_b: str,
    *,
    foam_a: str = "",
    foam_b: str = "",
) -> tuple[CouplingKind, CouplingMethod]:
    _ = (bc_a, bc_b, zone_a, zone_b)
    same = bool(foam_a and foam_b and foam_a == foam_b)
    return classify_interface_method(type_a, type_b, same_foam_region=same)


def _bc_face_keys(
    zone: CGNSZone,
    zt: _ZoneTopo,
    bc_name: str,
    tol: float,
) -> set[tuple[tuple[int, ...], ...]]:
    ids = zt.bc_face_lists.get(bc_name)
    if ids is None or ids.size == 0:
        return set()
    pts = zone.coords
    return {
        _face_coord_key(pts, zt.face_offsets, zt.face_vertices, int(fi), tol)
        for fi in ids
    }


def prepare_zone_topos(
    case: CGNSCase,
    *,
    point_tol: float = 1e-4,
    trim: bool = True,
) -> list[_ZoneTopo]:
    """Build (and optionally BC-trim) per-zone topologies once."""
    if not case.zones:
        raise ValueError("CGNS case contains no zones")
    zone_topos = [_build_zone_topology(z) for z in case.zones]
    if trim:
        _trim_cross_zone_bc_overlaps(case.zones, zone_topos, point_tol=point_tol)
    return zone_topos


def scan_couplings(
    case: CGNSCase,
    *,
    source: str = "",
    point_tol: float = 1e-4,
    solid_patterns: list[str] | None = None,
    fluid_patterns: list[str] | None = None,
    regions_config: RegionsConfig | None = None,
    trim: bool = True,
    min_overlap_faces: int = 1,
    zone_topos: list[_ZoneTopo] | None = None,
) -> CouplingReport:
    """Scan a loaded CGNS case for regions and coupling interface pairs.

    Steps:

    1. Build per-zone topology and optionally apply cross-zone BC trim
       (skipped when *zone_topos* is supplied).
    2. Classify each zone as fluid / solid (JSON sidecar preferred).
    3. Pair FaceCenter BCs that share geometric faces across distinct zones.
    4. Assign interface method: fluid–fluid → cyclicAMI; else mappedWall
       (except same-region solid → stitch).
    """
    if not case.zones:
        raise ValueError("CGNS case contains no zones")

    if zone_topos is None:
        zone_topos = prepare_zone_topos(case, point_tol=point_tol, trim=trim)
    elif len(zone_topos) != len(case.zones):
        raise ValueError("zone_topos length must match number of CGNS zones")

    regions: list[RegionInfo] = []
    zone_types: list[str] = []
    foam_names: list[str] = []
    for zone, zt in zip(case.zones, zone_topos):
        rtype = classify_zone_type(
            zone.name,
            solid_patterns=solid_patterns,
            fluid_patterns=fluid_patterns,
            regions_config=regions_config,
        )
        fname = foam_name_for_zone(zone.name, regions_config=regions_config)
        zone_types.append(rtype)
        foam_names.append(fname)
        bc_names = sorted(zt.bc_face_lists.keys())
        regions.append(
            RegionInfo(
                zone_name=zone.name,
                foam_name=fname,
                region_type=rtype,
                n_cells=int(zt.n_cells),
                n_vertices=int(zt.n_vertices),
                bc_names=bc_names,
            )
        )

    # Entry: (zi, bc_name) → face keys
    entries: list[tuple[int, str, set]] = []
    for zi, (zone, zt) in enumerate(zip(case.zones, zone_topos)):
        for bc_name, ids in zt.bc_face_lists.items():
            if ids.size == 0:
                continue
            keys = _bc_face_keys(zone, zt, bc_name, point_tol)
            if keys:
                entries.append((zi, bc_name, keys))

    couplings: list[CouplingPair] = []
    seen: set[frozenset[str]] = set()

    for i in range(len(entries)):
        zi_a, bc_a, keys_a = entries[i]
        for j in range(i + 1, len(entries)):
            zi_b, bc_b, keys_b = entries[j]
            if zi_a == zi_b:
                continue
            same_name = _sanitize_patch_name(bc_a) == _sanitize_patch_name(bc_b)
            inter = keys_a & keys_b
            n_ov = len(inter)
            if n_ov < min_overlap_faces:
                continue
            if not same_name:
                ratio_a = n_ov / max(len(keys_a), 1)
                ratio_b = n_ov / max(len(keys_b), 1)
                if ratio_a < 0.9 or ratio_b < 0.9:
                    continue

            zone_a = case.zones[zi_a].name
            zone_b = case.zones[zi_b].name
            key = frozenset(
                {
                    f"{zone_a}/{bc_a}",
                    f"{zone_b}/{bc_b}",
                }
            )
            if key in seen:
                continue
            seen.add(key)

            type_a, type_b = zone_types[zi_a], zone_types[zi_b]
            foam_a, foam_b = foam_names[zi_a], foam_names[zi_b]
            kind, method = _classify_pair(
                type_a, type_b, bc_a, bc_b, zone_a, zone_b,
                foam_a=foam_a, foam_b=foam_b,
            )

            # Stable master/slave: fluid before solid; else lexicographic zone name
            if type_a == type_b:
                if zone_a <= zone_b:
                    master_z, slave_z, master_bc, slave_bc = zone_a, zone_b, bc_a, bc_b
                    master_r, slave_r = foam_a, foam_b
                else:
                    master_z, slave_z, master_bc, slave_bc = zone_b, zone_a, bc_b, bc_a
                    master_r, slave_r = foam_b, foam_a
            elif type_a == "fluid":
                master_z, slave_z, master_bc, slave_bc = zone_a, zone_b, bc_a, bc_b
                master_r, slave_r = foam_a, foam_b
            else:
                master_z, slave_z, master_bc, slave_bc = zone_b, zone_a, bc_b, bc_a
                master_r, slave_r = foam_b, foam_a

            note = "same_bc_name" if same_name else "geometric_match"
            if master_r == slave_r:
                note += ";same_region"
            couplings.append(
                CouplingPair(
                    kind=kind,
                    method=method,
                    master_zone=master_z,
                    slave_zone=slave_z,
                    master_bc=master_bc,
                    slave_bc=slave_bc,
                    n_faces=n_ov,
                    note=note,
                    master_region=master_r,
                    slave_region=slave_r,
                )
            )

    couplings.sort(key=lambda c: (c.kind.value, -c.n_faces, c.master_zone, c.slave_zone))

    if regions_config is not None:
        fluid = regions_config.fluid_regions()
        solid = regions_config.solid_regions()
        # Include any unmatched zones that still appear in the scan
        for r in regions:
            if r.region_type == "fluid" and r.foam_name not in fluid:
                fluid.append(r.foam_name)
            if r.region_type == "solid" and r.foam_name not in solid:
                solid.append(r.foam_name)
    else:
        fluid = []
        solid = []
        for r in regions:
            if r.region_type == "fluid" and r.foam_name not in fluid:
                fluid.append(r.foam_name)
            elif r.region_type == "solid" and r.foam_name not in solid:
                solid.append(r.foam_name)

    return CouplingReport(
        regions=regions,
        couplings=couplings,
        fluid_regions=fluid,
        solid_regions=solid,
        source=source,
        regions_json=str(regions_config.path) if regions_config else "",
    )


def format_coupling_summary(report: CouplingReport) -> str:
    """Human-readable multi-line summary for CLI stdout."""
    lines = [
        f"source: {report.source or '(in-memory)'}",
    ]
    if report.regions_json:
        lines.append(f"regions_json: {report.regions_json}")
    lines.extend(
        [
            f"regions: {len(report.regions)} CGNS zones → "
            f"OpenFOAM fluid={len(report.fluid_regions)}, "
            f"solid={len(report.solid_regions)}",
            "",
            "=== Regions (CGNS zone → OpenFOAM) ===",
        ]
    )
    for r in report.regions:
        lines.append(
            f"  [{r.region_type:5}] {r.zone_name}  →  {r.foam_name}  "
            f"cells={r.n_cells}  bcs={len(r.bc_names)}"
        )
    lines.append("")
    for kind in (
        CouplingKind.FLUID_FLUID,
        CouplingKind.FLUID_SOLID,
        CouplingKind.SOLID_SOLID,
    ):
        items = [c for c in report.couplings if c.kind == kind]
        lines.append(f"=== {kind.value} ({len(items)}) ===")
        if not items:
            lines.append("  (none)")
            continue
        for c in items:
            lines.append(
                f"  {c.master_region or c.master_zone}:{c.master_bc}  <->  "
                f"{c.slave_region or c.slave_zone}:{c.slave_bc}  "
                f"nFaces={c.n_faces}  method={c.method.value}"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"
