"""Scan CGNS zones for fluid/solid regions and coupling interface pairs.

Couplings are inferred from FaceCenter BC geometry: same-named (or geometrically
coincident) boundary faces shared by two zones form an interface.  Zone type is
heuristic (``solid_region`` / ``solid.`` → solid, else fluid) and can be refined
via explicit solid/fluid name patterns.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any

import numpy as np

from .reader import CGNSCase, CGNSZone
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
            "fluid_regions": self.fluid_regions,
            "solid_regions": self.solid_regions,
            "regions": [asdict(r) for r in self.regions],
            "couplings_by_kind": by_kind,
            "coupling_pairs": [c.to_dict() for c in self.couplings],
            "n_couplings": len(self.couplings),
        }


_AMI_RE = re.compile(r"ami_rot\d+|.*rotation\d*", re.I)


def classify_zone_type(
    zone_name: str,
    *,
    solid_patterns: list[str] | None = None,
    fluid_patterns: list[str] | None = None,
) -> str:
    """Return ``\"solid\"`` or ``\"fluid\"`` for a CGNS zone name."""
    name = zone_name.lower()
    for pat in solid_patterns or []:
        if re.search(pat, zone_name, re.I) or re.search(pat, name, re.I):
            return "solid"
    for pat in fluid_patterns or []:
        if re.search(pat, zone_name, re.I) or re.search(pat, name, re.I):
            return "fluid"
    return "solid" if _zone_is_solid(zone_name) else "fluid"


def _looks_ami(name: str) -> bool:
    return bool(_AMI_RE.fullmatch(name) or _AMI_RE.search(name))


def _classify_pair(
    type_a: str,
    type_b: str,
    bc_a: str,
    bc_b: str,
    zone_a: str,
    zone_b: str,
) -> tuple[CouplingKind, CouplingMethod]:
    ami = (
        _looks_ami(bc_a)
        or _looks_ami(bc_b)
        or _looks_ami(zone_a)
        or _looks_ami(zone_b)
        or "rotation" in zone_a.lower()
        or "rotation" in zone_b.lower()
    )
    if type_a == "fluid" and type_b == "fluid":
        kind = CouplingKind.FLUID_FLUID
        method = CouplingMethod.CYCLIC_AMI if ami else CouplingMethod.MAPPED_WALL
    elif type_a == "solid" and type_b == "solid":
        kind = CouplingKind.SOLID_SOLID
        method = CouplingMethod.MAPPED_WALL
    else:
        kind = CouplingKind.FLUID_SOLID
        method = CouplingMethod.MAPPED_WALL
    return kind, method


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


def scan_couplings(
    case: CGNSCase,
    *,
    source: str = "",
    point_tol: float = 1e-4,
    solid_patterns: list[str] | None = None,
    fluid_patterns: list[str] | None = None,
    trim: bool = True,
    min_overlap_faces: int = 1,
) -> CouplingReport:
    """Scan a loaded CGNS case for regions and coupling interface pairs.

    Steps:

    1. Build per-zone topology and optionally apply cross-zone BC trim.
    2. Classify each zone as fluid / solid.
    3. Pair FaceCenter BCs that share geometric faces across distinct zones.
    """
    if not case.zones:
        raise ValueError("CGNS case contains no zones")

    zone_topos = [_build_zone_topology(z) for z in case.zones]
    if trim:
        _trim_cross_zone_bc_overlaps(case.zones, zone_topos, point_tol=point_tol)

    regions: list[RegionInfo] = []
    zone_types: list[str] = []
    for zone, zt in zip(case.zones, zone_topos):
        rtype = classify_zone_type(
            zone.name,
            solid_patterns=solid_patterns,
            fluid_patterns=fluid_patterns,
        )
        zone_types.append(rtype)
        bc_names = sorted(zt.bc_face_lists.keys())
        regions.append(
            RegionInfo(
                zone_name=zone.name,
                foam_name=_sanitize_patch_name(zone.name),
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
            # Same-named BC preferred; otherwise require substantial geometric overlap
            # and matching stem after sanitize.
            same_name = _sanitize_patch_name(bc_a) == _sanitize_patch_name(bc_b)
            inter = keys_a & keys_b
            n_ov = len(inter)
            if n_ov < min_overlap_faces:
                continue
            if not same_name:
                # Geometric-only match: both sides must be mostly overlapping
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
            kind, method = _classify_pair(type_a, type_b, bc_a, bc_b, zone_a, zone_b)

            # Stable master/slave: fluid before solid; else lexicographic zone name
            if type_a == type_b:
                if zone_a <= zone_b:
                    master_z, slave_z, master_bc, slave_bc = zone_a, zone_b, bc_a, bc_b
                else:
                    master_z, slave_z, master_bc, slave_bc = zone_b, zone_a, bc_b, bc_a
            elif type_a == "fluid":
                master_z, slave_z, master_bc, slave_bc = zone_a, zone_b, bc_a, bc_b
            else:
                master_z, slave_z, master_bc, slave_bc = zone_b, zone_a, bc_b, bc_a

            note = "same_bc_name" if same_name else "geometric_match"
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
                )
            )

    couplings.sort(key=lambda c: (c.kind.value, -c.n_faces, c.master_zone, c.slave_zone))

    fluid = [r.foam_name for r in regions if r.region_type == "fluid"]
    solid = [r.foam_name for r in regions if r.region_type == "solid"]
    return CouplingReport(
        regions=regions,
        couplings=couplings,
        fluid_regions=fluid,
        solid_regions=solid,
        source=source,
    )


def format_coupling_summary(report: CouplingReport) -> str:
    """Human-readable multi-line summary for CLI stdout."""
    lines = [
        f"source: {report.source or '(in-memory)'}",
        f"regions: {len(report.regions)} "
        f"(fluid={len(report.fluid_regions)}, solid={len(report.solid_regions)})",
        "",
        "=== Regions ===",
    ]
    for r in report.regions:
        lines.append(
            f"  [{r.region_type:5}] {r.foam_name}  "
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
                f"  {c.master_zone}:{c.master_bc}  <->  "
                f"{c.slave_zone}:{c.slave_bc}  "
                f"nFaces={c.n_faces}  method={c.method.value}"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"
