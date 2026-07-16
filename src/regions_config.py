"""Load fluid/solid region definitions from a CGNS-sidecar JSON.

Expected path: ``<same-basename-as-cgns>.json`` (e.g. ``foo.cgns`` → ``foo.json``).

**Canonical (minimal) format** — see
``tests/laptop_thermal_steady_scaled_v3_orig_BCs_fix.json``::

    {
      "fluid_regions": [
        "laptop_3d_geom.air.air_domain",
        "FPHPARTS.rotation1",
        "FPHPARTS.rotation2"
      ],
      "solid_regions": [
        "laptop_3d_geom.fan2.case2",
        "solid_region.Cu_block"
      ]
    }

Each entry is a CGNS zone name (matched flexibly). **All fluid zones are
placed into one OpenFOAM region named ``air``**
(``constant/air/polyMesh``) — zones are concatenated into that directory,
not renamed to a new ``fluid`` region. Each solid zone stays its own region
(sanitized zone name).

Optional foam2thermal-style ``regions`` with ``name`` / ``type`` / ``cellZones``
is still accepted; fluids are likewise coalesced into ``air``.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .topology import _sanitize_patch_name


@dataclass
class RegionSpec:
    name: str
    region_type: str  # "fluid" | "solid"
    cell_zones: list[str] = field(default_factory=list)


@dataclass
class RegionsConfig:
    path: Path
    specs: list[RegionSpec]
    #: CGNS zone name → (OpenFOAM region name, fluid|solid)
    zone_map: dict[str, tuple[str, str]]

    def foam_name_for(self, zone_name: str) -> str | None:
        hit = self.zone_map.get(zone_name)
        return hit[0] if hit else None

    def region_type_for(self, zone_name: str) -> str | None:
        hit = self.zone_map.get(zone_name)
        return hit[1] if hit else None

    def fluid_regions(self) -> list[str]:
        seen: list[str] = []
        for s in self.specs:
            if s.region_type == "fluid" and s.name not in seen:
                seen.append(s.name)
        return seen

    def solid_regions(self) -> list[str]:
        seen: list[str] = []
        for s in self.specs:
            if s.region_type == "solid" and s.name not in seen:
                seen.append(s.name)
        return seen


def sidecar_json_path(cgns_path: str | Path) -> Path:
    """``foo.cgns`` / ``foo.CGNS`` → ``foo.json`` beside the file."""
    p = Path(cgns_path)
    return p.with_suffix(".json")


def find_regions_json(cgns_path: str | Path) -> Path | None:
    cand = sidecar_json_path(cgns_path)
    return cand if cand.is_file() else None


def _norm(name: str) -> str:
    return _sanitize_patch_name(name).lower()


def _tokens(name: str) -> set[str]:
    parts = re.split(r"[._\-\s]+", name.lower())
    noise = {
        "partsurface", "laptop", "3d", "geom", "solid", "region",
        "fphparts", "domain", "block",
    }
    return {p for p in parts if p and not p.isdigit() and p not in noise}


def match_cell_zone_to_cgns(
    cell_zone: str,
    zone_names: list[str],
) -> str | None:
    """Map a JSON zone / cellZone entry onto a CGNS zone name."""
    if cell_zone in zone_names:
        return cell_zone
    n_cz = _norm(cell_zone)
    by_norm = {_norm(z): z for z in zone_names}
    if n_cz in by_norm:
        return by_norm[n_cz]

    for z in zone_names:
        nz = _norm(z)
        if nz.endswith(n_cz) or n_cz.endswith(nz):
            return z
        if n_cz in nz or nz in n_cz:
            return z

    tz = _tokens(cell_zone)
    best: tuple[int, str] | None = None
    for z in zone_names:
        score = len(tz & _tokens(z))
        if score <= 0:
            continue
        if best is None or score > best[0]:
            best = (score, z)
    return best[1] if best else None


_TRAILING_COMMA_RE = re.compile(r",(\s*[}\]])")


def _loads_json_relaxed(text: str) -> Any:
    """``json.loads`` with trailing-comma tolerance (common in hand-edited files)."""
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        cleaned = _TRAILING_COMMA_RE.sub(r"\1", text)
        return json.loads(cleaned)


# OpenFOAM directory / regionProperties name for all fluid zones' polyMesh.
MERGED_FLUID_REGION = "air"


def _dedupe_preserve(items: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for x in items:
        if x in seen:
            continue
        seen.add(x)
        out.append(x)
    return out


def _merge_fluid_specs(specs: list[RegionSpec]) -> list[RegionSpec]:
    """Collapse every fluid RegionSpec into one region named ``air``."""
    fluid_zones: list[str] = []
    solids: list[RegionSpec] = []
    for s in specs:
        if s.region_type == "fluid":
            fluid_zones.extend(s.cell_zones)
        else:
            solids.append(s)
    out: list[RegionSpec] = []
    fluid_zones = _dedupe_preserve(fluid_zones)
    if fluid_zones:
        out.append(
            RegionSpec(
                name=MERGED_FLUID_REGION,
                region_type="fluid",
                cell_zones=fluid_zones,
            )
        )
    out.extend(solids)
    return out


def _specs_from_fluid_solid_lists(data: dict[str, Any]) -> list[RegionSpec]:
    """Minimal format: ``fluid_regions`` / ``solid_regions`` (or ``fluid`` / ``solid``)."""
    fluid = data.get("fluid_regions")
    if fluid is None:
        fluid = data.get("fluid")
    solid = data.get("solid_regions")
    if solid is None:
        solid = data.get("solid")
    if not fluid and not solid:
        return []
    specs: list[RegionSpec] = []
    fluid_zones: list[str] = []
    for z in fluid or []:
        z = str(z).strip()
        if z:
            fluid_zones.append(z)
    if fluid_zones:
        specs.append(
            RegionSpec(
                name=MERGED_FLUID_REGION,
                region_type="fluid",
                cell_zones=_dedupe_preserve(fluid_zones),
            )
        )
    for z in solid or []:
        z = str(z).strip()
        if not z:
            continue
        specs.append(
            RegionSpec(
                name=_sanitize_patch_name(z),
                region_type="solid",
                cell_zones=[z],
            )
        )
    return specs


def _specs_from_regions_list(data: dict[str, Any]) -> list[RegionSpec]:
    regions = data.get("regions")
    if not isinstance(regions, list):
        return []
    specs: list[RegionSpec] = []
    for item in regions:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name", "")).strip()
        rtype = str(item.get("type", "")).strip().lower()
        if not name or rtype not in ("fluid", "solid"):
            continue
        cz = item.get("cellZones") or item.get("cell_zones") or [name]
        if isinstance(cz, str):
            cz = [cz]
        specs.append(
            RegionSpec(
                name=_sanitize_patch_name(name),
                region_type=rtype,
                cell_zones=[str(c) for c in cz],
            )
        )
    return specs


def _specs_from_regions_dict(data: dict[str, Any]) -> list[RegionSpec]:
    regions = data.get("regions")
    if not isinstance(regions, dict):
        return []
    specs: list[RegionSpec] = []
    for zname, rtype in regions.items():
        rt = str(rtype).strip().lower()
        if rt not in ("fluid", "solid"):
            continue
        specs.append(
            RegionSpec(
                name=_sanitize_patch_name(str(zname)),
                region_type=rt,
                cell_zones=[str(zname)],
            )
        )
    return specs


def _parse_specs(data: dict[str, Any]) -> list[RegionSpec]:
    # Prefer the minimal fluid_regions / solid_regions layout.
    specs = _specs_from_fluid_solid_lists(data)
    if not specs:
        specs = _specs_from_regions_list(data)
    if not specs:
        specs = _specs_from_regions_dict(data)
    # Always coalesce fluids into constant/air/polyMesh.
    return _merge_fluid_specs(specs)


def load_regions_config(
    path: str | Path,
    zone_names: list[str],
) -> RegionsConfig:
    """Parse *path* and bind JSON zone lists to CGNS *zone_names*."""
    p = Path(path)
    data = _loads_json_relaxed(p.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"regions JSON must be an object: {p}")

    specs = _parse_specs(data)
    if not specs:
        raise ValueError(
            f"no fluid/solid regions found in {p}; "
            "expected 'fluid_regions' / 'solid_regions' "
            "(see tests/laptop_thermal_steady_scaled_v3_orig_BCs_fix.json)"
        )

    zone_map: dict[str, tuple[str, str]] = {}
    unmatched: list[str] = []
    for spec in specs:
        for cz in spec.cell_zones:
            hit = match_cell_zone_to_cgns(cz, zone_names)
            if hit is None:
                unmatched.append(f"{spec.name}/{cz}")
                continue
            prev = zone_map.get(hit)
            if prev and prev[0] != spec.name:
                raise ValueError(
                    f"CGNS zone {hit!r} mapped to both {prev[0]!r} and "
                    f"{spec.name!r} in {p}"
                )
            zone_map[hit] = (spec.name, spec.region_type)

    if not zone_map:
        raise ValueError(
            f"no CGNS zones matched entries in {p}; "
            f"zones={zone_names!r}, unmatched={unmatched!r}"
        )

    return RegionsConfig(path=p, specs=specs, zone_map=zone_map)


def load_sidecar_regions(
    cgns_path: str | Path,
    zone_names: list[str],
    *,
    required: bool = False,
) -> RegionsConfig | None:
    """Load ``<cgns>.json`` if present (required for ``--cht`` / ``--cht-direct``)."""
    path = find_regions_json(cgns_path)
    if path is None:
        if required:
            expect = sidecar_json_path(cgns_path)
            raise FileNotFoundError(
                f"CHT mode requires a regions JSON beside the CGNS file: "
                f"{expect} "
                "(minimal format: fluid_regions / solid_regions; "
                "see tests/laptop_thermal_steady_scaled_v3_orig_BCs_fix.json)"
            )
        return None
    return load_regions_config(path, zone_names)
