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

Optional MRF (written to ``constant/air/MRFProperties``)::

    "mrf_regions": [
      {
        "cellZone": "FPHPARTS.rotation1",
        "origin": [-0.0678, -0.003, 0.081],
        "axis": [0.0, 1.0, 0.0],
        "omega": 100
      }
    ]

``origin`` may be ``"centroid"`` (zone vertex mean). foam2thermal nested
``regions[].mrf`` is also accepted.

Other optional keys:

* ``"g"`` / ``"gravity"``: gravity vector, e.g. ``[0, -9.81, 0]``
  (written to ``constant/g``; default ``(0 0 -9.81)``).
* ``"materials"``: per-region property overrides.  Keys match regions the
  same way ``heat_sources`` does; solid values accept ``rho`` / ``Cp`` /
  ``kappa`` / ``molWeight``, the fluid region ``air`` accepts ``mu`` /
  ``Pr`` / ``Cp`` / ``molWeight``::

      "materials": {
        "solid_region.Cu_block": {"rho": 8960, "Cp": 385, "kappa": 390},
        "air": {"mu": 1.846e-5, "Pr": 0.706, "Cp": 1006.43}
      }

* ``"external_convection"``: outer-wall convection BC
  (``externalWallHeatFluxTemperature``, mode ``coefficient``)::

      "external_convection": {"patches": ["Cover_outer", ".*_outer"],
                              "Ta": 300, "h": 8}

  ``patches`` entries are regexes matched against generated patch names.
* ``"initial_conditions"``: ``{"T": 300.0, "p": 101325.0}``.
* ``"heat_sources"``: ``{"<region or zone key>": <total watts>}``; multiple
  keys mapping to one region are summed.
* ``"n_procs"`` (or ``"parallel": {"nProcs": N}``): MPI ranks for
  ``decomposeParDict`` and ``Allrun`` (default 8).
* ``"endTime"`` / ``"writeInterval"`` / ``"purgeWrite"``: controlDict
  overrides (defaults 500 / 50 / 0).
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .topology import _sanitize_patch_name


@dataclass
class RegionSpec:
    name: str
    region_type: str  # "fluid" | "solid"
    cell_zones: list[str] = field(default_factory=list)


@dataclass
class MrfRegionSpec:
    """One rotating cellZone entry for OpenFOAM MRFProperties."""

    cell_zone: str  # matched CGNS zone name
    foam_cell_zone: str  # sanitized name used in polyMesh/cellZones
    omega: float
    axis: tuple[float, float, float]
    #: Explicit origin, or None → compute zone centroid at write time
    origin: tuple[float, float, float] | None = None
    non_rotating_patches: list[str] | None = None

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["axis"] = list(self.axis)
        if self.origin is not None:
            d["origin"] = list(self.origin)
        return d


@dataclass
class HeatSourceSpec:
    """Volumetric heat source for a solid region (watts)."""
    region_name: str  # sanitized OpenFOAM region name
    power: float       # total power in watts


@dataclass
class MaterialSpec:
    """Per-region material / thermophysical property overrides.

    Solid regions use ``rho`` / ``cp`` / ``kappa``; the fluid region uses
    ``mu`` / ``pr`` / ``cp``.  ``None`` keeps the built-in default.
    """

    region_name: str  # sanitized OpenFOAM region name
    rho: float | None = None
    cp: float | None = None
    kappa: float | None = None
    mu: float | None = None
    pr: float | None = None
    mol_weight: float | None = None


@dataclass
class ExternalConvectionSpec:
    """Outer-wall convection BC (externalWallHeatFluxTemperature)."""

    patterns: list[str]  # regexes matched against generated patch names
    ta: float = 300.0    # ambient temperature [K]
    h: float = 10.0      # heat transfer coefficient [W/(m2 K)]

    def matches(self, patch_name: str) -> bool:
        return any(re.search(p, patch_name) for p in self.patterns)


@dataclass
class RegionsConfig:
    path: Path
    specs: list[RegionSpec]
    #: CGNS zone name → (OpenFOAM region name, fluid|solid)
    zone_map: dict[str, tuple[str, str]]
    mrf_regions: list[MrfRegionSpec] = field(default_factory=list)
    heat_sources: list[HeatSourceSpec] = field(default_factory=list)
    materials: list[MaterialSpec] = field(default_factory=list)
    external_convection: ExternalConvectionSpec | None = None
    gravity: tuple[float, float, float] | None = None
    initial_t: float | None = None
    initial_p: float | None = None
    n_procs: int | None = None
    end_time: int | None = None
    write_interval: int | None = None
    purge_write: int | None = None

    def material_for(self, foam_name: str) -> MaterialSpec | None:
        for m in self.materials:
            if m.region_name == foam_name:
                return m
        return None

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


def _as_vec3(value: Any, *, label: str) -> tuple[float, float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{label} must be a length-3 list, got {value!r}")
    return (float(value[0]), float(value[1]), float(value[2]))


def _default_mrf_axis(zone_name: str) -> tuple[float, float, float]:
    name = zone_name.lower()
    if "rotation1" in name:
        return (0.0, 1.0, 0.0)
    if "rotation2" in name:
        return (0.0, -1.0, 0.0)
    return (0.0, 0.0, 1.0)


def _parse_mrf_regions(
    data: dict[str, Any],
    zone_names: list[str],
) -> list[MrfRegionSpec]:
    """Parse ``mrf_regions`` and/or foam2thermal ``regions[].mrf``."""
    out: list[MrfRegionSpec] = []
    seen_cz: set[str] = set()

    def _add(
        cz_raw: str,
        *,
        omega: float,
        axis: Any,
        origin: Any,
        non_rotating: list[str] | None,
    ) -> None:
        hit = match_cell_zone_to_cgns(str(cz_raw), zone_names)
        if hit is None:
            raise ValueError(f"MRF cellZone {cz_raw!r} not found in CGNS zones")
        if hit in seen_cz:
            return
        seen_cz.add(hit)
        if axis is None:
            ax = _default_mrf_axis(hit)
        else:
            ax = _as_vec3(axis, label=f"mrf axis for {cz_raw}")
        org: tuple[float, float, float] | None
        if origin is None or origin == "centroid":
            org = None
        else:
            org = _as_vec3(origin, label=f"mrf origin for {cz_raw}")
        out.append(
            MrfRegionSpec(
                cell_zone=hit,
                foam_cell_zone=_sanitize_patch_name(hit),
                omega=float(omega),
                axis=ax,
                origin=org,
                non_rotating_patches=list(non_rotating) if non_rotating else None,
            )
        )

    for item in data.get("mrf_regions") or []:
        if not isinstance(item, dict):
            continue
        cz = item.get("cellZone") or item.get("cell_zone") or item.get("zone")
        if not cz:
            continue
        _add(
            str(cz),
            omega=float(item.get("omega", 100)),
            axis=item.get("axis"),
            origin=item.get("origin", "centroid"),
            non_rotating=item.get("nonRotatingPatches")
            or item.get("non_rotating_patches"),
        )

    # foam2thermal: regions[].mrf { cellZones, omega, origin, axis/axes }
    for reg in data.get("regions") or []:
        if not isinstance(reg, dict):
            continue
        mrf = reg.get("mrf")
        if not isinstance(mrf, dict):
            continue
        zones = mrf.get("cellZones") or mrf.get("cell_zones") or []
        if isinstance(zones, str):
            zones = [zones]
        omega = float(mrf.get("omega", 100))
        origin = mrf.get("origin", "centroid")
        axes_cfg = mrf.get("axes")
        default_axis = mrf.get("axis")
        nr = mrf.get("nonRotatingPatches") or mrf.get("non_rotating_patches")
        for i, cz in enumerate(zones):
            axis = default_axis
            if isinstance(axes_cfg, dict):
                axis = axes_cfg.get(cz, axis)
            elif isinstance(axes_cfg, list) and i < len(axes_cfg):
                axis = axes_cfg[i]
            org = origin
            if isinstance(origin, list) and origin and isinstance(origin[0], (list, tuple)):
                org = origin[i] if i < len(origin) else origin[0]
            _add(
                str(cz),
                omega=omega,
                axis=axis,
                origin=org,
                non_rotating=nr if isinstance(nr, list) else None,
            )

    return out


def _region_key_matcher(
    specs: list[RegionSpec],
    *,
    solids_only: bool = True,
):
    """Build a ``key → RegionSpec`` matcher for JSON region references.

    Matches by sanitized foam name, raw/sanitized cell-zone name, and a
    last-resort substring match (same rules for ``heat_sources`` and
    ``materials``).
    """
    by_foam_name: dict[str, RegionSpec] = {}
    by_zone: dict[str, RegionSpec] = {}
    for s in specs:
        if solids_only and s.region_type != "solid":
            continue
        by_foam_name[s.name] = s
        for cz in s.cell_zones:
            by_zone[_sanitize_patch_name(cz).lower()] = s
            by_zone[cz.lower()] = s

    def match(key: str) -> RegionSpec | None:
        k = key.strip()
        if not k:
            return None
        spec = by_foam_name.get(_sanitize_patch_name(k)) or by_foam_name.get(k)
        if spec is None:
            spec = by_zone.get(_sanitize_patch_name(k).lower())
        if spec is None:
            spec = by_zone.get(k.lower())
        if spec is None:
            for fname, s in by_foam_name.items():
                if k.lower() in fname.lower() or fname.lower() in k.lower():
                    return s
        return spec

    return match


def _parse_heat_sources(
    data: dict[str, Any],
    specs: list[RegionSpec],
) -> list[HeatSourceSpec]:
    """Parse ``heat_sources`` from JSON.

    Format::

        "heat_sources": {
            "laptop_3d_geom.solid_region.CPU": 20,
            "solid_region.Cu_block": 15
        }

    Keys are matched to solid RegionSpec names (sanitized) or their CGNS
    cell-zone names.  Values are total power in watts; several keys may
    resolve to the same region (powers are summed by the writer).
    """
    raw = data.get("heat_sources")
    if not isinstance(raw, dict):
        return []
    match = _region_key_matcher(specs, solids_only=True)

    out: list[HeatSourceSpec] = []
    for key, val in raw.items():
        k = str(key).strip()
        if not k:
            continue
        power = float(val)
        if power <= 0:
            continue
        spec = match(k)
        if spec is None:
            raise ValueError(
                f"heat_sources key {k!r} does not match any solid region"
            )
        out.append(HeatSourceSpec(region_name=spec.name, power=power))
    return out


def _parse_materials(
    data: dict[str, Any],
    specs: list[RegionSpec],
) -> list[MaterialSpec]:
    """Parse ``materials`` property overrides.

    Format::

        "materials": {
            "solid_region.Cu_block": {"rho": 8960, "Cp": 385, "kappa": 390},
            "air": {"mu": 1.846e-5, "Pr": 0.706, "Cp": 1006.43}
        }

    Keys match regions like ``heat_sources`` (both fluid and solid
    regions are eligible).  Any subset of properties may be given;
    missing ones keep the built-in defaults.
    """
    raw = data.get("materials")
    if not isinstance(raw, dict):
        return []
    match = _region_key_matcher(specs, solids_only=False)

    out: list[MaterialSpec] = []
    for key, val in raw.items():
        k = str(key).strip()
        if not k or not isinstance(val, dict):
            continue
        spec = match(k)
        if spec is None:
            raise ValueError(
                f"materials key {k!r} does not match any region"
            )

        def _f(*names: str) -> float | None:
            for n in names:
                if n in val and val[n] is not None:
                    return float(val[n])
            return None

        out.append(
            MaterialSpec(
                region_name=spec.name,
                rho=_f("rho", "density"),
                cp=_f("Cp", "cp"),
                kappa=_f("kappa", "k", "conductivity"),
                mu=_f("mu", "viscosity"),
                pr=_f("Pr", "pr"),
                mol_weight=_f("molWeight", "mol_weight"),
            )
        )
    return out


def _parse_external_convection(data: dict[str, Any]) -> ExternalConvectionSpec | None:
    """Parse ``external_convection`` outer-wall convection BC.

    Format::

        "external_convection": {
            "patches": ["Cover_outer", ".*_outer"],
            "Ta": 300,
            "h": 8
        }

    ``patches`` entries are regexes matched against generated patch names.
    """
    raw = data.get("external_convection")
    if not isinstance(raw, dict):
        return None
    pats = raw.get("patches") or raw.get("patterns") or []
    if isinstance(pats, str):
        pats = [pats]
    pats = [str(p) for p in pats if str(p).strip()]
    if not pats:
        return None
    return ExternalConvectionSpec(
        patterns=pats,
        ta=float(raw.get("Ta", raw.get("ta", 300.0))),
        h=float(raw.get("h", raw.get("htc", 10.0))),
    )


def _parse_optional_vec3(data: dict[str, Any], *keys: str) -> tuple[float, float, float] | None:
    for k in keys:
        v = data.get(k)
        if v is not None:
            return _as_vec3(v, label=k)
    return None


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

    mrf_regions = _parse_mrf_regions(data, zone_names)
    heat_sources = _parse_heat_sources(data, specs)
    materials = _parse_materials(data, specs)
    external_convection = _parse_external_convection(data)
    gravity = _parse_optional_vec3(data, "g", "gravity")

    initial_t: float | None = None
    initial_p: float | None = None
    ic = data.get("initial_conditions") or data.get("initial")
    if isinstance(ic, dict):
        if ic.get("T") is not None:
            initial_t = float(ic["T"])
        if ic.get("p") is not None:
            initial_p = float(ic["p"])

    n_procs: int | None = None
    par = data.get("parallel")
    if isinstance(par, dict) and par.get("nProcs") is not None:
        n_procs = int(par["nProcs"])
    elif isinstance(par, (int, float)):
        n_procs = int(par)
    if data.get("n_procs") is not None:
        n_procs = int(data["n_procs"])
    if data.get("nProcs") is not None:
        n_procs = int(data["nProcs"])

    def _opt_int(*keys: str) -> int | None:
        for k in keys:
            if data.get(k) is not None:
                return int(data[k])
        return None

    return RegionsConfig(
        path=p,
        specs=specs,
        zone_map=zone_map,
        mrf_regions=mrf_regions,
        heat_sources=heat_sources,
        materials=materials,
        external_convection=external_convection,
        gravity=gravity,
        initial_t=initial_t,
        initial_p=initial_p,
        n_procs=n_procs,
        end_time=_opt_int("endTime", "end_time"),
        write_interval=_opt_int("writeInterval", "write_interval"),
        purge_write=_opt_int("purgeWrite", "purge_write"),
    )


def load_sidecar_regions(
    cgns_path: str | Path,
    zone_names: list[str],
    *,
    required: bool = False,
) -> RegionsConfig | None:
    """Load ``<cgns>.json`` if present (required for ``--cht-direct``)."""
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
