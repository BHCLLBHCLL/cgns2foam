"""One-step CGNS → multi-region chtMultiRegionSimpleFoam case writer.

Requires a sidecar ``<cgns>.json`` that lists fluid/solid OpenFOAM regions
and their CGNS ``cellZones``.  Zones that share one JSON region are merged
into a single ``constant/<region>/polyMesh``.

Coupling BCs follow region types:

* fluid–fluid → ``cyclicAMI`` (same mesh after merge; ``neighbourPatch``)
* fluid–solid / solid–solid → ``mappedWall`` (``sampleRegion`` / ``samplePatch``)
"""

from __future__ import annotations

import json
import os
import shutil
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from .cht_case import (
    _allclean,
    _control_dict_cht,
    _decompose_par_dict,
    _field_T,
    _field_U,
    _field_p,
    _field_p_rgh,
    _fv_options_fluid,
    _fv_options_solid_heat,
    _fv_schemes_fluid,
    _fv_schemes_solid,
    _fv_solution_fluid,
    _fv_solution_solid,
    _gravity,
    _radiation_none,
    _region_properties,
    _thermophysical_fluid,
    _thermophysical_solid,
    _turbulence_laminar,
    _write_text,
    mrf_non_rotating_patches,
    mrf_properties,
)
from .couplings import (
    CouplingMethod,
    CouplingReport,
    format_coupling_summary,
    prepare_zone_topos,
    scan_couplings,
)
from .reader import CGNSCase, read_cgns
from .regions_config import (
    MERGED_FLUID_REGION,
    MrfRegionSpec,
    RegionsConfig,
    HeatSourceSpec,
    load_sidecar_regions,
)
from .topology import (
    CellZone,
    Mesh,
    _ZoneTopo,
    _assemble_ordered_mesh,
    _bc_type_to_foam,
    _sanitize_patch_name,
    foam_patch_type_for_name,
)
from .writer import WriteOptions, write_poly_mesh


def _zone_vertex_centroid(zone) -> tuple[float, float, float]:
    c = zone.coords.mean(axis=0)
    return (float(c[0]), float(c[1]), float(c[2]))


def resolve_mrf_entries(
    mrf_specs: list[MrfRegionSpec],
    case: CGNSCase,
    patch_names: list[str],
) -> list[dict[str, Any]]:
    """Turn JSON MRF specs into MRFProperties entry dicts."""
    zone_by_name = {z.name: z for z in case.zones}
    default_nr = mrf_non_rotating_patches(patch_names)
    entries: list[dict[str, Any]] = []
    for i, spec in enumerate(mrf_specs):
        zone = zone_by_name.get(spec.cell_zone)
        if zone is None:
            raise ValueError(f"MRF cellZone {spec.cell_zone!r} missing from CGNS case")
        origin = spec.origin if spec.origin is not None else _zone_vertex_centroid(zone)
        nr = spec.non_rotating_patches if spec.non_rotating_patches is not None else default_nr
        entries.append(
            {
                "name": f"MRF{i + 1}" if len(mrf_specs) > 1 else "MRF",
                "cellZone": spec.foam_cell_zone,
                "origin": origin,
                "axis": list(spec.axis),
                "omega": spec.omega,
                "nonRotatingPatches": nr,
                "cgnsZone": spec.cell_zone,
            }
        )
    return entries



def zone_name_stem(zone_name: str) -> str:
    """Last ``.``-separated token of a CGNS zone / cellZone name (OF-safe).

    Examples: ``solid_region.Cu_block`` → ``Cu_block``,
    ``laptop_3d_geom.solid_region.CPU`` → ``CPU``.
    """
    tail = zone_name.rsplit(".", 1)[-1]
    return _sanitize_patch_name(tail)


def coupling_stem(zone_name: str, foam_name: str) -> str:
    """Short label for coupling patch pairing.

    Merged fluid region always uses ``air``; solids use :func:`zone_name_stem`.
    """
    if foam_name == MERGED_FLUID_REGION or foam_name == "air":
        return "air"
    return zone_name_stem(zone_name)


def coupling_patch_name(local_stem: str, remote_stem: str) -> str:
    """mappedWall inter-region patch name, e.g. ``CPU_to_Cover``."""
    return f"{local_stem}_to_{remote_stem}"


def ami_patch_name(local_zone: str, remote_zone: str) -> str:
    """cyclicAMI patch name from zone stems (no ``_to_``; field stubs use ami_)."""
    a = zone_name_stem(local_zone)
    b = zone_name_stem(remote_zone)
    return f"ami_{a}_{b}"


def _region_patch_plan(
    zone_indices: list[int],
    case: CGNSCase,
    zone_topos: list[_ZoneTopo],
    foam_name: str,
    report: CouplingReport,
    *,
    face_offset: list[int],
    default_exterior_name: str = "default_exterior",
) -> list[tuple]:
    """Build patch plan for a (possibly multi-zone) OpenFOAM region.

    Face ids are in the merged face index space (global within the region).
    Coupling patches use short stems (``CPU_to_Cover``); ``sampleRegion``
    still points at the full OpenFOAM region directory name.
    """
    zones_in = {case.zones[i].name for i in zone_indices}
    foam_by_zone = {r.zone_name: r.foam_name for r in report.regions}
    zi_by_name = {case.zones[i].name: i for i in zone_indices}

    assigned: dict[int, np.ndarray] = {
        i: np.zeros(zone_topos[i].n_faces, dtype=bool) for i in zone_indices
    }
    patches: list[tuple] = []
    used: set[str] = set()

    def global_ids(zi: int, local_ids: np.ndarray) -> np.ndarray:
        return local_ids.astype(np.int64) + face_offset[zi]

    # --- Couplings first ---
    for c in report.couplings:
        if c.master_zone in zones_in and c.slave_zone in zones_in:
            # Both sides in this merged region → cyclicAMI pair
            if c.method != CouplingMethod.CYCLIC_AMI:
                # same-region solid stitch: leave as ordinary walls for now
                continue
            pairs = (
                (c.master_zone, c.master_bc, c.slave_zone),
                (c.slave_zone, c.slave_bc, c.master_zone),
            )
            names = (
                ami_patch_name(c.master_zone, c.slave_zone),
                ami_patch_name(c.slave_zone, c.master_zone),
            )
            for (zname, bc, _remote), pname, nbr in zip(pairs, names, names[::-1]):
                zi = zi_by_name[zname]
                zt = zone_topos[zi]
                ids = zt.bc_face_lists.get(bc)
                if ids is None or ids.size == 0:
                    continue
                boundary_mask = zt.neighbour < 0
                keep = boundary_mask[ids] & ~assigned[zi][ids]
                ids = ids[keep]
                if ids.size == 0:
                    continue
                if pname in used:
                    for i, entry in enumerate(patches):
                        if entry[0] == pname:
                            patches[i] = (
                                pname,
                                entry[1],
                                np.unique(np.concatenate([entry[2], global_ids(zi, ids)])),
                                entry[3],
                            )
                            assigned[zi][ids] = True
                            break
                    continue
                used.add(pname)
                extras = {"neighbour_patch": nbr}
                patches.append((pname, "cyclicAMI", global_ids(zi, ids), extras))
                assigned[zi][ids] = True
            continue

        # One side in this region → mappedWall (or skip if stitch-only)
        if c.master_zone in zones_in:
            local_zone, local_bc, remote_zone = c.master_zone, c.master_bc, c.slave_zone
        elif c.slave_zone in zones_in:
            local_zone, local_bc, remote_zone = c.slave_zone, c.slave_bc, c.master_zone
        else:
            continue

        remote = foam_by_zone.get(remote_zone, _sanitize_patch_name(remote_zone))
        if remote == foam_name:
            continue

        if c.method == CouplingMethod.STITCH:
            continue

        zi = zi_by_name[local_zone]
        zt = zone_topos[zi]
        ids = zt.bc_face_lists.get(local_bc)
        if ids is None or ids.size == 0:
            continue
        boundary_mask = zt.neighbour < 0
        keep = boundary_mask[ids] & ~assigned[zi][ids]
        ids = ids[keep]
        if ids.size == 0:
            continue

        if c.method not in (
            CouplingMethod.CYCLIC_AMI,
            CouplingMethod.MAPPED_WALL,
        ):
            continue

        local_stem = coupling_stem(local_zone, foam_name)
        remote_stem = coupling_stem(remote_zone, remote)
        pname = coupling_patch_name(local_stem, remote_stem)
        extras = {
            "sample_mode": "nearestPatchFace",
            "sample_region": remote,
            "sample_patch": coupling_patch_name(remote_stem, local_stem),
        }

        if pname in used:
            for i, entry in enumerate(patches):
                if entry[0] == pname:
                    patches[i] = (
                        pname,
                        entry[1],
                        np.unique(np.concatenate([entry[2], global_ids(zi, ids)])),
                        entry[3],
                    )
                    assigned[zi][ids] = True
                    break
            continue
        used.add(pname)
        patches.append((pname, "mappedWall", global_ids(zi, ids), extras))
        assigned[zi][ids] = True

    # --- Remaining CGNS BCs ---
    for zi in zone_indices:
        zone = case.zones[zi]
        zt = zone_topos[zi]
        boundary_mask = zt.neighbour < 0
        for bc in zone.bcs:
            if bc.name not in zt.bc_face_lists:
                continue
            ids = zt.bc_face_lists[bc.name]
            keep = boundary_mask[ids] & ~assigned[zi][ids]
            ids = ids[keep]
            if ids.size == 0:
                continue
            base = _sanitize_patch_name(bc.name)
            pname = base
            i = 1
            while pname in used:
                pname = f"{base}_{i}"
                i += 1
            used.add(pname)
            patches.append(
                (
                    pname,
                    foam_patch_type_for_name(pname, bc.bc_type),
                    global_ids(zi, ids),
                    {"cgns_bc_type": bc.bc_type},
                )
            )
            assigned[zi][ids] = True

        remaining = np.where(boundary_mask & ~assigned[zi])[0]
        if remaining.size:
            pname = default_exterior_name
            i = 1
            while pname in used:
                pname = f"{default_exterior_name}_{i}"
                i += 1
            used.add(pname)
            # Name-based override: "open*" must be patch (total-pressure opening),
            # never wall — even when faces were not tagged as a CGNS BC.
            patches.append(
                (pname, foam_patch_type_for_name(pname, "wall"), global_ids(zi, remaining), {})
            )
            assigned[zi][remaining] = True

    return patches


# CGNS BCType strings that denote a flow opening (inlet / outlet / farfield).
# Field BCs for these patches follow the "open" total-pressure template.
_OPENING_BC_TOKENS = (
    "bcinflow",
    "bcoutflow",
    "bcfarfield",
    "bcextrapolate",
    "inlet",
    "outlet",
    "farfield",
)


def is_opening_bc_type(cgns_bc_type: str | None) -> bool:
    """True when a CGNS BCType marks an inlet / outlet / farfield patch."""
    if not cgns_bc_type:
        return False
    bt = cgns_bc_type.lower()
    return any(tok in bt for tok in _OPENING_BC_TOKENS)


def build_merged_region_mesh(
    case: CGNSCase,
    zone_indices: list[int],
    zone_topos: list[_ZoneTopo],
    foam_name: str,
    report: CouplingReport,
) -> Mesh:
    """Concatenate selected zones into one OpenFOAM region mesh."""
    indices = list(zone_indices)
    if not indices:
        raise ValueError(f"no zones for region {foam_name!r}")

    # Local face offsets within the merged region (index by global zi)
    face_offset_map: dict[int, int] = {}
    vtx_off = 0
    cell_off = 0
    face_off = 0
    n_fv = 0
    for zi in indices:
        face_offset_map[zi] = face_off
        zt = zone_topos[zi]
        vtx_off += zt.n_vertices
        cell_off += zt.n_cells
        face_off += zt.n_faces
        n_fv += int(zt.face_vertices.size)

    n_points = vtx_off
    n_cells = cell_off
    n_faces = face_off

    points = np.empty((n_points, 3), dtype=np.float64)
    owner = np.empty(n_faces, dtype=np.int64)
    neighbour = np.empty(n_faces, dtype=np.int64)
    flip = np.empty(n_faces, dtype=bool)
    face_offsets = np.empty(n_faces + 1, dtype=np.int64)
    face_offsets[0] = 0
    face_vertices = np.empty(n_fv, dtype=np.int64)

    vtx_cursor = 0
    cell_cursor = 0
    face_cursor = 0
    fv_cursor = 0
    # Per-source-zone cellZones (useful for MRF) plus the merged region name.
    cell_zones: list[CellZone] = []
    all_labels: list[np.ndarray] = []

    for zi in indices:
        zone = case.zones[zi]
        zt = zone_topos[zi]
        n_v = zt.n_vertices
        n_c = zt.n_cells
        n_f = zt.n_faces
        n_fvl = int(zt.face_vertices.size)

        points[vtx_cursor:vtx_cursor + n_v] = zone.coords
        face_vertices[fv_cursor:fv_cursor + n_fvl] = zt.face_vertices + vtx_cursor
        face_offsets[face_cursor + 1:face_cursor + n_f + 1] = (
            zt.face_offsets[1:] + fv_cursor
        )
        owner[face_cursor:face_cursor + n_f] = zt.owner + cell_cursor
        nb = zt.neighbour.copy()
        nb[nb >= 0] += cell_cursor
        neighbour[face_cursor:face_cursor + n_f] = nb
        flip[face_cursor:face_cursor + n_f] = zt.flip
        labels = np.arange(cell_cursor, cell_cursor + n_c, dtype=np.int64)
        all_labels.append(labels)
        cz_name = _sanitize_patch_name(zone.name)
        if cz_name != foam_name:
            cell_zones.append(CellZone(name=cz_name, cell_labels=labels))

        vtx_cursor += n_v
        cell_cursor += n_c
        face_cursor += n_f
        fv_cursor += n_fvl

    face_offset_list = [0] * len(case.zones)
    for zi, off in face_offset_map.items():
        face_offset_list[zi] = off

    plan = _region_patch_plan(
        indices, case, zone_topos, foam_name, report,
        face_offset=face_offset_list,
    )
    # Prefer per-source-zone cellZones (needed by MRF). Only add a region-wide
    # zone when this region is a single CGNS zone (typical solids).
    if not cell_zones:
        cell_zones = [
            CellZone(
                name=foam_name,
                cell_labels=(
                    np.concatenate(all_labels)
                    if all_labels
                    else np.empty(0, dtype=np.int64)
                ),
            )
        ]
    return _assemble_ordered_mesh(
        points,
        face_offsets,
        face_vertices,
        owner,
        neighbour,
        flip,
        plan,
        n_cells=n_cells,
        cell_zones=cell_zones,
    )


def build_region_meshes(
    case: CGNSCase,
    report: CouplingReport,
    zone_topos: list[_ZoneTopo],
) -> dict[str, Mesh]:
    """Build one OpenFOAM :class:`Mesh` per OpenFOAM region (JSON foam name)."""
    groups: dict[str, list[int]] = defaultdict(list)
    for i, reg in enumerate(report.regions):
        groups[reg.foam_name].append(i)

    out: dict[str, Mesh] = {}
    for foam_name in list(report.fluid_regions) + list(report.solid_regions):
        indices = groups.get(foam_name)
        if not indices:
            continue
        out[foam_name] = build_merged_region_mesh(
            case, indices, zone_topos, foam_name, report,
        )
    # Any foam names only present on zone infos
    for foam_name, indices in groups.items():
        if foam_name not in out:
            out[foam_name] = build_merged_region_mesh(
                case, indices, zone_topos, foam_name, report,
            )
    return out


def _allrun_direct(*, n_procs: int = 8) -> str:
    return f"""#!/bin/sh
set -e
cd "${{0%/*}}" || exit
. ${{WM_PROJECT_DIR:?}}/bin/tools/RunFunctions
#------------------------------------------------------------------------------
# cgns2foam --cht-direct: regions already split; decompose + parallel run
runApplication -o -s decomposePar decomposePar -allRegions -copyZero -force
runParallel -o -np {n_procs} $(getApplication)
runApplication -o -s reconstructParMesh reconstructParMesh -allRegions -constant
runApplication -o -s reconstructPar reconstructPar -allRegions
#------------------------------------------------------------------------------
"""


def write_cht_direct_case(
    out_dir: str | Path,
    case: CGNSCase,
    report: CouplingReport,
    zone_topos: list[_ZoneTopo],
    *,
    source_path: str,
    regions_config: RegionsConfig | None = None,
    end_time: int = 500,
    write_interval: int = 50,
    purge_write: int = 0,
    n_procs: int = 8,
    gravity: list[float] | tuple[float, float, float] | None = None,
    initial_t: float = 300.0,
    initial_p: float = 101325.0,
) -> dict[str, Any]:
    """Write a ready-to-run multi-region chtMultiRegionSimpleFoam case."""
    out = Path(out_dir)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    region_meshes = build_region_meshes(case, report, zone_topos)
    fluid = list(report.fluid_regions)
    solid = list(report.solid_regions)
    region_type = {}
    for r in report.regions:
        region_type[r.foam_name] = r.region_type

    warnings: list[str] = []
    mrf_specs = list(regions_config.mrf_regions) if regions_config else []
    mrf_summary: list[dict[str, Any]] = []
    # Aggregate heat sources per region (total watts); several JSON keys may
    # point at the same region — sum instead of silently dropping them.
    heat_by_region: dict[str, float] = defaultdict(float)
    if regions_config:
        for hs in regions_config.heat_sources:
            heat_by_region[hs.region_name] += hs.power
    ext_conv = regions_config.external_convection if regions_config else None

    if gravity is None and regions_config is not None:
        gravity = regions_config.gravity

    _write_text(out / "constant" / "regionProperties", _region_properties(fluid, solid))
    _write_text(out / "constant" / "g", _gravity(gravity))
    _write_text(
        out / "system" / "controlDict",
        _control_dict_cht(
            end_time=end_time,
            write_interval=write_interval,
            purge_write=purge_write,
        ),
    )
    _write_text(out / "system" / "fvSchemes", _fv_schemes_fluid())
    _write_text(out / "system" / "fvSolution", _fv_solution_fluid())
    _write_text(out / "system" / "decomposeParDict", _decompose_par_dict(n_procs))

    opening_patches_all: dict[str, list[str]] = {}
    convection_applied: dict[str, dict[str, list[float]]] = {}
    materials_applied: dict[str, dict[str, float]] = {}

    for foam_name, mesh in region_meshes.items():
        rtype = region_type.get(foam_name)
        if rtype is None:
            rtype = "solid" if foam_name in solid else "fluid"
            warnings.append(
                f"region {foam_name!r} not in region scan; assumed {rtype}"
            )
        loc = f"constant/{foam_name}/polyMesh"
        opts = WriteOptions(
            mesh_format="binary",
            mesh_location=loc,
            full_neighbour=False,
            ansa_headers=False,
        )
        write_poly_mesh(
            str(out / "constant" / foam_name / "polyMesh"),
            mesh,
            source_path,
            options=opts,
        )

        cdir = out / "constant" / foam_name
        sdir = out / "system" / foam_name
        zdir = out / "0" / foam_name
        patch_names = [p.name for p in mesh.patches]
        patch_types = {p.name: p.bc_type for p in mesh.patches}

        # Opening patches: "open*" by name, or a CGNS inlet/outlet/farfield BC.
        openings: set[str] = set()
        bc_typed_openings: list[str] = []
        for p in mesh.patches:
            if p.name == "open" or p.name.startswith("open"):
                openings.add(p.name)
            elif is_opening_bc_type(p.cgns_bc_type):
                openings.add(p.name)
                bc_typed_openings.append(f"{p.name} ({p.cgns_bc_type})")
        if openings:
            opening_patches_all[foam_name] = sorted(openings)
        if bc_typed_openings:
            warnings.append(
                f"region {foam_name!r}: CGNS inlet/outlet BC(s) mapped to "
                f"total-pressure openings: {', '.join(bc_typed_openings)}; "
                "values default to p0/T0 - edit 0/ as needed"
            )

        # External convection on outer walls (regex match from JSON).
        conv: dict[str, tuple[float, float]] = {}
        if ext_conv is not None:
            for pname in patch_names:
                if ext_conv.matches(pname):
                    conv[pname] = (ext_conv.ta, ext_conv.h)
            if conv:
                convection_applied[foam_name] = {
                    p: [ta, h] for p, (ta, h) in sorted(conv.items())
                }

        mat = regions_config.material_for(foam_name) if regions_config else None
        if mat is not None:
            applied = {
                k: v
                for k, v in (
                    ("rho", mat.rho),
                    ("Cp", mat.cp),
                    ("kappa", mat.kappa),
                    ("mu", mat.mu),
                    ("Pr", mat.pr),
                    ("molWeight", mat.mol_weight),
                )
                if v is not None
            }
            if applied:
                materials_applied[foam_name] = applied

        if rtype == "fluid":
            _write_text(cdir / "thermophysicalProperties", _thermophysical_fluid(mat))
            _write_text(cdir / "turbulenceProperties", _turbulence_laminar())
            _write_text(cdir / "radiationProperties", _radiation_none())
            if mrf_specs and foam_name == MERGED_FLUID_REGION:
                valid_cz = {cz.name for cz in mesh.cell_zones}
                for spec in mrf_specs:
                    if spec.foam_cell_zone not in valid_cz:
                        raise ValueError(
                            f"MRF cellZone {spec.foam_cell_zone!r} is not part of "
                            f"the {MERGED_FLUID_REGION!r} region mesh "
                            f"(cellZones: {sorted(valid_cz)}); list the CGNS zone "
                            "under fluid_regions in the sidecar JSON"
                        )
                mrf_entries = resolve_mrf_entries(mrf_specs, case, patch_names)
                _write_text(
                    cdir / "MRFProperties",
                    mrf_properties(
                        mrf_entries,
                        location=f"constant/{foam_name}",
                    ),
                )
                mrf_summary = mrf_entries
            _write_text(sdir / "fvSchemes", _fv_schemes_fluid())
            _write_text(sdir / "fvSolution", _fv_solution_fluid())
            _write_text(sdir / "fvOptions", _fv_options_fluid())
            _write_text(sdir / "decomposeParDict", _decompose_par_dict(n_procs, location="system"))
            _write_text(
                zdir / "T",
                _field_T(
                    "fluid",
                    patch_names,
                    initial_t,
                    patch_types=patch_types,
                    opening_patches=openings,
                    convection_patches=conv,
                ),
            )
            _write_text(
                zdir / "U",
                _field_U(
                    patch_names,
                    patch_types=patch_types,
                    moving_wall_patches=[
                        p for p in patch_names if "impeller" in p.lower()
                    ] if mrf_specs else None,
                    opening_patches=openings,
                ),
            )
            _write_text(
                zdir / "p",
                _field_p(
                    patch_names,
                    initial_p,
                    patch_types=patch_types,
                    opening_patches=openings,
                ),
            )
            _write_text(
                zdir / "p_rgh",
                _field_p_rgh(
                    patch_names,
                    initial_p,
                    patch_types=patch_types,
                    opening_patches=openings,
                ),
            )
        else:
            _write_text(cdir / "thermophysicalProperties", _thermophysical_solid(mat))
            _write_text(cdir / "radiationProperties", _radiation_none())
            _write_text(sdir / "fvSchemes", _fv_schemes_solid())
            _write_text(sdir / "fvSolution", _fv_solution_solid())
            _write_text(sdir / "decomposeParDict", _decompose_par_dict(n_procs, location="system"))
            # Volumetric heat source if this solid region has one
            power = heat_by_region.get(foam_name)
            if power:
                _write_text(sdir / "fvOptions", _fv_options_solid_heat(power))
            _write_text(
                zdir / "T",
                _field_T(
                    "solid",
                    patch_names,
                    initial_t,
                    patch_types=patch_types,
                    opening_patches=openings,
                    convection_patches=conv,
                ),
            )
            _write_text(
                zdir / "p",
                _field_p(
                    patch_names,
                    initial_p,
                    patch_types=patch_types,
                    opening_patches=openings,
                ),
            )

    summary = report.to_dict()
    summary["mode"] = "cht-direct"
    summary["solver"] = "chtMultiRegionSimpleFoam"
    if mrf_summary:
        summary["mrf"] = mrf_summary
    summary["settings"] = {
        "endTime": end_time,
        "writeInterval": write_interval,
        "purgeWrite": purge_write,
        "nProcs": n_procs,
        "gravity": list(gravity) if gravity is not None else [0.0, 0.0, -9.81],
        "initialT": initial_t,
        "initialP": initial_p,
    }
    if heat_by_region:
        summary["heat_sources"] = {k: v for k, v in sorted(heat_by_region.items())}
    if materials_applied:
        summary["materials_applied"] = materials_applied
    if convection_applied:
        summary["external_convection"] = convection_applied
    if opening_patches_all:
        summary["opening_patches"] = opening_patches_all
    if warnings:
        summary["warnings"] = warnings
    summary["regions_written"] = {
        name: {
            "n_cells": mesh.n_cells,
            "n_faces": int(mesh.owner.size),
            "n_internal_faces": mesh.n_internal_faces,
            "patches": [
                {
                    "name": p.name,
                    "type": p.bc_type,
                    "nFaces": p.n_faces,
                    "sampleRegion": p.sample_region,
                    "samplePatch": p.sample_patch,
                    "neighbourPatch": p.neighbour_patch,
                }
                for p in mesh.patches
            ],
        }
        for name, mesh in region_meshes.items()
    }
    text = json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    _write_text(out / "coupling_scan.json", text)
    _write_text(out / "setup_report.json", text)
    _write_text(out / "Allrun", _allrun_direct(n_procs=n_procs))
    _write_text(out / "Allclean", _allclean())
    for name in ("Allrun", "Allclean"):
        path = out / name
        try:
            path.chmod(path.stat().st_mode | 0o111)
        except OSError:
            pass
    return summary


def convert_cht_direct(
    cgns_path: str,
    out_dir: str,
    *,
    verbose: bool = True,
    solid_patterns: list[str] | None = None,
    fluid_patterns: list[str] | None = None,
    point_tol: float = 1e-4,
    regions_config: RegionsConfig | None = None,
) -> CouplingReport:
    """CGNS → multi-region chtMultiRegionSimpleFoam case in one pass."""
    t0 = time.perf_counter()
    case = read_cgns(cgns_path)
    if verbose:
        print(f"[cgns2foam] loaded {cgns_path} (cht-direct)")
        print(f"            nZones={len(case.zones)}")

    if regions_config is None:
        regions_config = load_sidecar_regions(
            cgns_path,
            [z.name for z in case.zones],
            required=True,
        )
    if verbose and regions_config is not None:
        print(f"[cgns2foam] regions from {regions_config.path}")

    t1 = time.perf_counter()
    zone_topos = prepare_zone_topos(case, point_tol=point_tol, trim=True)
    report = scan_couplings(
        case,
        source=os.path.abspath(cgns_path),
        solid_patterns=solid_patterns,
        fluid_patterns=fluid_patterns,
        regions_config=regions_config,
        trim=False,
        zone_topos=zone_topos,
    )
    if verbose:
        print(f"[cgns2foam] topology+scan [{time.perf_counter() - t1:.2f}s]")
        print(format_coupling_summary(report))
        cross_ff = [
            c for c in report.couplings
            if c.kind.value == "fluid_fluid"
            and c.master_region
            and c.slave_region
            and c.master_region != c.slave_region
        ]
        if cross_ff:
            print(
                "[cgns2foam] warning: "
                f"{len(cross_ff)} fluid-fluid pair(s) span different OpenFOAM "
                "regions; cyclicAMI needs one polyMesh - put those cellZones "
                "under one fluid region in the sidecar JSON."
            )

    if verbose and regions_config is not None and regions_config.mrf_regions:
        print(
            f"[cgns2foam] MRF: {len(regions_config.mrf_regions)} rotating "
            f"cellZone(s) -> constant/air/MRFProperties"
        )
        for m in regions_config.mrf_regions:
            org = "centroid" if m.origin is None else list(m.origin)
            print(
                f"            - {m.foam_cell_zone}  omega={m.omega}  "
                f"axis={list(m.axis)}  origin={org}"
            )

    def _cfg(name: str, default):
        return getattr(regions_config, name) if regions_config is not None else None

    t2 = time.perf_counter()
    summary = write_cht_direct_case(
        out_dir,
        case,
        report,
        zone_topos,
        source_path=os.path.abspath(cgns_path),
        regions_config=regions_config,
        end_time=_cfg("end_time", None) or 500,
        write_interval=_cfg("write_interval", None) or 50,
        purge_write=_cfg("purge_write", None) or 0,
        n_procs=_cfg("n_procs", None) or 8,
        gravity=_cfg("gravity", None),
        initial_t=_cfg("initial_t", None) or 300.0,
        initial_p=_cfg("initial_p", None) or 101325.0,
    )
    if verbose:
        for w in summary.get("warnings", []):
            print(f"[cgns2foam] warning: {w}")
        print(
            f"[cgns2foam] cht-direct case written to {out_dir} "
            f"({len(summary.get('fluid_regions', []))} fluid, "
            f"{len(summary.get('solid_regions', []))} solid, "
            f"{summary.get('n_couplings', 0)} couplings"
            f"{', mrf=' + str(len(summary.get('mrf', []))) if summary.get('mrf') else ''}) "
            f"[{time.perf_counter() - t2:.2f}s write, "
            f"{time.perf_counter() - t0:.2f}s total]"
        )
        print(f"            next: cd {out_dir} && ./Allrun  (OpenFOAM v2412)")
    return report
