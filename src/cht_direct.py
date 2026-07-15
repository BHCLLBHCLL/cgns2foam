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
    _field_T,
    _field_U,
    _field_p,
    _field_p_rgh,
    _fv_schemes_fluid,
    _fv_schemes_solid,
    _fv_solution_fluid,
    _fv_solution_solid,
    _gravity,
    _region_properties,
    _thermophysical_fluid,
    _thermophysical_solid,
    _turbulence_laminar,
    _write_text,
)
from .couplings import (
    CouplingMethod,
    CouplingReport,
    format_coupling_summary,
    prepare_zone_topos,
    scan_couplings,
)
from .reader import CGNSCase, read_cgns
from .regions_config import RegionsConfig, load_sidecar_regions
from .topology import (
    CellZone,
    Mesh,
    _ZoneTopo,
    _assemble_ordered_mesh,
    _bc_type_to_foam,
    _sanitize_patch_name,
)
from .writer import WriteOptions, write_poly_mesh


def coupling_patch_name(local: str, remote: str) -> str:
    """mappedWall inter-region patch name."""
    return f"{local}_to_{remote}"


def ami_patch_name(local_zone: str, remote_zone: str) -> str:
    """cyclicAMI patch name (no ``_to_`` so field stubs treat it as AMI)."""
    a = _sanitize_patch_name(local_zone).replace(".", "_")
    b = _sanitize_patch_name(remote_zone).replace(".", "_")
    # shorten to last path-like token if very long
    if len(a) > 40:
        a = a.split("_")[-1] or a[-40:]
    if len(b) > 40:
        b = b.split("_")[-1] or b[-40:]
    return f"ami_{a}__{b}"


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

        # fluid–fluid across distinct OpenFOAM regions: cyclicAMI needs one
        # polyMesh. Emit mappedWall and rely on JSON merging coupled fluids
        # into one region when true AMI flow coupling is required.
        if c.method == CouplingMethod.CYCLIC_AMI:
            pname = coupling_patch_name(foam_name, remote)
            foam_type = "mappedWall"
            extras = {
                "sample_mode": "nearestPatchFace",
                "sample_region": remote,
                "sample_patch": coupling_patch_name(remote, foam_name),
            }
        elif c.method == CouplingMethod.MAPPED_WALL:
            pname = coupling_patch_name(foam_name, remote)
            foam_type = "mappedWall"
            extras = {
                "sample_mode": "nearestPatchFace",
                "sample_region": remote,
                "sample_patch": coupling_patch_name(remote, foam_name),
            }
        else:
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
        patches.append((pname, foam_type, global_ids(zi, ids), extras))
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
                (pname, _bc_type_to_foam(bc.bc_type), global_ids(zi, ids), {})
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
            patches.append((pname, "wall", global_ids(zi, remaining), {}))
            assigned[zi][remaining] = True

    return patches


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
    cell_labels_all: list[np.ndarray] = []

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
        cell_labels_all.append(
            np.arange(cell_cursor, cell_cursor + n_c, dtype=np.int64)
        )

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
    return _assemble_ordered_mesh(
        points,
        face_offsets,
        face_vertices,
        owner,
        neighbour,
        flip,
        plan,
        n_cells=n_cells,
        cell_zones=[
            CellZone(
                name=foam_name,
                cell_labels=np.concatenate(cell_labels_all),
            )
        ],
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


def _allrun_direct() -> str:
    return """#!/bin/sh
set -e
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions
#------------------------------------------------------------------------------
# cgns2foam --cht-direct: regions already split; just run the solver
runApplication $(getApplication)
#------------------------------------------------------------------------------
"""


def write_cht_direct_case(
    out_dir: str | Path,
    case: CGNSCase,
    report: CouplingReport,
    zone_topos: list[_ZoneTopo],
    *,
    source_path: str,
    end_time: int = 500,
    write_interval: int = 50,
    gravity: list[float] | None = None,
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

    _write_text(out / "constant" / "regionProperties", _region_properties(fluid, solid))
    _write_text(out / "constant" / "g", _gravity(gravity))
    _write_text(
        out / "system" / "controlDict",
        _control_dict_cht(end_time=end_time, write_interval=write_interval),
    )
    _write_text(out / "system" / "fvSchemes", _fv_schemes_fluid())
    _write_text(out / "system" / "fvSolution", _fv_solution_fluid())

    for foam_name, mesh in region_meshes.items():
        rtype = region_type.get(foam_name, "fluid")
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

        if rtype == "fluid":
            _write_text(cdir / "thermophysicalProperties", _thermophysical_fluid())
            _write_text(cdir / "turbulenceProperties", _turbulence_laminar())
            _write_text(sdir / "fvSchemes", _fv_schemes_fluid())
            _write_text(sdir / "fvSolution", _fv_solution_fluid())
            _write_text(zdir / "T", _field_T("fluid", patch_names, patch_types=patch_types))
            _write_text(zdir / "U", _field_U(patch_names, patch_types=patch_types))
            _write_text(zdir / "p", _field_p(patch_names, patch_types=patch_types))
            _write_text(zdir / "p_rgh", _field_p_rgh(patch_names, patch_types=patch_types))
        else:
            _write_text(cdir / "thermophysicalProperties", _thermophysical_solid())
            _write_text(sdir / "fvSchemes", _fv_schemes_solid())
            _write_text(sdir / "fvSolution", _fv_solution_solid())
            _write_text(zdir / "T", _field_T("solid", patch_names, patch_types=patch_types))
            _write_text(zdir / "p", _field_p(patch_names, patch_types=patch_types))

    summary = report.to_dict()
    summary["mode"] = "cht-direct"
    summary["solver"] = "chtMultiRegionSimpleFoam"
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
    _write_text(out / "Allrun", _allrun_direct())
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
                f"{len(cross_ff)} fluid–fluid pair(s) span different OpenFOAM "
                "regions; cyclicAMI needs one polyMesh — merge those cellZones "
                "under one fluid region in the sidecar JSON."
            )

    t2 = time.perf_counter()
    summary = write_cht_direct_case(
        out_dir,
        case,
        report,
        zone_topos,
        source_path=os.path.abspath(cgns_path),
    )
    if verbose:
        print(
            f"[cgns2foam] cht-direct case written to {out_dir} "
            f"({len(summary.get('fluid_regions', []))} fluid, "
            f"{len(summary.get('solid_regions', []))} solid, "
            f"{summary.get('n_couplings', 0)} couplings) "
            f"[{time.perf_counter() - t2:.2f}s write, "
            f"{time.perf_counter() - t0:.2f}s total]"
        )
        print(f"            next: cd {out_dir} && ./Allrun  (OpenFOAM v2412)")
    return report
