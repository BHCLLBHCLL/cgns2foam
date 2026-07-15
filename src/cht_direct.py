"""One-step CGNS → multi-region chtMultiRegionSimpleFoam case writer.

Each CGNS zone becomes ``constant/<region>/polyMesh`` with coupling BCs
rewritten as ``mappedWall`` (``local_to_remote``).  No mono-block merge and
no OpenFOAM ``splitMeshRegions`` step is required.
"""

from __future__ import annotations

import json
import os
import shutil
import time
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
    CouplingReport,
    format_coupling_summary,
    prepare_zone_topos,
    scan_couplings,
)
from .reader import CGNSCase, CGNSZone, read_cgns
from .topology import (
    Mesh,
    _ZoneTopo,
    _bc_type_to_foam,
    _sanitize_patch_name,
    build_zone_mesh,
)
from .writer import WriteOptions, write_poly_mesh


def coupling_patch_name(local: str, remote: str) -> str:
    return f"{local}_to_{remote}"


def _region_patch_plan(
    zone: CGNSZone,
    zt: _ZoneTopo,
    foam_name: str,
    report: CouplingReport,
    *,
    default_exterior_name: str = "default_exterior",
) -> list[tuple]:
    """Build patch plan for one zone (coupling mappedWall first, then BCs)."""
    foam_by_zone = {r.zone_name: r.foam_name for r in report.regions}
    boundary_mask = zt.neighbour < 0
    assigned = np.zeros(zt.n_faces, dtype=bool)
    patches: list[tuple] = []
    used: set[str] = set()

    for c in report.couplings:
        if c.master_zone == zone.name:
            local_bc, remote_zone = c.master_bc, c.slave_zone
        elif c.slave_zone == zone.name:
            local_bc, remote_zone = c.slave_bc, c.master_zone
        else:
            continue
        remote = foam_by_zone.get(remote_zone, _sanitize_patch_name(remote_zone))
        ids = zt.bc_face_lists.get(local_bc)
        if ids is None or ids.size == 0:
            continue
        keep = boundary_mask[ids] & ~assigned[ids]
        ids = ids[keep]
        if ids.size == 0:
            continue
        pname = coupling_patch_name(foam_name, remote)
        if pname in used:
            for i, entry in enumerate(patches):
                if entry[0] == pname:
                    patches[i] = (
                        pname,
                        entry[1],
                        np.unique(np.concatenate([entry[2], ids])),
                        entry[3],
                    )
                    assigned[ids] = True
                    break
            continue
        used.add(pname)
        extras = {
            "sample_mode": "nearestPatchFace",
            "sample_region": remote,
            "sample_patch": coupling_patch_name(remote, foam_name),
        }
        patches.append((pname, "mappedWall", ids.copy(), extras))
        assigned[ids] = True

    for bc in zone.bcs:
        if bc.name not in zt.bc_face_lists:
            continue
        ids = zt.bc_face_lists[bc.name]
        keep = boundary_mask[ids] & ~assigned[ids]
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
        patches.append((pname, _bc_type_to_foam(bc.bc_type), ids.copy(), {}))
        assigned[ids] = True

    remaining = np.where(boundary_mask & ~assigned)[0]
    if remaining.size:
        pname = default_exterior_name
        i = 1
        while pname in used:
            pname = f"{default_exterior_name}_{i}"
            i += 1
        patches.append((pname, "wall", remaining.copy(), {}))

    return patches


def build_region_meshes(
    case: CGNSCase,
    report: CouplingReport,
    zone_topos: list[_ZoneTopo],
) -> dict[str, Mesh]:
    """Build one OpenFOAM :class:`Mesh` per CGNS zone / CHT region."""
    out: dict[str, Mesh] = {}
    for zone, zt, reg in zip(case.zones, zone_topos, report.regions):
        foam = reg.foam_name
        plan = _region_patch_plan(zone, zt, foam, report)
        out[foam] = build_zone_mesh(zone, zt, plan, cell_zone_name=foam)
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
    region_type = {r.foam_name: r.region_type for r in report.regions}

    _write_text(out / "constant" / "regionProperties", _region_properties(fluid, solid))
    _write_text(out / "constant" / "g", _gravity(gravity))
    _write_text(
        out / "system" / "controlDict",
        _control_dict_cht(end_time=end_time, write_interval=write_interval),
    )
    _write_text(out / "system" / "fvSchemes", _fv_schemes_fluid())
    _write_text(out / "system" / "fvSolution", _fv_solution_fluid())

    for foam_name, mesh in region_meshes.items():
        rtype = region_type[foam_name]
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

        if rtype == "fluid":
            _write_text(cdir / "thermophysicalProperties", _thermophysical_fluid())
            _write_text(cdir / "turbulenceProperties", _turbulence_laminar())
            _write_text(sdir / "fvSchemes", _fv_schemes_fluid())
            _write_text(sdir / "fvSolution", _fv_solution_fluid())
            _write_text(zdir / "T", _field_T("fluid", patch_names))
            _write_text(zdir / "U", _field_U(patch_names))
            _write_text(zdir / "p", _field_p(patch_names))
            _write_text(zdir / "p_rgh", _field_p_rgh(patch_names))
        else:
            _write_text(cdir / "thermophysicalProperties", _thermophysical_solid())
            _write_text(sdir / "fvSchemes", _fv_schemes_solid())
            _write_text(sdir / "fvSolution", _fv_solution_solid())
            _write_text(zdir / "T", _field_T("solid", patch_names))
            _write_text(zdir / "p", _field_p(patch_names))

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
) -> CouplingReport:
    """CGNS → multi-region chtMultiRegionSimpleFoam case in one pass."""
    t0 = time.perf_counter()
    case = read_cgns(cgns_path)
    if verbose:
        print(f"[cgns2foam] loaded {cgns_path} (cht-direct)")
        print(f"            nZones={len(case.zones)}")

    t1 = time.perf_counter()
    zone_topos = prepare_zone_topos(case, point_tol=point_tol, trim=True)
    report = scan_couplings(
        case,
        source=os.path.abspath(cgns_path),
        solid_patterns=solid_patterns,
        fluid_patterns=fluid_patterns,
        trim=False,
        zone_topos=zone_topos,
    )
    if verbose:
        print(f"[cgns2foam] topology+scan [{time.perf_counter() - t1:.2f}s]")
        print(format_coupling_summary(report))

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
