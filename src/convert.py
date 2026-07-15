"""High-level orchestration: glue reader → topology → writer together."""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

from .cht_case import write_cht_case
from .cht_direct import convert_cht_direct
from .couplings import CouplingReport, format_coupling_summary, scan_couplings
from .reader import read_cgns
from .regions_config import load_sidecar_regions
from .topology import Mesh, build_mesh
from .writer import WriteOptions, write_case


def convert_file(
    cgns_path: str,
    out_dir: str,
    *,
    verbose: bool = True,
    write_options: WriteOptions | None = None,
    cht: bool = False,
    cht_direct: bool = False,
    solid_patterns: list[str] | None = None,
    fluid_patterns: list[str] | None = None,
) -> Mesh | CouplingReport:
    """Convert ``cgns_path`` to an OpenFOAM case rooted at ``out_dir``.

    Modes:

    * default – mono-block polyMesh
    * ``cht=True`` – mono polyMesh + CHT scaffolding (``Allrun.pre`` runs
      ``splitMeshRegions``); requires sidecar ``<cgns>.json`` regions
    * ``cht_direct=True`` – one-step multi-region
      ``chtMultiRegionSimpleFoam`` case (no mono mesh / no split);
      requires sidecar ``<cgns>.json``

    Returns :class:`~src.topology.Mesh` for mono/cht modes, or
    :class:`~src.couplings.CouplingReport` for ``cht_direct``.
    """
    if cht_direct:
        return convert_cht_direct(
            cgns_path,
            out_dir,
            verbose=verbose,
            solid_patterns=solid_patterns,
            fluid_patterns=fluid_patterns,
        )

    t0 = time.perf_counter()
    case = read_cgns(cgns_path)
    if verbose:
        print(f"[cgns2foam] loaded {cgns_path}")
        print(f"            CellDim={case.cell_dim} PhysDim={case.phys_dim} "
              f"nZones={len(case.zones)}")
        for z in case.zones:
            print(f"            - zone {z.name!r}: "
                  f"{z.n_vertices} vertices, {z.n_cells} cells, "
                  f"{len(z.bcs)} BCs")

    coupling_report: CouplingReport | None = None
    regions_config = None
    zone_cell_zone_map: dict[str, str] | None = None
    if cht:
        regions_config = load_sidecar_regions(
            cgns_path,
            [z.name for z in case.zones],
            required=True,
        )
        if verbose and regions_config is not None:
            print(f"[cgns2foam] regions from {regions_config.path}")
        zone_cell_zone_map = {
            z.name: regions_config.foam_name_for(z.name) or z.name
            for z in case.zones
            if regions_config.foam_name_for(z.name)
        }
        unmatched = [
            z.name for z in case.zones
            if regions_config.foam_name_for(z.name) is None
        ]
        if unmatched and verbose:
            print(
                f"[cgns2foam] warning: {len(unmatched)} CGNS zone(s) not listed "
                f"in regions JSON: "
                f"{unmatched[:5]}{'...' if len(unmatched) > 5 else ''}"
            )
        t_scan = time.perf_counter()
        coupling_report = scan_couplings(
            case,
            source=os.path.abspath(cgns_path),
            solid_patterns=solid_patterns,
            fluid_patterns=fluid_patterns,
            regions_config=regions_config,
        )
        if verbose:
            print(f"[cgns2foam] coupling scan "
                  f"({time.perf_counter() - t_scan:.2f}s):")
            print(format_coupling_summary(coupling_report))

    t1 = time.perf_counter()
    mesh = build_mesh(case, zone_cell_zone_map=zone_cell_zone_map)
    if verbose:
        print(f"[cgns2foam] mesh assembled: "
              f"{mesh.points.shape[0]} points, "
              f"{mesh.owner.size} faces "
              f"({mesh.n_internal_faces} internal), "
              f"{mesh.n_cells} cells, "
              f"{len(mesh.patches)} patches, "
              f"{len(mesh.cell_zones)} cellZones "
              f"[{t1 - t0:.2f}s read + {time.perf_counter() - t1:.2f}s build]")

    t2 = time.perf_counter()
    write_case(out_dir, mesh, source_path=os.path.abspath(cgns_path),
               options=write_options)
    if verbose:
        print(f"[cgns2foam] case written to {out_dir} "
              f"[{time.perf_counter() - t2:.2f}s]")

    if cht and coupling_report is not None:
        t3 = time.perf_counter()
        summary = write_cht_case(out_dir, mesh, coupling_report)
        if verbose:
            print(f"[cgns2foam] CHT scaffolding written "
                  f"({len(summary.get('fluid_regions', []))} fluid, "
                  f"{len(summary.get('solid_regions', []))} solid, "
                  f"{summary.get('n_couplings', 0)} couplings) "
                  f"[{time.perf_counter() - t3:.2f}s]")
            print(f"            next: cd {out_dir} && ./Allrun.pre  "
                  f"(inside OpenFOAM v2412)")
    return mesh


def scan_file(
    cgns_path: str,
    *,
    report_path: str | None = None,
    verbose: bool = True,
    solid_patterns: list[str] | None = None,
    fluid_patterns: list[str] | None = None,
    use_regions_json: bool = True,
) -> CouplingReport:
    """Scan CGNS structure for regions and coupling pairs (no mesh write).

    When a sidecar ``<cgns>.json`` exists and *use_regions_json* is True, zone
    types / OpenFOAM names are taken from that file.
    """
    t0 = time.perf_counter()
    case = read_cgns(cgns_path)
    regions_config = None
    if use_regions_json:
        regions_config = load_sidecar_regions(
            cgns_path,
            [z.name for z in case.zones],
            required=False,
        )
        if verbose and regions_config is not None:
            print(f"[cgns2foam] regions from {regions_config.path}")
    report = scan_couplings(
        case,
        source=os.path.abspath(cgns_path),
        solid_patterns=solid_patterns,
        fluid_patterns=fluid_patterns,
        regions_config=regions_config,
    )
    if verbose:
        print(f"[cgns2foam] coupling scan of {cgns_path} "
              f"[{time.perf_counter() - t0:.2f}s]")
        print(format_coupling_summary(report))
    if report_path:
        path = Path(report_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(report.to_dict(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        if verbose:
            print(f"[cgns2foam] report written to {path}")
    return report
