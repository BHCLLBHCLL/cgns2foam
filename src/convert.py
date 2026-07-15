"""High-level orchestration: glue reader → topology → writer together."""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

from .cht_case import write_cht_case
from .couplings import CouplingReport, format_coupling_summary, scan_couplings
from .reader import read_cgns
from .topology import Mesh, build_mesh
from .writer import WriteOptions, write_case


def convert_file(
    cgns_path: str,
    out_dir: str,
    *,
    verbose: bool = True,
    write_options: WriteOptions | None = None,
    cht: bool = False,
    solid_patterns: list[str] | None = None,
    fluid_patterns: list[str] | None = None,
) -> Mesh:
    """Convert ``cgns_path`` to an OpenFOAM case rooted at ``out_dir``.

    When *cht* is True, also write ``chtMultiRegionSimpleFoam`` scaffolding
    (regionProperties, per-region thermo / 0.orig, Allrun.pre) based on an
    automatic coupling scan.

    Returns the intermediate :class:`~src.topology.Mesh`.
    """
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
    if cht:
        t_scan = time.perf_counter()
        coupling_report = scan_couplings(
            case,
            source=os.path.abspath(cgns_path),
            solid_patterns=solid_patterns,
            fluid_patterns=fluid_patterns,
        )
        if verbose:
            print(f"[cgns2foam] coupling scan "
                  f"({time.perf_counter() - t_scan:.2f}s):")
            print(format_coupling_summary(coupling_report))

    t1 = time.perf_counter()
    mesh = build_mesh(case)
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
) -> CouplingReport:
    """Scan CGNS structure for regions and coupling pairs (no mesh write)."""
    t0 = time.perf_counter()
    case = read_cgns(cgns_path)
    report = scan_couplings(
        case,
        source=os.path.abspath(cgns_path),
        solid_patterns=solid_patterns,
        fluid_patterns=fluid_patterns,
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
