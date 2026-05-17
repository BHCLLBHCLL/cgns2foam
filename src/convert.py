"""High-level orchestration: glue reader → topology → writer together."""

from __future__ import annotations

import os
import time

from .reader import read_cgns
from .topology import Mesh, build_mesh
from .writer import write_case


def convert_file(cgns_path: str, out_dir: str, *, verbose: bool = True) -> Mesh:
    """Convert ``cgns_path`` to an OpenFOAM case rooted at ``out_dir``.

    Returns the intermediate :class:`~src.topology.Mesh` so that
    callers (tests, CLI) can introspect the result.
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
    write_case(out_dir, mesh, source_path=os.path.abspath(cgns_path))
    if verbose:
        print(f"[cgns2foam] case written to {out_dir} "
              f"[{time.perf_counter() - t2:.2f}s]")
    return mesh
