"""End-to-end test driver.

For each test case under ``cases/``:

1. Run the cgns2foam converter to produce an OpenFOAM project.
2. Compare gross statistics with the ANSA-produced reference unzipped
   from the matching ``.zip`` file (number of points, faces, internal
   faces, cells, sum of patch face counts).
3. Optionally run OpenFOAM's ``checkMesh`` and capture the verdict.

Usage::

    python3 tests/run_all.py [--with-checkmesh] [--out-root /tmp/cgns2foam_out]
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from src.convert import convert_file  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(HERE))
from validate import (                       # noqa: E402
    read_boundary,
    read_faces_any,
    read_owner,
    read_neighbour,
    read_points,
)


CASES = [
    ("box_ansa",          "cases/box_ansa/box_ansa_orig_fix"),
    ("tr03",              "cases/tr03/tr03_orig_fix"),
    ("laptop_simplified", "cases/laptop_simplified/laptop_simplified_voxel_less_orig_fix"),
]


def _stats(case_dir: str) -> dict:
    poly = os.path.join(case_dir, "constant", "polyMesh")
    pts = read_points(os.path.join(poly, "points"))
    own = read_owner(os.path.join(poly, "owner"))
    nei = read_neighbour(os.path.join(poly, "neighbour"))
    ofs, conn = read_faces_any(os.path.join(poly, "faces"))
    bnd = read_boundary(os.path.join(poly, "boundary"))
    nei_int = nei[nei >= 0]
    n_cells = int(max(own.max(), nei.max())) + 1
    return {
        "nPoints": int(pts.shape[0]),
        "nFaces":  int(own.size),
        "nInternalFaces": int(nei_int.size),
        "nCells":  n_cells,
        "nPatches": len(bnd),
        "nBoundaryFaces": int(own.size - nei_int.size),
        "patchTotalNFaces": int(sum(p["nFaces"] for p in bnd)),
        "faceConnSum": int(conn.sum()),
    }


def _extract_reference(zip_path: str, target_root: str) -> str:
    """Extract the ANSA-produced zip and return the OpenFOAM case dir.

    The zips contain a single top-level directory whose name matches
    the basename of the .zip file; inside that directory sits the
    actual OpenFOAM case (the dir holding ``constant/polyMesh``).
    """
    name = Path(zip_path).stem
    out = Path(target_root) / "ref" / name
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(out)
    # Locate the case dir (the one containing constant/polyMesh).
    for root, _dirs, files in os.walk(out):
        if (Path(root) / "constant" / "polyMesh" / "points").is_file():
            return root
    raise FileNotFoundError(f"polyMesh not found under {out}")


_OPENFOAM_BASHRC_CANDIDATES = [
    # OpenFOAM v2412 / v2506 / … from openfoam.com (apt or tarball)
    "/usr/lib/openfoam/openfoam2412/etc/bashrc",
    "/usr/lib/openfoam/openfoam2506/etc/bashrc",
    "/opt/openfoam2412/etc/bashrc",
    "/opt/OpenFOAM-v2412/etc/bashrc",
    # Generic environment variable
    os.environ.get("FOAM_BASHRC", "") or "/dev/null",
]


def _foam_bashrc() -> str | None:
    for c in _OPENFOAM_BASHRC_CANDIDATES:
        if c and os.path.isfile(c):
            return c
    return None


def _run_check_mesh(case_dir: str, log_path: str) -> str:
    """Run ``checkMesh`` and return a one-line verdict."""
    fbash = _foam_bashrc()
    if fbash is None:
        return "checkMesh not available (no OpenFOAM v2412 bashrc found)"
    cmd = f"source {fbash} > /dev/null 2>&1 && cd {case_dir} && checkMesh -allTopology"
    try:
        r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True,
                           timeout=900)
    except subprocess.TimeoutExpired:
        return "checkMesh timed out"
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as fh:
        fh.write(r.stdout)
        fh.write("\n--- stderr ---\n")
        fh.write(r.stderr)
    text = r.stdout
    if "Mesh OK." in text:
        return "Mesh OK"
    # Count failed checks
    for line in reversed(text.splitlines()):
        if "Failed" in line and "mesh checks" in line:
            return line.strip()
    return "checkMesh finished (see log)"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--out-root", default="/tmp/cgns2foam_out")
    p.add_argument("--with-checkmesh", action="store_true")
    args = p.parse_args(argv)

    out_root = Path(args.out_root)
    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir(parents=True)

    overall_ok = True
    summary: list[str] = []
    summary.append("=" * 78)
    summary.append(f"cgns2foam end-to-end test run at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    summary.append("=" * 78)

    for name, rel in CASES:
        cgns = ROOT / f"{rel}.cgns"
        zip_ref = ROOT / f"{rel}.zip"
        out_case = out_root / "ours" / name
        log_path = out_root / "logs" / f"{name}.checkMesh.log"

        summary.append("")
        summary.append(f"### Case: {name}")
        summary.append(f"   input : {cgns}")
        summary.append(f"   output: {out_case}")

        t0 = time.perf_counter()
        convert_file(str(cgns), str(out_case), verbose=False)
        conv_s = time.perf_counter() - t0
        summary.append(f"   convert time: {conv_s:.2f}s")

        ours_stats = _stats(str(out_case))

        ref_dir = ""
        ref_stats = None
        if zip_ref.exists():
            ref_dir = _extract_reference(str(zip_ref), str(out_root))
            ref_stats = _stats(ref_dir)

        # Topology-invariant metrics – these must match the reference
        # exactly regardless of zone or patch ordering.
        topo_keys = ["nPoints", "nFaces", "nInternalFaces", "nCells",
                     "nBoundaryFaces", "patchTotalNFaces"]
        # Informational metrics: depends on point ordering / patch split
        info_keys = ["faceConnSum"]
        col = max(len(k) for k in topo_keys + info_keys) + 2
        summary.append(f"   {'metric'.ljust(col)} {'ours':>15} {'ref':>15} match")
        all_match = True
        for k in topo_keys:
            ov = ours_stats[k]
            rv = ref_stats[k] if ref_stats else None
            match = (rv is None) or (ov == rv)
            if not match:
                all_match = False
            summary.append(
                f"   {k.ljust(col)} {ov:>15} {rv if rv is not None else '-':>15} {'OK' if match else 'DIFF'}"
            )
        for k in info_keys:
            ov = ours_stats[k]
            rv = ref_stats[k] if ref_stats else None
            tag = "info" if (rv is None or ov == rv) else "info-diff"
            summary.append(
                f"   {k.ljust(col)} {ov:>15} {rv if rv is not None else '-':>15} {tag}"
            )

        summary.append(f"   nPatches: ours={ours_stats['nPatches']} "
                       f"ref={ref_stats['nPatches'] if ref_stats else '-'} "
                       f"(may differ when ANSA splits/coalesces patches)")
        if not all_match:
            overall_ok = False

        if args.with_checkmesh:
            verdict = _run_check_mesh(str(out_case), str(log_path))
            summary.append(f"   checkMesh ours: {verdict}  (log: {log_path})")
            if ref_dir:
                # The ANSA reference doesn't ship fvSchemes/fvSolution but
                # OpenFOAM v2412's checkMesh requires them, so we copy ours
                # into the reference tree before running checkMesh.
                ref_sys = Path(ref_dir) / "system"
                ref_sys.mkdir(exist_ok=True)
                for f in ("fvSchemes", "fvSolution"):
                    src = out_case / "system" / f
                    dst = ref_sys / f
                    if src.is_file() and not dst.is_file():
                        shutil.copy2(src, dst)
                ref_log = out_root / "logs" / f"{name}.checkMesh.ref.log"
                ref_verdict = _run_check_mesh(ref_dir, str(ref_log))
                summary.append(f"   checkMesh ref : {ref_verdict}  (log: {ref_log})")

    summary.append("")
    summary.append("=" * 78)
    summary.append(f"Overall result: {'PASS' if overall_ok else 'FAIL'}")
    summary.append("=" * 78)

    print("\n".join(summary))

    # also write summary
    (out_root / "summary.txt").write_text("\n".join(summary))
    return 0 if overall_ok else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
