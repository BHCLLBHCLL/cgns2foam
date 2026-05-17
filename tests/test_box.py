"""Small unit-style test for the box case.

Verifies that converting ``cases/box_ansa/box_ansa_orig_fix.cgns`` yields
an OpenFOAM mesh with the same gross statistics as the ANSA-produced
reference.

Run with::

    python3 -m unittest tests.test_box

(no pytest dependency required).
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

from cgns2foam.convert import convert_file        # noqa: E402
from validate import (                            # noqa: E402
    read_boundary,
    read_faces_any,
    read_owner,
    read_neighbour,
    read_points,
)


CGNS_FILE = ROOT / "cases" / "box_ansa" / "box_ansa_orig_fix.cgns"
ZIP_REF = ROOT / "cases" / "box_ansa" / "box_ansa_orig_fix.zip"


def _stats(case_dir: Path) -> dict:
    poly = case_dir / "constant" / "polyMesh"
    pts = read_points(str(poly / "points"))
    own = read_owner(str(poly / "owner"))
    nei = read_neighbour(str(poly / "neighbour"))
    ofs, conn = read_faces_any(str(poly / "faces"))
    bnd = read_boundary(str(poly / "boundary"))
    nei_int = nei[nei >= 0]
    return {
        "nPoints": int(pts.shape[0]),
        "nFaces": int(own.size),
        "nInternalFaces": int(nei_int.size),
        "nCells": int(max(own.max(), nei.max())) + 1,
        "nBoundaryFaces": int(own.size - nei_int.size),
        "patchTotalNFaces": sum(p["nFaces"] for p in bnd),
        "patchNames": [p["name"] for p in bnd],
    }


class TestBoxCase(unittest.TestCase):
    def test_convert_matches_reference(self):
        self.assertTrue(CGNS_FILE.is_file(), f"missing input {CGNS_FILE}")
        self.assertTrue(ZIP_REF.is_file(), f"missing reference {ZIP_REF}")
        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td) / "out"
            convert_file(str(CGNS_FILE), str(out_dir), verbose=False)
            ours = _stats(out_dir)

            ref_root = Path(td) / "ref"
            ref_root.mkdir()
            with zipfile.ZipFile(ZIP_REF) as zf:
                zf.extractall(ref_root)
            ref_case = None
            for root, _d, files in os.walk(ref_root):
                if "points" in files and Path(root).name == "polyMesh":
                    ref_case = Path(root).parent.parent
                    break
            self.assertIsNotNone(ref_case)
            ref = _stats(ref_case)

            for key in ("nPoints", "nFaces", "nInternalFaces", "nCells",
                        "nBoundaryFaces", "patchTotalNFaces"):
                self.assertEqual(
                    ours[key], ref[key],
                    f"mismatch on {key}: ours={ours[key]} ref={ref[key]}",
                )
            self.assertEqual(ours["patchNames"], ref["patchNames"])


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
