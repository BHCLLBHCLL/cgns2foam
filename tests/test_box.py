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

from src.convert import convert_file              # noqa: E402
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


def _binary_faces_layout(faces_path: Path) -> None:
    """Default export: binary faceCompactList with ANSA banner."""
    data = faces_path.read_bytes()
    assert b"ANSA_VERSION: 25.1.0" in data
    assert b"faceCompactList" in data
    assert b"format binary" in data
    ofs, conn = read_faces_any(str(faces_path))
    assert ofs.size > 1
    assert conn.size > 0


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

    def test_faces_binary_compact_layout(self):
        self.assertTrue(CGNS_FILE.is_file(), f"missing input {CGNS_FILE}")
        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td) / "out"
            convert_file(str(CGNS_FILE), str(out_dir), verbose=False)
            _binary_faces_layout(out_dir / "constant" / "polyMesh" / "faces")

    def test_case_file_ansa_headers(self):
        """system/constant/0 use ANSA banner, location \"\" and format binary."""
        self.assertTrue(CGNS_FILE.is_file(), f"missing input {CGNS_FILE}")
        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td) / "out"
            convert_file(str(CGNS_FILE), str(out_dir), verbose=False)
            for rel in (
                "system/controlDict",
                "system/fvSchemes",
                "constant/turbulenceProperties",
                "0/U",
            ):
                data = (out_dir / rel).read_bytes()
                self.assertIn(b"ANSA_VERSION: 25.1.0", data, rel)
                self.assertIn(b'location "";', data, rel)
                self.assertIn(b"format binary;", data, rel)
            ctrl = (out_dir / "system" / "controlDict").read_text(encoding="ascii")
            self.assertIn("writeCompression\tuncompressed;", ctrl)
            self.assertIn("application UserSolver;", ctrl)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
