"""Validation utility comparing a generated OpenFOAM case to a reference.

The comparison is *topological*, not byte-exact: face order within a patch
or vertex order around a face may legitimately differ between converters.
We check that:

* point sets are identical up to floating-point tolerance and a permutation,
* the boundary patch dictionary lists the same patches with the same
  number of faces (names may be sanitized differently, so we compare on
  ``n_faces`` multisets and on individual patches by name when possible),
* the number of cells, internal faces and total faces match,
* the cell -> face incidence (as multisets of (sorted face vertices))
  match per cell, after a point-permutation that maps the two point
  arrays.
"""

from __future__ import annotations

import os
import re
import struct
import sys
from collections import Counter
from typing import Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Low-level binary OpenFOAM polyMesh parser (just enough for validation)
# ---------------------------------------------------------------------------


def _read_file(path: str) -> bytes:
    with open(path, "rb") as fh:
        return fh.read()


def _skip_header(data: bytes) -> int:
    """Return the byte offset right after the FoamFile dictionary."""
    end = data.find(b"}\n", data.find(b"FoamFile"))
    if end < 0:
        raise ValueError("malformed FoamFile header")
    # Skip following comment lines until we hit the size token line.
    idx = end + 2
    # Skip comment lines starting with /*
    while True:
        nl = data.find(b"\n", idx)
        line = data[idx:nl]
        s = line.strip()
        if s.startswith(b"/*") or s == b"":
            idx = nl + 1
            continue
        break
    return idx


def _read_binary_label_list(data: bytes, start: int) -> tuple[np.ndarray, int]:
    """Parse ``<N>\\n(<N×int32>)`` from *start*. Returns (array, new_idx)."""
    paren = data.find(b"(", start)
    n = int(data[start:paren].strip())
    arr = np.frombuffer(data[paren + 1:paren + 1 + 4 * n], dtype="<i4")
    return arr.copy(), paren + 1 + 4 * n + 1   # skip closing ')'


def _read_binary_scalar_list(data: bytes, start: int) -> tuple[np.ndarray, int]:
    paren = data.find(b"(", start)
    n = int(data[start:paren].strip())
    arr = np.frombuffer(data[paren + 1:paren + 1 + 8 * n], dtype="<f8")
    return arr.copy(), paren + 1 + 8 * n + 1


def read_points(path: str) -> np.ndarray:
    data = _read_file(path)
    idx = _skip_header(data)
    paren = data.find(b"(", idx)
    n = int(data[idx:paren].strip())
    arr = np.frombuffer(data[paren + 1:paren + 1 + 8 * 3 * n], dtype="<f8")
    return arr.reshape(n, 3).copy()


def read_owner(path: str) -> np.ndarray:
    data = _read_file(path)
    idx = _skip_header(data)
    arr, _ = _read_binary_label_list(data, idx)
    return arr


def read_neighbour(path: str) -> np.ndarray:
    return read_owner(path)  # same encoding


def read_faces_compact(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (offsets, connectivity) for a compact-binary faces file.

    Format::

        <nOffsets>\n(<int32×nOffsets>)<nConn>\n(<int32×nConn>)
    """
    data = _read_file(path)
    idx = _skip_header(data)
    paren1 = data.find(b"(", idx)
    n_ofs = int(data[idx:paren1].strip())
    ofs = np.frombuffer(data[paren1 + 1:paren1 + 1 + 4 * n_ofs], dtype="<i4").copy()
    # right after closing ')'
    pos = paren1 + 1 + 4 * n_ofs + 1
    paren2 = data.find(b"(", pos)
    n_conn = int(data[pos:paren2].strip())
    conn = np.frombuffer(data[paren2 + 1:paren2 + 1 + 4 * n_conn], dtype="<i4").copy()
    return ofs, conn


def read_faces_classic(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Read the classic per-face binary face list ``N(k(vlist))``.

    Returns ``(offsets, connectivity)`` like the compact reader so the
    caller does not care which format the file uses.
    """
    data = _read_file(path)
    idx = _skip_header(data)
    paren = data.find(b"(", idx)
    n_faces = int(data[idx:paren].strip())
    pos = paren + 1
    offsets = [0]
    conn: list[int] = []
    for _ in range(n_faces):
        # ASCII size, then '(' then size×int32 then ')'
        # but for classic binary the size precedes each face as ascii
        # Try ascii decode for the size up to '('
        p2 = data.find(b"(", pos)
        sz = int(data[pos:p2].strip())
        verts = np.frombuffer(data[p2 + 1:p2 + 1 + 4 * sz], dtype="<i4")
        conn.extend(int(v) for v in verts)
        offsets.append(len(conn))
        pos = p2 + 1 + 4 * sz + 1   # skip past ')'
    return np.asarray(offsets, dtype=np.int64), np.asarray(conn, dtype=np.int64)


def read_faces_any(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Auto-detect the faces file flavour and return (offsets, conn)."""
    data = _read_file(path)
    cls = re.search(rb"class\s+(\w+)\s*;", data)
    if cls and cls.group(1) in (b"faceCompactList", b"faceCompactIOList"):
        return read_faces_compact(path)
    # ANSA also writes a compact format under ``faceList``. Heuristic:
    # parse the header count; if it equals nFaces+1 we have the compact form.
    idx = _skip_header(data)
    paren = data.find(b"(", idx)
    n_hdr = int(data[idx:paren].strip())
    # Read first int32 inside parens; for compact form it must be 0.
    first = struct.unpack("<i", data[paren + 1:paren + 5])[0]
    if first == 0:
        return read_faces_compact(path)
    return read_faces_classic(path)


def read_boundary(path: str) -> list[dict]:
    with open(path) as fh:
        text = fh.read()
    body = text[text.index("}\n", text.index("FoamFile")) + 2:]
    # Strip header comment block
    body = re.sub(r"/\*.*?\*/", "", body, flags=re.S)
    m = re.search(r"(\d+)\s*\(", body)
    n = int(m.group(1))
    start = m.end()
    patches = []
    pat = re.compile(
        r"(\w[\w.]*)\s*\{(.*?)\}", re.S
    )
    for pm in pat.finditer(body, start):
        name = pm.group(1)
        block = pm.group(2)
        d = {"name": name}
        for k in ("type", "nFaces", "startFace"):
            mm = re.search(rf"{k}\s+([^;]+);", block)
            if mm:
                v = mm.group(1).strip()
                d[k] = int(v) if k in ("nFaces", "startFace") else v
        patches.append(d)
        if len(patches) == n:
            break
    return patches


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------


def _face_keys(offsets: np.ndarray, conn: np.ndarray, point_perm: np.ndarray) -> list[tuple[int, ...]]:
    """Return a sorted-tuple key (in the *reference* point indexing) per face."""
    keys: list[tuple[int, ...]] = []
    for i in range(len(offsets) - 1):
        verts = conn[offsets[i]:offsets[i + 1]]
        keys.append(tuple(sorted(int(point_perm[v]) for v in verts)))
    return keys


def _point_permutation(pts_ours: np.ndarray, pts_ref: np.ndarray,
                       tol: float = 1e-9) -> np.ndarray:
    """Return ``perm`` so that ``pts_ref[perm[i]] ≈ pts_ours[i]``.

    Built by lexicographic sort with a tolerance via rounding.
    """
    scale = 1.0 / max(tol, 1e-15)
    a = np.round(pts_ours * scale).astype(np.int64)
    b = np.round(pts_ref * scale).astype(np.int64)
    if a.shape != b.shape:
        raise ValueError(f"point count mismatch: {a.shape[0]} vs {b.shape[0]}")
    # Build a dict from rounded coord tuple -> index in ref
    ref_map: dict[tuple[int, int, int], int] = {}
    for i in range(b.shape[0]):
        ref_map.setdefault((int(b[i, 0]), int(b[i, 1]), int(b[i, 2])), i)
    perm = np.empty(a.shape[0], dtype=np.int64)
    for i in range(a.shape[0]):
        k = (int(a[i, 0]), int(a[i, 1]), int(a[i, 2]))
        idx = ref_map.get(k)
        if idx is None:
            raise ValueError(f"point {i} = {pts_ours[i]} has no match in reference")
        perm[i] = idx
    return perm


def compare_cases(our_dir: str, ref_dir: str) -> dict:
    our_poly = os.path.join(our_dir, "constant", "polyMesh")
    ref_poly = os.path.join(ref_dir, "constant", "polyMesh")

    pts_ours = read_points(os.path.join(our_poly, "points"))
    pts_ref = read_points(os.path.join(ref_poly, "points"))
    own_ours = read_owner(os.path.join(our_poly, "owner"))
    own_ref = read_owner(os.path.join(ref_poly, "owner"))
    nei_ours = read_neighbour(os.path.join(our_poly, "neighbour"))
    nei_ref = read_neighbour(os.path.join(ref_poly, "neighbour"))
    ofs_ours, conn_ours = read_faces_any(os.path.join(our_poly, "faces"))
    ofs_ref, conn_ref = read_faces_any(os.path.join(ref_poly, "faces"))
    bnd_ours = read_boundary(os.path.join(our_poly, "boundary"))
    bnd_ref = read_boundary(os.path.join(ref_poly, "boundary"))

    n_cells_ours = int(max(own_ours.max(), nei_ours.max())) + 1
    n_cells_ref = int(max(own_ref.max(), nei_ref.max())) + 1

    # ANSA stores -1 for boundary faces in the neighbour file; OpenFOAM
    # stores only internal-face entries.  Filter accordingly.
    nei_ours_int = nei_ours[nei_ours >= 0]
    nei_ref_int = nei_ref[nei_ref >= 0]
    n_int_ours = int(nei_ours_int.size)
    n_int_ref = int(nei_ref_int.size)

    # Point permutation
    perm_ours_to_ref = _point_permutation(pts_ours, pts_ref)

    # Build cell → set of face-vertex-keys (in the reference point space)
    def cell_face_signature(own, nei_full, ofs, conn, perm, n_cells):
        cells: list[list[tuple[int, ...]]] = [[] for _ in range(n_cells)]
        keys = _face_keys(ofs, conn, perm)
        nf = own.size
        for f in range(nf):
            cells[int(own[f])].append(keys[f])
            n = int(nei_full[f]) if f < len(nei_full) else -1
            if n >= 0:
                cells[n].append(keys[f])
        return [tuple(sorted(s)) for s in cells]

    # Reconstruct full neighbour array for both
    nei_full_ours = np.full(own_ours.size, -1, dtype=np.int64)
    nei_full_ours[:n_int_ours] = nei_ours_int
    nei_full_ref = np.full(own_ref.size, -1, dtype=np.int64)
    # For ANSA reference, neighbour file already has full nFaces length
    if nei_ref.size == own_ref.size:
        nei_full_ref[:] = nei_ref
    else:
        nei_full_ref[:n_int_ref] = nei_ref_int

    sig_ours = cell_face_signature(own_ours, nei_full_ours, ofs_ours, conn_ours,
                                   perm_ours_to_ref, n_cells_ours)
    sig_ref = cell_face_signature(own_ref, nei_full_ref, ofs_ref, conn_ref,
                                  np.arange(pts_ref.shape[0]), n_cells_ref)

    cells_match = Counter(sig_ours) == Counter(sig_ref)

    # Patch comparison (by name)
    patches_ours_by_name = {p["name"]: p for p in bnd_ours}
    patches_ref_by_name = {p["name"]: p for p in bnd_ref}
    common = set(patches_ours_by_name) & set(patches_ref_by_name)
    only_ours = set(patches_ours_by_name) - common
    only_ref = set(patches_ref_by_name) - common
    patch_face_counts_match = (
        Counter(p["nFaces"] for p in bnd_ours)
        == Counter(p["nFaces"] for p in bnd_ref)
    )

    return {
        "n_points_ours": pts_ours.shape[0],
        "n_points_ref":  pts_ref.shape[0],
        "n_cells_ours":  n_cells_ours,
        "n_cells_ref":   n_cells_ref,
        "n_faces_ours":  own_ours.size,
        "n_faces_ref":   own_ref.size,
        "n_int_faces_ours": n_int_ours,
        "n_int_faces_ref":  n_int_ref,
        "cells_topology_match": cells_match,
        "boundary_only_in_ours": sorted(only_ours),
        "boundary_only_in_ref":  sorted(only_ref),
        "boundary_common": sorted(common),
        "patch_face_count_multisets_match": patch_face_counts_match,
        "patches_ours": [(p["name"], p.get("type"), p.get("nFaces")) for p in bnd_ours],
        "patches_ref":  [(p["name"], p.get("type"), p.get("nFaces")) for p in bnd_ref],
    }


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: validate.py <our_case_dir> <reference_case_dir>", file=sys.stderr)
        return 2
    res = compare_cases(argv[0], argv[1])
    width = max(len(k) for k in res)
    for k, v in res.items():
        print(f"{k.ljust(width)} : {v}")
    return 0 if res["cells_topology_match"] else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main(sys.argv[1:]))
