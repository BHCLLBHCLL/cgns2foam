"""Topology layer: turn CGNS NGON/NFACE data into an OpenFOAM polyMesh.

Responsibilities:

* Build a per-zone face-vertex (NGON) array and a per-zone cell-face
  (NFACE) connectivity.
* Optionally derive an NGON/NFACE representation from fixed-shape
  element sections (TETRA_4, HEXA_8, …) when the CGNS file does not
  ship polyhedral element sections.
* Merge multiple CGNS zones into a single OpenFOAM mesh, keeping a
  cellZone per source zone.  Faces shared between two zones are *not*
  coalesced – the converter does not perform a geometric interface
  search; the user is expected to do that via OpenFOAM's
  ``mergeMeshes`` / ``stitchMesh`` if interfaces are required.
* Derive ``owner`` / ``neighbour`` arrays and apply the OpenFOAM
  upper-triangular ordering: internal faces first, sorted by
  ``(owner, neighbour)`` ascending; boundary faces grouped per patch.
* Flip the vertex order of faces whose normal would otherwise point
  from neighbour to owner.

The intermediate :class:`Mesh` object produced here is consumed by
:mod:`cgns2foam.writer`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .reader import (
    CGNSCase,
    CGNSElements,
    CGNSZone,
    HEXA_8,
    NFACE_n,
    NGON_n,
    PENTA_6,
    PYRA_5,
    QUAD_4,
    TETRA_4,
    TRI_3,
)


# ---------------------------------------------------------------------------
# Polyhedral mesh data container
# ---------------------------------------------------------------------------


@dataclass
class Patch:
    name: str
    bc_type: str          # OpenFOAM patch type, e.g. "wall", "patch"
    start_face: int = 0
    n_faces: int = 0


@dataclass
class CellZone:
    name: str
    cell_labels: np.ndarray   # 0-based cell indices


@dataclass
class Mesh:
    """An OpenFOAM-ready polyhedral mesh."""

    points: np.ndarray                  # (nPoints, 3) float64
    # Face vertex list stored in compact form (offsets + connectivity)
    face_offsets: np.ndarray            # (nFaces + 1,) int32
    face_vertices: np.ndarray           # (sum_face_size,) int32, 0-based
    owner: np.ndarray                   # (nFaces,) int32, 0-based
    neighbour: np.ndarray               # (nInternalFaces,) int32, 0-based
    n_internal_faces: int
    n_cells: int
    patches: list[Patch] = field(default_factory=list)
    cell_zones: list[CellZone] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Per-zone topology
# ---------------------------------------------------------------------------


# Canonical face decomposition of fixed-shape cells.
# Each entry is a list of (face_vertices) tuples in CGNS local-vertex
# numbering (1-based to match CGNS), with normals pointing outwards.
# Reference: CGNS SIDS §11.2, OpenFOAM cell shape conventions.
_FIXED_CELL_FACES: dict[int, list[list[int]]] = {
    TETRA_4: [
        [1, 3, 2],
        [1, 2, 4],
        [2, 3, 4],
        [3, 1, 4],
    ],
    PYRA_5: [
        [1, 4, 3, 2],
        [1, 2, 5],
        [2, 3, 5],
        [3, 4, 5],
        [4, 1, 5],
    ],
    PENTA_6: [
        [1, 2, 5, 4],
        [2, 3, 6, 5],
        [3, 1, 4, 6],
        [1, 3, 2],
        [4, 5, 6],
    ],
    HEXA_8: [
        [1, 4, 3, 2],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 4, 8, 7],
        [1, 5, 8, 4],
        [5, 6, 7, 8],
    ],
}


def _ngon_from_fixed(elements: list[CGNSElements]) -> tuple[CGNSElements, CGNSElements]:
    """Build synthetic NGON / NFACE arrays from fixed-shape cells.

    A unique face dictionary is built; identical faces (same sorted
    vertex set) are merged so internal faces are shared between cells.
    """
    if not elements:
        raise ValueError("Cannot build NGON/NFACE: no element sections provided")

    face_offsets: list[int] = [0]
    face_conn: list[int] = []
    cell_face_offsets: list[int] = [0]
    cell_face_conn: list[int] = []
    face_dict: dict[tuple[int, ...], int] = {}   # sorted vertex key -> face id (1-based)

    next_face_id = 1
    total_cells = 0

    for sec in sorted(elements, key=lambda e: e.erange[0]):
        if sec.etype not in _FIXED_CELL_FACES:
            # Treat 2D fixed shapes (TRI_3, QUAD_4) as faces – emitted
            # by NGON construction routine elsewhere if needed.  Cell
            # sections only.
            continue
        nvtx = {TETRA_4: 4, PYRA_5: 5, PENTA_6: 6, HEXA_8: 8}[sec.etype]
        n_cell = sec.erange[1] - sec.erange[0] + 1
        face_specs = _FIXED_CELL_FACES[sec.etype]
        conn = sec.connectivity.reshape(n_cell, nvtx)
        for ic in range(n_cell):
            cell_signed: list[int] = []
            for fspec in face_specs:
                fverts = [int(conn[ic, j - 1]) for j in fspec]
                key = tuple(sorted(fverts))
                fid = face_dict.get(key)
                if fid is None:
                    fid = next_face_id
                    face_dict[key] = fid
                    next_face_id += 1
                    face_conn.extend(fverts)
                    face_offsets.append(len(face_conn))
                    cell_signed.append(fid)               # outward normal
                else:
                    # Compare stored ordering – flip sign if reversed
                    start = face_offsets[fid - 1]
                    end = face_offsets[fid]
                    stored = face_conn[start:end]
                    sign = +1 if stored == fverts else -1
                    cell_signed.append(sign * fid)
            cell_face_conn.extend(cell_signed)
            cell_face_offsets.append(len(cell_face_conn))
            total_cells += 1

    n_faces = next_face_id - 1
    ngon = CGNSElements(
        name="GridElements_Faces_synth",
        etype=NGON_n,
        erange=(1, n_faces),
        connectivity=np.asarray(face_conn, dtype=np.int64),
        start_offset=np.asarray(face_offsets, dtype=np.int64),
    )
    nface = CGNSElements(
        name="Cells_synth",
        etype=NFACE_n,
        erange=(n_faces + 1, n_faces + total_cells),
        connectivity=np.asarray(cell_face_conn, dtype=np.int64),
        start_offset=np.asarray(cell_face_offsets, dtype=np.int64),
    )
    return ngon, nface


@dataclass
class _ZoneTopo:
    n_vertices: int
    n_cells: int
    n_faces: int
    face_offsets: np.ndarray         # (n_faces+1,) int64
    face_vertices: np.ndarray        # (sum_face_size,) int64, 0-based zone-local
    owner: np.ndarray                # (n_faces,) int64, 0-based zone-local cell id
    neighbour: np.ndarray            # (n_faces,) int64, 0-based or -1
    flip: np.ndarray                 # (n_faces,) bool – True ⇒ face must be reversed
    bc_face_lists: dict[str, np.ndarray]   # patch_name -> 0-based face ids (zone-local)


def _build_zone_topology(zone: CGNSZone) -> _ZoneTopo:
    ngon, nface = zone.ngon, zone.nface
    if ngon is None or nface is None:
        # Try to build from fixed-shape sections.
        ngon, nface = _ngon_from_fixed(zone.fixed_elements)

    face_offsets = ngon.start_offset.astype(np.int64, copy=False)
    face_conn = ngon.connectivity.astype(np.int64, copy=False) - 1  # 0-based vertex ids
    n_faces = ngon.erange[1] - ngon.erange[0] + 1
    face_first_id = ngon.erange[0]   # CGNS 1-based id of the first face

    nface_offsets = nface.start_offset.astype(np.int64, copy=False)
    nface_conn = nface.connectivity.astype(np.int64, copy=False)
    n_cells = nface.erange[1] - nface.erange[0] + 1
    cell_first_id = nface.erange[0]

    # For every reference (cell, face) in NFACE, the absolute value is a
    # CGNS 1-based face id and the sign indicates orientation
    # (positive = outward from this cell).
    abs_face = np.abs(nface_conn) - face_first_id            # 0-based face id
    sign = np.where(nface_conn > 0, 1, -1).astype(np.int8)

    counts = np.diff(nface_offsets).astype(np.int64)
    cell_for_ref = np.repeat(np.arange(n_cells, dtype=np.int64), counts)

    # Sort references by face id.  Argsort is stable so the order of the
    # two cells of an internal face is reproducible.
    order = np.argsort(abs_face, kind="stable")
    f_sorted = abs_face[order]
    c_sorted = cell_for_ref[order]
    s_sorted = sign[order]

    # Locate the first reference of every face.
    face_first = np.searchsorted(f_sorted, np.arange(n_faces))
    face_next = np.searchsorted(f_sorted, np.arange(1, n_faces + 1))
    n_refs = face_next - face_first

    owner = np.full(n_faces, -1, dtype=np.int64)
    neighbour = np.full(n_faces, -1, dtype=np.int64)
    flip = np.zeros(n_faces, dtype=bool)

    # Boundary faces (1 reference)
    mask1 = n_refs == 1
    if mask1.any():
        idx = face_first[mask1]
        c = c_sorted[idx]
        s = s_sorted[idx]
        owner[mask1] = c
        flip[mask1] = s < 0

    # Internal faces (2 references)
    mask2 = n_refs == 2
    if mask2.any():
        a = face_first[mask2]
        b = a + 1
        ca, cb = c_sorted[a], c_sorted[b]
        sa = s_sorted[a]
        # cell with outward normal (sign == +1)
        pos_cell = np.where(sa > 0, ca, cb)
        neg_cell = np.where(sa > 0, cb, ca)
        mn = np.minimum(pos_cell, neg_cell)
        mx = np.maximum(pos_cell, neg_cell)
        owner[mask2] = mn
        neighbour[mask2] = mx
        flip[mask2] = mn != pos_cell  # flip when OpenFOAM owner ≠ CGNS "+" cell

    # Faces with neither 1 nor 2 references are degenerate (should not happen).
    if (~(mask1 | mask2)).any():
        bad = np.where(~(mask1 | mask2))[0]
        raise ValueError(
            f"Zone {zone.name!r}: {bad.size} faces have an unexpected reference count "
            f"(must be 1 or 2); first offending CGNS face id = {int(bad[0]) + face_first_id}"
        )

    # Boundary condition → face id mapping (0-based, zone-local).
    bc_faces: dict[str, np.ndarray] = {}
    for bc in zone.bcs:
        if bc.grid_location != "FaceCenter":
            # Skip vertex / EdgeCenter / CellCenter BCs – we cannot
            # straightforwardly turn them into polyMesh patches.
            continue
        ids = bc.point_list.astype(np.int64, copy=False) - face_first_id
        # Only keep ids that are valid face references and lie on the boundary
        # (i.e. correspond to a single-cell face).
        valid = (ids >= 0) & (ids < n_faces) & (n_refs[np.clip(ids, 0, n_faces - 1)] == 1)
        ids = ids[valid]
        if ids.size:
            bc_faces[bc.name] = ids

    return _ZoneTopo(
        n_vertices=zone.n_vertices,
        n_cells=n_cells,
        n_faces=n_faces,
        face_offsets=face_offsets,
        face_vertices=face_conn,
        owner=owner,
        neighbour=neighbour,
        flip=flip,
        bc_face_lists=bc_faces,
    )


# ---------------------------------------------------------------------------
# Multi-zone merging and OpenFOAM ordering
# ---------------------------------------------------------------------------


def _unique_patch_name(base: str, used: set[str]) -> str:
    """Return a patch name not yet present in *used*, derived from *base*."""
    if base not in used:
        return base
    i = 1
    while True:
        cand = f"{base}_{i}"
        if cand not in used:
            return cand
        i += 1


def _sanitize_patch_name(name: str) -> str:
    """Make a patch name compatible with OpenFOAM dictionary keys."""
    out: list[str] = []
    for ch in name:
        if ch.isalnum() or ch in "_.":
            out.append(ch)
        else:
            out.append("_")
    result = "".join(out) or "patch"
    if result[0].isdigit():
        result = "p_" + result
    return result


def _bc_type_to_foam(bc_type: str) -> str:
    """Map a CGNS BC type string to an OpenFOAM patch type."""
    bt = (bc_type or "").lower()
    if bt.startswith("bcwall") or "wall" in bt:
        return "wall"
    if bt.startswith("bcsymmetryplane") or "symmetry" in bt:
        return "symmetryPlane"
    if "axis" in bt:
        return "empty"
    return "patch"


def build_mesh(case: CGNSCase, *, default_exterior_name: str = "default_exterior") -> Mesh:
    """Convert a :class:`CGNSCase` into an OpenFOAM-ready :class:`Mesh`."""
    if not case.zones:
        raise ValueError("CGNS case contains no zones")

    zone_topos = [_build_zone_topology(z) for z in case.zones]

    # ------------------------------------------------------------------
    # Concatenate zones (no interface stitching).
    # ------------------------------------------------------------------
    vtx_offsets = np.cumsum([0] + [zt.n_vertices for zt in zone_topos])
    cell_offsets = np.cumsum([0] + [zt.n_cells for zt in zone_topos])
    face_offsets_per_zone = np.cumsum([0] + [zt.n_faces for zt in zone_topos])

    n_points_total = int(vtx_offsets[-1])
    n_cells_total = int(cell_offsets[-1])
    n_faces_total = int(face_offsets_per_zone[-1])
    n_face_verts_total = int(sum(zt.face_vertices.size for zt in zone_topos))

    points = np.empty((n_points_total, 3), dtype=np.float64)
    owner = np.empty(n_faces_total, dtype=np.int64)
    neighbour = np.empty(n_faces_total, dtype=np.int64)
    flip = np.empty(n_faces_total, dtype=bool)

    face_offsets_local = np.empty(n_faces_total + 1, dtype=np.int64)
    face_offsets_local[0] = 0
    face_vertices = np.empty(n_face_verts_total, dtype=np.int64)

    fv_cursor = 0
    f_cursor = 0
    patches: list[tuple[str, str, np.ndarray]] = []   # (name, foam_type, global face ids)
    used_patch_names: set[str] = set()
    cell_zones: list[CellZone] = []

    for zi, (zone, zt) in enumerate(zip(case.zones, zone_topos)):
        # Points -------------------------------------------------------
        points[vtx_offsets[zi]:vtx_offsets[zi + 1], :] = zone.coords

        # Faces --------------------------------------------------------
        n_f = zt.n_faces
        n_fv = zt.face_vertices.size
        face_vertices[fv_cursor:fv_cursor + n_fv] = zt.face_vertices + vtx_offsets[zi]
        # rebuild offsets accumulated globally
        face_offsets_local[f_cursor + 1:f_cursor + n_f + 1] = (
            zt.face_offsets[1:] + fv_cursor
        )

        # Owner / neighbour shifted into global cell index space
        owner[f_cursor:f_cursor + n_f] = zt.owner + cell_offsets[zi]
        nb = zt.neighbour.copy()
        nb[nb >= 0] += cell_offsets[zi]
        neighbour[f_cursor:f_cursor + n_f] = nb
        flip[f_cursor:f_cursor + n_f] = zt.flip

        # Patches: assign per-BC face lists with disambiguated names.
        # A single face may legitimately appear in multiple CGNS BC nodes
        # (typical for rotor/stator interfaces in ANSA exports); OpenFOAM
        # requires each boundary face to belong to exactly one patch, so
        # we assign on a first-come, first-served basis.
        boundary_mask = nb == -1
        bc_assigned = np.zeros(n_f, dtype=bool)
        for bc in zone.bcs:
            if bc.name not in zt.bc_face_lists:
                continue
            local_ids = zt.bc_face_lists[bc.name]
            # ensure these ids really are boundary faces and not already
            # claimed by an earlier BC of this zone
            keep = boundary_mask[local_ids] & ~bc_assigned[local_ids]
            local_ids = local_ids[keep]
            if local_ids.size == 0:
                continue
            global_ids = local_ids + f_cursor
            base_name = _sanitize_patch_name(bc.name)
            patch_name = _unique_patch_name(base_name, used_patch_names)
            used_patch_names.add(patch_name)
            patches.append((patch_name, _bc_type_to_foam(bc.bc_type), global_ids))
            bc_assigned[local_ids] = True

        # Remaining boundary faces (not covered by any BC) → default patch
        remaining = np.where(boundary_mask & ~bc_assigned)[0]
        if remaining.size:
            default_name = _unique_patch_name(default_exterior_name, used_patch_names)
            used_patch_names.add(default_name)
            patches.append((default_name, "wall", remaining + f_cursor))

        # cellZone covering all cells from this zone
        cz_labels = np.arange(cell_offsets[zi], cell_offsets[zi + 1], dtype=np.int64)
        cell_zones.append(CellZone(name=_sanitize_patch_name(zone.name), cell_labels=cz_labels))

        fv_cursor += n_fv
        f_cursor += n_f

    # ------------------------------------------------------------------
    # Re-order faces: internal first (sorted by owner/neighbour),
    # then per patch (preserve patch face order).
    # ------------------------------------------------------------------
    is_internal = neighbour >= 0
    internal_ids = np.where(is_internal)[0]
    # Stable sort by (owner, neighbour)
    sort_key = owner[internal_ids].astype(np.int64) * (n_cells_total + 1) + neighbour[internal_ids]
    int_order = internal_ids[np.argsort(sort_key, kind="stable")]

    boundary_order_chunks: list[np.ndarray] = []
    final_patches: list[Patch] = []
    cursor = int_order.size

    for name, foam_type, face_ids in patches:
        if face_ids.size == 0:
            continue
        boundary_order_chunks.append(face_ids)
        final_patches.append(Patch(name=name, bc_type=foam_type,
                                   start_face=cursor, n_faces=int(face_ids.size)))
        cursor += int(face_ids.size)

    if boundary_order_chunks:
        bnd_order = np.concatenate(boundary_order_chunks)
    else:
        bnd_order = np.empty(0, dtype=np.int64)

    new_order = np.concatenate([int_order, bnd_order])
    assert new_order.size == n_faces_total, (
        f"face reorder mismatch: {new_order.size} vs {n_faces_total}"
    )

    # Apply reorder to owner / neighbour / flip / face_vertices.
    new_owner = owner[new_order].astype(np.int32, copy=False)
    new_neighbour = neighbour[new_order]
    new_flip = flip[new_order]

    # Rebuild face_vertices compactly in the new order.  Compute new
    # offsets by gathering the face sizes.
    face_sizes = (face_offsets_local[1:] - face_offsets_local[:-1])[new_order]
    new_face_offsets = np.empty(n_faces_total + 1, dtype=np.int64)
    new_face_offsets[0] = 0
    np.cumsum(face_sizes, out=new_face_offsets[1:])

    new_face_vertices = np.empty(int(new_face_offsets[-1]), dtype=np.int32)
    for new_fi, old_fi in enumerate(new_order):
        s = int(face_offsets_local[old_fi])
        e = int(face_offsets_local[old_fi + 1])
        verts = face_vertices[s:e]
        if new_flip[new_fi]:
            verts = verts[::-1]
        new_face_vertices[int(new_face_offsets[new_fi]):int(new_face_offsets[new_fi + 1])] = verts

    # neighbour array stored only for internal faces in OpenFOAM
    n_internal_faces = int(int_order.size)
    new_neighbour_internal = new_neighbour[:n_internal_faces].astype(np.int32, copy=False)

    return Mesh(
        points=points,
        face_offsets=new_face_offsets.astype(np.int32, copy=False),
        face_vertices=new_face_vertices,
        owner=new_owner,
        neighbour=new_neighbour_internal,
        n_internal_faces=n_internal_faces,
        n_cells=n_cells_total,
        patches=final_patches,
        cell_zones=cell_zones,
    )
