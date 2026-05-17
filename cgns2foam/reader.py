"""CGNS (HDF5) reader built on top of :mod:`h5py`.

The CGNS standard stores nodes in HDF5 with the following conventions
(see CPEX 0001 / SIDS-to-HDF5):

* Every CGNS node maps to an HDF5 group.
* The data payload of a node, when present, is stored in a child dataset
  named ``" data"`` (a literal name starting with a single space).
* The CGNS *label* (e.g. ``Zone_t``, ``Elements_t`` …) is stored in the
  group attribute ``label``.
* The CGNS *name* is the HDF5 link name.
* The data *type* (``MT``, ``I4``, ``R8``, ``C1`` …) is stored in the
  attribute ``type``.

This module exposes a small, narrowly-scoped reader that returns plain
Python dataclasses ready to be consumed by the topology / writer layers.
Only the subset of the standard that is needed to convert a typical
volume mesh (NGON_n / NFACE_n unstructured grids, boundary conditions,
optional flow solutions) is supported.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import h5py
import numpy as np

# CGNS element type enum values (subset used here).
NGON_n = 22
NFACE_n = 23
NODE = 2
BAR_2 = 3
TRI_3 = 5
QUAD_4 = 7
TETRA_4 = 10
PYRA_5 = 12
PENTA_6 = 14
HEXA_8 = 17

# Number of vertices per fixed-shape element type.
_ELEM_NVTX = {
    NODE: 1,
    BAR_2: 2,
    TRI_3: 3,
    QUAD_4: 4,
    TETRA_4: 4,
    PYRA_5: 5,
    PENTA_6: 6,
    HEXA_8: 8,
}

# Subnode name used by SIDS-to-HDF5 to store a node's primary data array.
_DATA = " data"


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------


def _label(node: h5py.Group) -> str:
    """Return the CGNS label of an HDF5 group (e.g. ``Zone_t``)."""
    a = node.attrs.get("label")
    if a is None:
        return ""
    if hasattr(a, "tobytes"):
        return a.tobytes().rstrip(b"\x00").decode("ascii", errors="replace")
    if isinstance(a, bytes):
        return a.rstrip(b"\x00").decode("ascii", errors="replace")
    return str(a)


def _data(node: h5py.Group) -> np.ndarray | None:
    """Return the payload array (subnode ``' data'``) of a CGNS group."""
    if isinstance(node, h5py.Group) and _DATA in node:
        return node[_DATA][()]
    return None


def _as_str(arr: np.ndarray | None) -> str:
    """Convert a CGNS ``C1`` char array to a python string."""
    if arr is None:
        return ""
    return arr.tobytes().rstrip(b"\x00").decode("ascii", errors="replace")


def _children(group: h5py.Group, label: str | None = None) -> Iterable[h5py.Group]:
    """Yield child *groups* of ``group``, optionally filtered by label.

    Datasets and CGNS internal payload nodes are skipped.
    """
    for name in group:
        # Skip raw HDF5 metadata datasets that occasionally sit beside groups.
        if name.startswith(" "):
            continue
        item = group[name]
        if not isinstance(item, h5py.Group):
            continue
        if label is None or _label(item) == label:
            yield item


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class CGNSBoundaryCondition:
    """A single CGNS BC node (``BC_t``)."""

    name: str
    bc_type: str
    grid_location: str  # e.g. "FaceCenter", "Vertex"
    # 1-based face / element indices that this BC applies to.
    point_list: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int32))


@dataclass
class CGNSElements:
    """A CGNS ``Elements_t`` block.

    ``etype`` is the element-type enum (see CGNS standard).  For
    fixed-shape element sections, ``start_offset`` is synthesised on
    the fly so the higher layers can iterate uniformly.
    """

    name: str
    etype: int
    erange: tuple[int, int]               # inclusive, 1-based
    connectivity: np.ndarray              # signed for NFACE_n, 1-based vertex/face ids
    start_offset: np.ndarray              # length = n_elem + 1


@dataclass
class CGNSZone:
    """A CGNS ``Zone_t`` node carrying unstructured mesh data."""

    name: str
    n_vertices: int
    n_cells: int
    coords: np.ndarray                     # shape (n_vertices, 3), float64
    ngon: CGNSElements | None              # face definitions (NGON_n)
    nface: CGNSElements | None             # cell definitions (NFACE_n)
    fixed_elements: list[CGNSElements] = field(default_factory=list)
    bcs: list[CGNSBoundaryCondition] = field(default_factory=list)


@dataclass
class CGNSCase:
    """A whole CGNS file, restricted to one ``CGNSBase_t`` node."""

    cell_dim: int
    phys_dim: int
    zones: list[CGNSZone] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Reader
# ---------------------------------------------------------------------------


def _read_elements(node: h5py.Group) -> CGNSElements:
    """Build a :class:`CGNSElements` from an ``Elements_t`` group."""
    name = node.name.rsplit("/", 1)[-1]
    etype, _ebnd = _data(node).tolist()
    erange = tuple(int(x) for x in _data(node["ElementRange"]))
    n_elem = erange[1] - erange[0] + 1
    conn = _data(node["ElementConnectivity"])
    if conn is None:
        raise ValueError(f"Elements {name!r} has no ElementConnectivity")
    if "ElementStartOffset" in node:
        # NGON_n / NFACE_n (and MIXED) use explicit offsets.
        so = _data(node["ElementStartOffset"])
    else:
        # Fixed-shape element type: synthesise the offset array.
        nvtx = _ELEM_NVTX.get(int(etype))
        if nvtx is None:
            raise NotImplementedError(
                f"Unsupported fixed-shape CGNS element type {etype}"
            )
        so = np.arange(0, (n_elem + 1) * nvtx, nvtx, dtype=np.int64)
    return CGNSElements(
        name=name,
        etype=int(etype),
        erange=erange,
        connectivity=np.ascontiguousarray(conn),
        start_offset=np.ascontiguousarray(so),
    )


def _read_bc(node: h5py.Group) -> CGNSBoundaryCondition:
    name = node.name.rsplit("/", 1)[-1]
    bc_type = _as_str(_data(node))
    grid_loc = "Vertex"
    if "GridLocation" in node:
        grid_loc = _as_str(_data(node["GridLocation"]))
    pl = np.empty(0, dtype=np.int32)
    if "PointList" in node:
        pl = _data(node["PointList"]).reshape(-1).astype(np.int64, copy=False)
    elif "PointRange" in node:
        pr = _data(node["PointRange"]).reshape(-1)
        pl = np.arange(int(pr[0]), int(pr[1]) + 1, dtype=np.int64)
    return CGNSBoundaryCondition(
        name=name, bc_type=bc_type, grid_location=grid_loc, point_list=pl
    )


def _read_zone(node: h5py.Group) -> CGNSZone | None:
    name = node.name.rsplit("/", 1)[-1]
    ztype_node = node.get("ZoneType")
    if ztype_node is None:
        return None
    ztype = _as_str(_data(ztype_node))
    if ztype != "Unstructured":
        # Structured zones are out of scope for this converter.
        raise NotImplementedError(
            f"Zone {name!r}: ZoneType {ztype!r} is not supported (only Unstructured)"
        )
    sz = _data(node).reshape(-1)
    n_vertices, n_cells = int(sz[0]), int(sz[1])

    # Coordinates ----------------------------------------------------------
    gc = node.get("GridCoordinates")
    if gc is None:
        raise ValueError(f"Zone {name!r}: missing GridCoordinates")
    cx = _data(gc.get("CoordinateX"))
    cy = _data(gc.get("CoordinateY"))
    cz_arr = _data(gc.get("CoordinateZ"))
    if cx is None or cy is None or cz_arr is None:
        raise ValueError(f"Zone {name!r}: missing CoordinateX/Y/Z")
    coords = np.column_stack(
        [cx.astype(np.float64, copy=False),
         cy.astype(np.float64, copy=False),
         cz_arr.astype(np.float64, copy=False)]
    )

    # Elements -------------------------------------------------------------
    ngon: CGNSElements | None = None
    nface: CGNSElements | None = None
    fixed: list[CGNSElements] = []
    for child in _children(node, label="Elements_t"):
        elem = _read_elements(child)
        if elem.etype == NGON_n:
            ngon = elem
        elif elem.etype == NFACE_n:
            nface = elem
        else:
            fixed.append(elem)

    # Boundary conditions --------------------------------------------------
    bcs: list[CGNSBoundaryCondition] = []
    zbc = node.get("ZoneBC")
    if zbc is not None:
        for child in _children(zbc, label="BC_t"):
            bcs.append(_read_bc(child))

    return CGNSZone(
        name=name,
        n_vertices=n_vertices,
        n_cells=n_cells,
        coords=coords,
        ngon=ngon,
        nface=nface,
        fixed_elements=fixed,
        bcs=bcs,
    )


def read_cgns(path: str) -> CGNSCase:
    """Read a CGNS/HDF5 file and return a :class:`CGNSCase`.

    Only the first ``CGNSBase_t`` node is consumed; additional bases
    raise :class:`NotImplementedError`.
    """
    with h5py.File(path, "r") as f:
        bases = [g for g in _children(f, label="CGNSBase_t")]
        if not bases:
            raise ValueError(f"{path!r}: no CGNSBase_t found")
        if len(bases) > 1:
            raise NotImplementedError(
                f"{path!r}: multiple CGNSBase_t nodes are not supported"
            )
        base = bases[0]
        dims = _data(base).reshape(-1)
        cell_dim, phys_dim = int(dims[0]), int(dims[1])
        zones: list[CGNSZone] = []
        for zn in _children(base, label="Zone_t"):
            zone = _read_zone(zn)
            if zone is not None:
                zones.append(zone)
    return CGNSCase(cell_dim=cell_dim, phys_dim=phys_dim, zones=zones)
