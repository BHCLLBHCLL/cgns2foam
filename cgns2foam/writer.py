"""OpenFOAM case writer.

Generates a complete OpenFOAM project directory (``constant/polyMesh/*``,
``constant/turbulenceProperties``, ``system/controlDict`` and a minimal
``0/`` directory).

Mesh files are written in **binary** format compatible with OpenFOAM's
default build (``WM_LABEL_SIZE = 32``, ``WM_PRECISION_OPTION = DP``);
all other files are written in ASCII.

Binary file layout (little-endian, matches what ANSA emits):

    points       : <header>  \n N \n ( <N×3×float64> ) \n
    faces        : <header>  \n (N+1) \n ( <(N+1)×int32 offsets> )
                              <S> \n ( <S×int32 connectivity> ) \n
    owner        : <header>  \n N \n ( <N×int32> ) \n
    neighbour    : <header>  \n M \n ( <M×int32> ) \n      (M = #internal faces)

The ``boundary``, ``cellZones`` and ``faceZones`` files are dictionaries.
``cellZones`` embeds a binary ``List<label>`` payload per entry.
"""

from __future__ import annotations

import datetime as _dt
import os
from typing import Iterable

import numpy as np

from .topology import Mesh, Patch

# ---------------------------------------------------------------------------
# Header helpers
# ---------------------------------------------------------------------------


def _banner(source_path: str) -> str:
    """Render the comment banner at the top of every emitted file.

    The opening and closing markers (``/*...*\\`` and ``\\*...*/``)
    match the OpenFOAM convention; this format is what ``checkMesh``
    and the rest of the OpenFOAM toolchain expect when scanning header
    comments.
    """
    ts = _dt.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Y UTC")
    pad = lambda s: f"|    {s}".ljust(85) + "| \n"
    return (
        "/*" + "-" * 84 + "*\\\n"
        + "|" + " " * 84 + "| \n"
        + pad("cgns2foam (h5py-based CGNS -> OpenFOAM converter)")
        + pad(f"file created at {ts}")
        + pad(f"source: {source_path}")
        + "|" + " " * 84 + "| \n"
        + "\\*" + "-" * 84 + "*/\n\n"
    )


def _foam_dict_header(class_name: str, obj_name: str, fmt: str = "binary",
                      location: str = "") -> str:
    return (
        "FoamFile\n"
        "{\n"
        "\tversion 2.0;\n"
        f"\tformat {fmt};\n"
        f"\tclass {class_name};\n"
        f"\tlocation \"{location}\";\n"
        f"\tobject {obj_name};\n"
        "}\n"
        "/*" + "-" * 75 + "*/\n"
        "/*" + "-" * 75 + "*/\n\n"
    )


def _full_header(source_path: str, class_name: str, obj_name: str,
                 fmt: str = "binary", location: str = "") -> bytes:
    return (
        _banner(source_path)
        + _foam_dict_header(class_name, obj_name, fmt=fmt, location=location)
    ).encode("ascii")


# ---------------------------------------------------------------------------
# Binary primitives
# ---------------------------------------------------------------------------


def _write_binary_label_list(fh, values: np.ndarray) -> None:
    """Write ``<N>\\n(<N×int32>)`` to *fh*."""
    arr = np.ascontiguousarray(values, dtype=np.int32)
    fh.write(f"\n{arr.size}\n(".encode("ascii"))
    fh.write(arr.tobytes(order="C"))
    fh.write(b")\n")


def _write_binary_scalar_list(fh, values: np.ndarray) -> None:
    arr = np.ascontiguousarray(values, dtype=np.float64)
    fh.write(f"\n{arr.size}\n(".encode("ascii"))
    fh.write(arr.tobytes(order="C"))
    fh.write(b")\n")


def _write_binary_vector_list(fh, values: np.ndarray) -> None:
    """Vector field: list of (x, y, z) stored as flat float64s."""
    arr = np.ascontiguousarray(values, dtype=np.float64).reshape(-1)
    n = arr.size // 3
    fh.write(f"\n{n}\n(".encode("ascii"))
    fh.write(arr.tobytes(order="C"))
    fh.write(b")\n")


def _write_binary_compact_label_list(fh, offsets: np.ndarray,
                                     connectivity: np.ndarray) -> None:
    """Write OpenFOAM's CompactListList<label> binary form.

    Layout::

        <nOffsets>\n(<nOffsets×int32>)
        <nConn>\n(<nConn×int32>)
    """
    ofs = np.ascontiguousarray(offsets, dtype=np.int32)
    con = np.ascontiguousarray(connectivity, dtype=np.int32)
    fh.write(f"\n{ofs.size}\n(".encode("ascii"))
    fh.write(ofs.tobytes(order="C"))
    fh.write(b")")
    fh.write(f"{con.size}\n(".encode("ascii"))
    fh.write(con.tobytes(order="C"))
    fh.write(b")\n")


# ---------------------------------------------------------------------------
# Mesh files
# ---------------------------------------------------------------------------


def _write_points(path: str, mesh: Mesh, source: str) -> None:
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "vectorField", "points",
                              location="constant/polyMesh"))
        _write_binary_vector_list(fh, mesh.points)


def _write_faces(path: str, mesh: Mesh, source: str) -> None:
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "faceCompactList", "faces",
                              location="constant/polyMesh"))
        _write_binary_compact_label_list(fh, mesh.face_offsets, mesh.face_vertices)


def _write_owner(path: str, mesh: Mesh, source: str) -> None:
    note = (
        f"nPoints:{mesh.points.shape[0]} "
        f"nCells:{mesh.n_cells} "
        f"nFaces:{mesh.owner.size} "
        f"nInternalFaces:{mesh.n_internal_faces}"
    )
    with open(path, "wb") as fh:
        fh.write(_banner(source).encode("ascii"))
        fh.write(
            (
                "FoamFile\n"
                "{\n"
                "\tversion 2.0;\n"
                "\tformat binary;\n"
                "\tclass labelList;\n"
                f"\tnote \"{note}\";\n"
                "\tlocation \"constant/polyMesh\";\n"
                "\tobject owner;\n"
                "}\n"
                "/*" + "-" * 75 + "*/\n"
                "/*" + "-" * 75 + "*/\n\n"
            ).encode("ascii")
        )
        _write_binary_label_list(fh, mesh.owner)


def _write_neighbour(path: str, mesh: Mesh, source: str) -> None:
    note = (
        f"nPoints:{mesh.points.shape[0]} "
        f"nCells:{mesh.n_cells} "
        f"nFaces:{mesh.owner.size} "
        f"nInternalFaces:{mesh.n_internal_faces}"
    )
    with open(path, "wb") as fh:
        fh.write(_banner(source).encode("ascii"))
        fh.write(
            (
                "FoamFile\n"
                "{\n"
                "\tversion 2.0;\n"
                "\tformat binary;\n"
                "\tclass labelList;\n"
                f"\tnote \"{note}\";\n"
                "\tlocation \"constant/polyMesh\";\n"
                "\tobject neighbour;\n"
                "}\n"
                "/*" + "-" * 75 + "*/\n"
                "/*" + "-" * 75 + "*/\n\n"
            ).encode("ascii")
        )
        _write_binary_label_list(fh, mesh.neighbour)


def _write_boundary(path: str, mesh: Mesh, source: str) -> None:
    lines: list[str] = [
        f"\n{len(mesh.patches)}\n(\n",
    ]
    for p in mesh.patches:
        lines.append(
            f"\n\t{p.name}\n\t{{\n"
            f"\t\ttype {p.bc_type};\n"
            f"\t\tnFaces {p.n_faces};\n"
            f"\t\tstartFace {p.start_face};\n"
            f"\t}}\n"
        )
    lines.append("\n)\n")
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "polyBoundaryMesh", "boundary", fmt="ascii",
                              location="constant/polyMesh"))
        fh.write("".join(lines).encode("ascii"))


def _write_cell_zones(path: str, mesh: Mesh, source: str) -> None:
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "cellZoneList", "cellZones",
                              location="constant/polyMesh"))
        # Filter out empty zones
        zones = [z for z in mesh.cell_zones if z.cell_labels.size > 0]
        fh.write(f"\n{len(zones)}\n(\n".encode("ascii"))
        for z in zones:
            fh.write(f"\t{z.name}\n\t{{\n".encode("ascii"))
            fh.write(b"\t\ttype cellZone;\n")
            fh.write(b"\t\tcellLabels\tList<label>")
            _write_binary_label_list(fh, z.cell_labels)
            fh.write(b"\t;\n\t}\n")
        fh.write(b")\n")


def _write_face_zones(path: str, source: str) -> None:
    """Empty faceZones – kept for consistency with reference output."""
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "faceZoneList", "faceZones", fmt="ascii",
                              location="constant/polyMesh"))
        fh.write(b"\n0\n(\n)\n")


# ---------------------------------------------------------------------------
# system / constant / 0
# ---------------------------------------------------------------------------


def _write_control_dict(path: str, source: str) -> None:
    body = (
        "\napplication UserSolver;\n\n"
        "startFrom startTime;\n\n"
        "startTime\t0.;\n\n"
        "stopAt endTime;\n\n"
        "endTime\t1.;\n\n"
        "deltaT \t1.;\n\n"
        "writeControl timeStep;\n\n"
        "writeInterval\t1.;\n\n"
        "purgeWrite\t0;\n\n"
        "writeFormat\tbinary;\n\n"
        "writePrecision\t6;\n\n"
        "writeCompression\tuncompressed;\n\n"
        "timeFormat\tgeneral;\n\n"
        "timePrecision\t6;\n\n"
        "graphFormat\traw;\n\n"
        "runTimeModifiable\tyes;\n\n"
        "adjustTimeStep\toff;\n\n"
        "maxCo\t1.;\n\n"
        "maxAlphaCo\t1.;\n\n"
        "maxDeltaT\t1.;\n\n"
        "functions {\n}\n"
    )
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "dictionary", "controlDict", fmt="ascii",
                              location="system"))
        fh.write(body.encode("ascii"))


def _write_turbulence_properties(path: str, source: str) -> None:
    body = (
        "\nsimulationType RAS;\n\n"
        "RAS\n{\n"
        "\tRASModel laminar;\n\n"
        "\tturbulence off;\n\n"
        "\tprintCoeffs off;\n\n"
        "}\n"
    )
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "dictionary", "turbulenceProperties",
                              fmt="ascii", location="constant"))
        fh.write(body.encode("ascii"))


def _write_initial_field(path: str, source: str, *, name: str,
                         dims: str, is_vector: bool,
                         internal_value: str, patches: Iterable[Patch],
                         patch_bc: dict | None = None) -> None:
    """Generic 0/ field writer."""
    if patch_bc is None:
        patch_bc = {}
    klass = "volVectorField" if is_vector else "volScalarField"
    lines = [
        f"\ndimensions {dims};\n\n",
        f"internalField uniform {internal_value};\n\n",
        "boundaryField\n{\n",
    ]
    for p in patches:
        rule = patch_bc.get(p.bc_type, patch_bc.get("__default__"))
        if rule is None:
            # Sensible defaults
            if is_vector:
                if p.bc_type == "wall":
                    rule = "type fixedValue;\n\t\tvalue\t uniform (0. 0. 0.);"
                elif p.bc_type == "symmetryPlane":
                    rule = "type symmetryPlane;"
                elif p.bc_type == "empty":
                    rule = "type empty;"
                else:
                    rule = "type zeroGradient;"
            else:
                if p.bc_type == "symmetryPlane":
                    rule = "type symmetryPlane;"
                elif p.bc_type == "empty":
                    rule = "type empty;"
                else:
                    rule = "type zeroGradient;"
        lines.append(f"\t{p.name}\n\t{{\n\t\t{rule}\n\t}}\n\n")
    lines.append("}\n")
    with open(path, "wb") as fh:
        fh.write(_full_header(source, klass, name, fmt="ascii", location="0"))
        fh.write("".join(lines).encode("ascii"))


def _write_initial_conditions(zero_dir: str, mesh: Mesh, source: str) -> None:
    _write_initial_field(
        os.path.join(zero_dir, "U"), source,
        name="U", dims="[0 1 -1 0 0 0 0]", is_vector=True,
        internal_value="( 0. 0. 0. )",
        patches=mesh.patches,
    )
    _write_initial_field(
        os.path.join(zero_dir, "p"), source,
        name="p", dims="[0 2 -2 0 0 0 0]", is_vector=False,
        internal_value="0.",
        patches=mesh.patches,
    )
    _write_initial_field(
        os.path.join(zero_dir, "p_rgh"), source,
        name="p_rgh", dims="[0 2 -2 0 0 0 0]", is_vector=False,
        internal_value="0.",
        patches=mesh.patches,
    )


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def write_case(out_dir: str, mesh: Mesh, source_path: str) -> None:
    """Write a full OpenFOAM case at *out_dir*.

    Creates the standard layout::

        out_dir/
            system/controlDict
            constant/turbulenceProperties
            constant/polyMesh/{points, faces, owner, neighbour, boundary,
                               cellZones, faceZones}
            0/{U, p, p_rgh}
    """
    poly = os.path.join(out_dir, "constant", "polyMesh")
    sysd = os.path.join(out_dir, "system")
    cstd = os.path.join(out_dir, "constant")
    zerd = os.path.join(out_dir, "0")
    for d in (poly, sysd, cstd, zerd):
        os.makedirs(d, exist_ok=True)

    _write_points(os.path.join(poly, "points"), mesh, source_path)
    _write_faces(os.path.join(poly, "faces"), mesh, source_path)
    _write_owner(os.path.join(poly, "owner"), mesh, source_path)
    _write_neighbour(os.path.join(poly, "neighbour"), mesh, source_path)
    _write_boundary(os.path.join(poly, "boundary"), mesh, source_path)
    _write_cell_zones(os.path.join(poly, "cellZones"), mesh, source_path)
    _write_face_zones(os.path.join(poly, "faceZones"), source_path)

    _write_control_dict(os.path.join(sysd, "controlDict"), source_path)
    _write_turbulence_properties(
        os.path.join(cstd, "turbulenceProperties"), source_path
    )

    _write_initial_conditions(zerd, mesh, source_path)
