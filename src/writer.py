"""OpenFOAM case writer.

Generates a complete OpenFOAM project directory (``constant/polyMesh/*``,
``constant/turbulenceProperties``, ``system/controlDict`` and a minimal
``0/`` directory).

By default, polyMesh geometry files use **binary** with ANSA 25.1-style
headers (``location ""``, full-length ``neighbour`` with ``-1`` on boundary
faces, ANSA ``note`` strings, extra banner spacing).  ``faces`` uses binary
``faceCompactList`` like ``points``, ``owner`` and ``neighbour``.

Use :class:`WriteOptions` (or CLI ``--openfoam-native``) for the legacy
OpenFOAM-native layout (``location "constant/polyMesh"``, internal-face
``neighbour`` only).

Binary file layout (little-endian, ``WM_LABEL_SIZE = 32``,
``WM_PRECISION_OPTION = DP``):

    points       : <header>  \n N \n ( <N×3×float64> ) \n
    faces        : <header>  \n (N+1) \n ( <(N+1)×int32 offsets> )
                              <S> \n ( <S×int32 connectivity> ) \n
    owner        : <header>  \n N \n ( <N×int32> ) \n
    neighbour    : <header>  \n M \n ( <M×int32> ) \n

The ``boundary``, ``cellZones`` and ``faceZones`` files are dictionaries.
"""

from __future__ import annotations

import dataclasses
import datetime as _dt
import os
from typing import Iterable

import numpy as np

from .topology import Mesh, Patch


@dataclasses.dataclass(frozen=True)
class WriteOptions:
    """Controls polyMesh output format and ANSA compatibility."""

    mesh_format: str = "binary"
    """``"binary"`` (default, ANSA-compatible) or ``"ascii"``."""

    mesh_location: str = ""
    """FoamFile ``location`` for polyMesh entries (ANSA uses ``""``)."""

    full_neighbour: bool = True
    """When ``True``, write ``neighbour`` with ``nFaces`` entries and ``-1`` on boundaries."""

    ansa_headers: bool = True
    """Match ANSA 25.1 banner spacing, ``note`` strings and dictionary headers."""

    @classmethod
    def openfoam_native(cls) -> WriteOptions:
        """OpenFOAM v2412 native binary mesh (internal-face neighbour only)."""
        return cls(
            mesh_format="binary",
            mesh_location="constant/polyMesh",
            full_neighbour=False,
            ansa_headers=False,
        )

    # Backward-compatible alias
    openfoam_binary = openfoam_native

# ---------------------------------------------------------------------------
# Header helpers
# ---------------------------------------------------------------------------


# ANSA 25.1 OpenFOAM export uses an 86-column banner; its importer also
# keys off ``ANSA_VERSION`` / ``Output from:`` in that banner.
_ANSA_BANNER_WIDTH = 86


def _pad_ansa_line(text: str) -> str:
    inner = f"    {text}".ljust(_ANSA_BANNER_WIDTH)
    return f"|{inner}| \n"


def _banner(source_path: str, *, ansa_spacing: bool = False) -> str:
    """Render the comment banner at the top of every emitted file."""
    if ansa_spacing:
        ts = _dt.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Y")
        src = os.path.abspath(source_path).replace("\\", "/")
        w = _ANSA_BANNER_WIDTH
        return (
            "/*" + "-" * w + "*\\\n"
            + _pad_ansa_line("")
            + _pad_ansa_line("ANSA_VERSION: 25.1.0")
            + _pad_ansa_line("")
            + _pad_ansa_line(f"file created by  cgns2foam  {ts}")
            + _pad_ansa_line("")
            + _pad_ansa_line(f"Output from: {src}")
            + _pad_ansa_line("")
            + "\\*" + "-" * w + "*/\n\n\n\n"
        )

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
                 fmt: str = "binary", location: str = "",
                 *, ansa_spacing: bool = False) -> bytes:
    return (
        _banner(source_path, ansa_spacing=ansa_spacing)
        + _foam_dict_header(class_name, obj_name, fmt=fmt, location=location)
    ).encode("ascii")


# ---------------------------------------------------------------------------
# Binary primitives
# ---------------------------------------------------------------------------


def _write_binary_label_list(fh, values: np.ndarray) -> None:
    """Write ``<N>\\n(<N×int32>)`` to *fh*.

    *fh* must already end with the header's trailing blank line; do not
    emit an extra ``\\n`` before the size token (ANSA's line parser keys
    off the ``2275435`` / ``(`` line numbers).
    """
    arr = np.ascontiguousarray(values, dtype=np.int32)
    fh.write(f"{arr.size}\n(".encode("ascii"))
    fh.write(arr.tobytes(order="C"))
    fh.write(b")\n")


def _write_binary_scalar_list(fh, values: np.ndarray) -> None:
    arr = np.ascontiguousarray(values, dtype=np.float64)
    fh.write(f"{arr.size}\n(".encode("ascii"))
    fh.write(arr.tobytes(order="C"))
    fh.write(b")\n")


def _write_binary_vector_list(fh, values: np.ndarray) -> None:
    """Vector field: list of (x, y, z) stored as flat float64s."""
    arr = np.ascontiguousarray(values, dtype=np.float64).reshape(-1)
    n = arr.size // 3
    fh.write(f"{n}\n(".encode("ascii"))
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
    fh.write(f"{ofs.size}\n(".encode("ascii"))
    fh.write(ofs.tobytes(order="C"))
    fh.write(b")")
    fh.write(f"{con.size}\n(".encode("ascii"))
    fh.write(con.tobytes(order="C"))
    fh.write(b")\n")


# ---------------------------------------------------------------------------
# ASCII primitives (ANSA-safe: no embedded 0x0A in numeric payloads)
# ---------------------------------------------------------------------------


def _write_ascii_label_list(fh, values: np.ndarray) -> None:
    """Write ``<N>\\n(\\n<N lines of int>\\n)``."""
    arr = np.ascontiguousarray(values, dtype=np.int64).ravel()
    fh.write(f"{arr.size}\n(\n".encode("ascii"))
    if arr.size:
        np.savetxt(fh, arr.reshape(-1, 1), fmt="%d")
    fh.write(b")\n")


def _write_ascii_vector_list(fh, values: np.ndarray) -> None:
    arr = np.ascontiguousarray(values, dtype=np.float64).reshape(-1, 3)
    n = arr.shape[0]
    fh.write(f"{n}\n(\n".encode("ascii"))
    for x, y, z in arr:
        fh.write(f"({x:g} {y:g} {z:g})\n".encode("ascii"))
    fh.write(b")\n")


def _write_ascii_compact_label_list(fh, offsets: np.ndarray,
                                    connectivity: np.ndarray) -> None:
    """CompactListList ASCII form: offsets list immediately followed by values."""
    ofs = np.ascontiguousarray(offsets, dtype=np.int64).ravel()
    con = np.ascontiguousarray(connectivity, dtype=np.int64).ravel()
    fh.write(f"{ofs.size}\n(\n".encode("ascii"))
    if ofs.size:
        np.savetxt(fh, ofs.reshape(-1, 1), fmt="%d")
    fh.write(b")\n")
    fh.write(f"{con.size}\n(\n".encode("ascii"))
    if con.size:
        np.savetxt(fh, con.reshape(-1, 1), fmt="%d")
    fh.write(b")\n")


def _neighbour_values(mesh: Mesh, options: WriteOptions) -> np.ndarray:
    if options.full_neighbour:
        full = np.full(mesh.owner.size, -1, dtype=np.int32)
        full[: mesh.n_internal_faces] = mesh.neighbour
        return full
    return mesh.neighbour


def _mesh_fmt(options: WriteOptions) -> str:
    return options.mesh_format


def _mesh_note(mesh: Mesh, options: WriteOptions) -> str:
    if options.ansa_headers:
        return (
            f"nCells:{mesh.n_cells} "
            f"nActiveFaces:{mesh.owner.size} "
            f"nActivePoints:{mesh.points.shape[0]}"
        )
    return (
        f"nPoints:{mesh.points.shape[0]} "
        f"nCells:{mesh.n_cells} "
        f"nFaces:{mesh.owner.size} "
        f"nInternalFaces:{mesh.n_internal_faces}"
    )


def _dict_header_fmt(options: WriteOptions, *, ascii_body: bool = False) -> str:
    """FoamFile ``format`` token for dictionary-like polyMesh entries."""
    if ascii_body and not options.ansa_headers:
        return "ascii"
    return options.mesh_format


def _case_file_location(options: WriteOptions, default: str) -> str:
    """FoamFile ``location`` for system/constant/0 entries."""
    return "" if options.ansa_headers else default


def _case_file_header_fmt(options: WriteOptions) -> str:
    """FoamFile ``format`` for non-polyMesh case files (ANSA uses ``binary``)."""
    return "binary" if options.ansa_headers else "ascii"


def _write_compression_token(options: WriteOptions) -> str:
    return "uncompressed" if options.ansa_headers else "off"


def _owner_neighbour_header(source: str, obj_name: str, mesh: Mesh,
                            options: WriteOptions) -> bytes:
    fmt = _mesh_fmt(options)
    note = _mesh_note(mesh, options)
    if options.ansa_headers:
        class_note = (
            "\tclass labelList;\n"
            "\n"
            f"\tnote \"{note}\";\n"
        )
    else:
        class_note = (
            "\tclass labelList;\n"
            f"\tnote \"{note}\";\n"
        )
    return (
        _banner(source, ansa_spacing=options.ansa_headers).encode("ascii")
        + (
            "FoamFile\n"
            "{\n"
            "\tversion 2.0;\n"
            f"\tformat {fmt};\n"
            f"{class_note}"
            f"\tlocation \"{options.mesh_location}\";\n"
            f"\tobject {obj_name};\n"
            "}\n"
            "/*" + "-" * 75 + "*/\n"
            "/*" + "-" * 75 + "*/\n\n"
        ).encode("ascii")
    )


# ---------------------------------------------------------------------------
# Mesh files
# ---------------------------------------------------------------------------


def _write_points(path: str, mesh: Mesh, source: str,
                  options: WriteOptions) -> None:
    fmt = _mesh_fmt(options)
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "vectorField", "points",
                              fmt=fmt, location=options.mesh_location,
                              ansa_spacing=options.ansa_headers))
        if fmt == "ascii":
            _write_ascii_vector_list(fh, mesh.points)
        else:
            _write_binary_vector_list(fh, mesh.points)


def _write_faces(path: str, mesh: Mesh, source: str,
                 options: WriteOptions) -> None:
    fmt = _mesh_fmt(options)
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "faceCompactList", "faces",
                              fmt=fmt, location=options.mesh_location,
                              ansa_spacing=options.ansa_headers))
        if fmt == "ascii":
            _write_ascii_compact_label_list(fh, mesh.face_offsets, mesh.face_vertices)
        else:
            _write_binary_compact_label_list(fh, mesh.face_offsets, mesh.face_vertices)


def _write_owner(path: str, mesh: Mesh, source: str,
                 options: WriteOptions) -> None:
    fmt = _mesh_fmt(options)
    with open(path, "wb") as fh:
        fh.write(_owner_neighbour_header(source, "owner", mesh, options))
        if fmt == "ascii":
            _write_ascii_label_list(fh, mesh.owner)
        else:
            _write_binary_label_list(fh, mesh.owner)


def _write_neighbour(path: str, mesh: Mesh, source: str,
                     options: WriteOptions) -> None:
    fmt = _mesh_fmt(options)
    values = _neighbour_values(mesh, options)
    with open(path, "wb") as fh:
        fh.write(_owner_neighbour_header(source, "neighbour", mesh, options))
        if fmt == "ascii":
            _write_ascii_label_list(fh, values)
        else:
            _write_binary_label_list(fh, values)


def _write_boundary(path: str, mesh: Mesh, source: str,
                    options: WriteOptions) -> None:
    lines: list[str] = [
        f"{len(mesh.patches)}\n(\n",
    ]
    for p in mesh.patches:
        if options.ansa_headers:
            lines.append(
                f"\n\t{p.name}\n\t{{\n"
                f"\t\ttype {p.bc_type};\n"
                f"\t\tstartFace {p.start_face};\n"
                f"\t\tnFaces {p.n_faces};\n"
                f"\t}}\n"
            )
        else:
            lines.append(
                f"\n\t{p.name}\n\t{{\n"
                f"\t\ttype {p.bc_type};\n"
                f"\t\tnFaces {p.n_faces};\n"
                f"\t\tstartFace {p.start_face};\n"
                f"\t}}\n"
            )
    lines.append("\n)\n")
    hdr_fmt = _dict_header_fmt(options, ascii_body=True)
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "polyBoundaryMesh", "boundary", fmt=hdr_fmt,
                              location=options.mesh_location,
                              ansa_spacing=options.ansa_headers))
        fh.write("".join(lines).encode("ascii"))


def _write_cell_zones(path: str, mesh: Mesh, source: str,
                      options: WriteOptions) -> None:
    fmt = _mesh_fmt(options)
    # OpenFOAM v2412 (openfoam.com) expects ``class regIOobject`` here;
    # the openfoam.org branch uses ``cellZoneList`` instead.  We target
    # v2412, which is also what ANSA emits.
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "regIOobject", "cellZones",
                              fmt=fmt, location=options.mesh_location,
                              ansa_spacing=options.ansa_headers))
        # Filter out empty zones
        zones = [z for z in mesh.cell_zones if z.cell_labels.size > 0]
        fh.write(f"{len(zones)}\n(\n".encode("ascii"))
        for z in zones:
            fh.write(f"\t{z.name}\n\t{{\n".encode("ascii"))
            fh.write(b"\t\ttype cellZone;\n")
            fh.write(b"\t\tcellLabels\tList<label>")
            if fmt == "ascii":
                _write_ascii_label_list(fh, z.cell_labels)
            else:
                _write_binary_label_list(fh, z.cell_labels)
            fh.write(b"\t;\n\t}\n")
        fh.write(b")\n")


def _write_face_zones(path: str, source: str, options: WriteOptions) -> None:
    """Empty faceZones – kept for consistency with the ANSA reference."""
    hdr_fmt = _dict_header_fmt(options, ascii_body=True)
    with open(path, "wb") as fh:
        fh.write(_full_header(source, "regIOobject", "faceZones", fmt=hdr_fmt,
                              location=options.mesh_location,
                              ansa_spacing=options.ansa_headers))
        fh.write(b"0\n(\n)\n")


# ---------------------------------------------------------------------------
# system / constant / 0
# ---------------------------------------------------------------------------


def _write_control_dict(path: str, source: str, options: WriteOptions) -> None:
    write_format = "binary" if options.ansa_headers else options.mesh_format
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
        f"writeFormat\t{write_format};\n\n"
        "writePrecision\t6;\n\n"
        f"writeCompression\t{_write_compression_token(options)};\n\n"
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
        fh.write(_full_header(
            source, "dictionary", "controlDict",
            fmt=_case_file_header_fmt(options),
            location=_case_file_location(options, "system"),
            ansa_spacing=options.ansa_headers,
        ))
        fh.write(body.encode("ascii"))


def _write_fv_schemes(path: str, source: str, options: WriteOptions) -> None:
    """Minimal ``system/fvSchemes`` accepted by OpenFOAM v2412.

    Required even for ``checkMesh``; users should replace the entries
    with solver-appropriate schemes before running a simulation.
    """
    body = (
        "\nddtSchemes\n{\n\tdefault\tsteadyState;\n}\n\n"
        "gradSchemes\n{\n\tdefault\tGauss linear;\n}\n\n"
        "divSchemes\n{\n\tdefault\tnone;\n}\n\n"
        "laplacianSchemes\n{\n\tdefault\tGauss linear corrected;\n}\n\n"
        "interpolationSchemes\n{\n\tdefault\tlinear;\n}\n\n"
        "snGradSchemes\n{\n\tdefault\tcorrected;\n}\n\n"
        "wallDist\n{\n\tmethod\tmeshWave;\n}\n"
    )
    with open(path, "wb") as fh:
        fh.write(_full_header(
            source, "dictionary", "fvSchemes",
            fmt=_case_file_header_fmt(options),
            location=_case_file_location(options, "system"),
            ansa_spacing=options.ansa_headers,
        ))
        fh.write(body.encode("ascii"))


def _write_fv_solution(path: str, source: str, options: WriteOptions) -> None:
    """Minimal ``system/fvSolution``."""
    body = (
        "\nsolvers\n{\n}\n\n"
        "SIMPLE\n{\n\tnNonOrthogonalCorrectors\t0;\n}\n"
    )
    with open(path, "wb") as fh:
        fh.write(_full_header(
            source, "dictionary", "fvSolution",
            fmt=_case_file_header_fmt(options),
            location=_case_file_location(options, "system"),
            ansa_spacing=options.ansa_headers,
        ))
        fh.write(body.encode("ascii"))


def _write_turbulence_properties(path: str, source: str,
                                 options: WriteOptions) -> None:
    # OpenFOAM v2412 accepts the long form ``simulationType RAS; RAS { … }``
    # which we keep here for parity with the ANSA reference.  The shorter
    # ``simulationType laminar;`` form is equally valid in v2412.
    body = (
        "\nsimulationType RAS;\n\n"
        "RAS\n{\n"
        "\tRASModel laminar;\n\n"
        "\tturbulence off;\n\n"
        "\tprintCoeffs off;\n\n"
        "}\n"
    )
    with open(path, "wb") as fh:
        fh.write(_full_header(
            source, "dictionary", "turbulenceProperties",
            fmt=_case_file_header_fmt(options),
            location=_case_file_location(options, "constant"),
            ansa_spacing=options.ansa_headers,
        ))
        fh.write(body.encode("ascii"))


def _write_initial_field(path: str, source: str, options: WriteOptions, *,
                         name: str, dims: str, is_vector: bool,
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
        fh.write(_full_header(
            source, klass, name,
            fmt=_case_file_header_fmt(options),
            location=_case_file_location(options, "0"),
            ansa_spacing=options.ansa_headers,
        ))
        fh.write("".join(lines).encode("ascii"))


def _write_initial_conditions(zero_dir: str, mesh: Mesh, source: str,
                              options: WriteOptions) -> None:
    _write_initial_field(
        os.path.join(zero_dir, "U"), source, options,
        name="U", dims="[0 1 -1 0 0 0 0]", is_vector=True,
        internal_value="( 0. 0. 0. )",
        patches=mesh.patches,
    )
    _write_initial_field(
        os.path.join(zero_dir, "p"), source, options,
        name="p", dims="[0 2 -2 0 0 0 0]", is_vector=False,
        internal_value="0.",
        patches=mesh.patches,
    )
    _write_initial_field(
        os.path.join(zero_dir, "p_rgh"), source, options,
        name="p_rgh", dims="[0 2 -2 0 0 0 0]", is_vector=False,
        internal_value="0.",
        patches=mesh.patches,
    )


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def write_case(out_dir: str, mesh: Mesh, source_path: str,
               options: WriteOptions | None = None) -> None:
    """Write a full OpenFOAM case at *out_dir*.

    Creates the standard layout::

        out_dir/
            system/controlDict
            constant/turbulenceProperties
            constant/polyMesh/{points, faces, owner, neighbour, boundary,
                               cellZones, faceZones}
            0/{U, p, p_rgh}

    By default, polyMesh files use binary format with ANSA-compatible
    headers (see :class:`WriteOptions`).
    """
    if options is None:
        options = WriteOptions()
    poly = os.path.join(out_dir, "constant", "polyMesh")
    sysd = os.path.join(out_dir, "system")
    cstd = os.path.join(out_dir, "constant")
    zerd = os.path.join(out_dir, "0")
    for d in (poly, sysd, cstd, zerd):
        os.makedirs(d, exist_ok=True)

    _write_points(os.path.join(poly, "points"), mesh, source_path, options)
    _write_faces(os.path.join(poly, "faces"), mesh, source_path, options)
    _write_owner(os.path.join(poly, "owner"), mesh, source_path, options)
    _write_neighbour(os.path.join(poly, "neighbour"), mesh, source_path, options)
    _write_boundary(os.path.join(poly, "boundary"), mesh, source_path, options)
    _write_cell_zones(os.path.join(poly, "cellZones"), mesh, source_path, options)
    _write_face_zones(os.path.join(poly, "faceZones"), source_path, options)

    _write_control_dict(os.path.join(sysd, "controlDict"), source_path, options)
    _write_fv_schemes(os.path.join(sysd, "fvSchemes"), source_path, options)
    _write_fv_solution(os.path.join(sysd, "fvSolution"), source_path, options)
    _write_turbulence_properties(
        os.path.join(cstd, "turbulenceProperties"), source_path, options
    )

    _write_initial_conditions(zerd, mesh, source_path, options)
