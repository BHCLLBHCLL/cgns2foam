"""CLI entry point: ``python -m src <input.cgns> [output_dir]``.

The package is imported as ``src`` because the source directory at
the repository root is named ``src/``; the CLI program name is still
``cgns2foam`` for user-facing messages.
"""

from __future__ import annotations

import argparse
import os
import sys

from .convert import convert_file


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="cgns2foam",
        description=(
            "Convert a CFD CGNS (HDF5) file to an OpenFOAM project directory. "
            "Built on top of h5py + numpy."
        ),
    )
    p.add_argument("cgns_file", help="Path to the input .cgns file")
    p.add_argument(
        "output_dir",
        nargs="?",
        default=None,
        help=(
            "Destination directory for the OpenFOAM case. Defaults to "
            "<dirname-of-input>/<basename-without-ext>/."
        ),
    )
    p.add_argument(
        "-q", "--quiet", action="store_true",
        help="Suppress progress messages on stdout.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    cgns_path = args.cgns_file
    if not os.path.isfile(cgns_path):
        print(f"error: not a file: {cgns_path}", file=sys.stderr)
        return 2
    out_dir = args.output_dir
    if out_dir is None:
        stem = os.path.splitext(os.path.basename(cgns_path))[0]
        out_dir = os.path.join(os.path.dirname(os.path.abspath(cgns_path)), stem)
    convert_file(cgns_path, out_dir, verbose=not args.quiet)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
