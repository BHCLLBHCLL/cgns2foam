"""CLI entry point: ``python -m src <input.cgns> [output_dir]``.

The package is imported as ``src`` because the source directory at
the repository root is named ``src/``; the CLI program name is still
``cgns2foam`` for user-facing messages.
"""

from __future__ import annotations

import argparse
import os
import sys

from .convert import convert_file, scan_file
from .writer import WriteOptions


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="cgns2foam",
        description=(
            "Convert a CFD CGNS (HDF5) file to an OpenFOAM project directory. "
            "Optionally scan fluid/solid couplings and emit "
            "chtMultiRegionSimpleFoam scaffolding. "
            "Built on top of h5py + numpy."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  %(prog)s mesh.cgns out_case
  %(prog)s --scan mesh.cgns
  %(prog)s --scan mesh.cgns --report couplings.json
  %(prog)s --cht mesh.cgns out_cht_case
""".rstrip(),
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
    p.add_argument(
        "--openfoam-native",
        "--binary-mesh",
        action="store_true",
        dest="openfoam_native",
        help=(
            "Write polyMesh as OpenFOAM-native binary "
            "(location constant/polyMesh, internal-face neighbour only). "
            "Default output matches ANSA 25.1 OpenFOAM export conventions."
        ),
    )
    p.add_argument(
        "--scan",
        action="store_true",
        help=(
            "Scan CGNS zones for fluid/solid regions and coupling interface "
            "pairs (fluid-fluid / fluid-solid / solid-solid). Does not write "
            "a polyMesh unless --cht is also set."
        ),
    )
    p.add_argument(
        "--cht",
        action="store_true",
        help=(
            "After conversion, write chtMultiRegionSimpleFoam scaffolding "
            "from an automatic coupling scan (regionProperties, "
            "per-region thermo/0.orig, Allrun.pre). Implies a coupling scan."
        ),
    )
    p.add_argument(
        "--report",
        metavar="PATH",
        default=None,
        help="Write coupling scan JSON to PATH (default: stdout only for --scan; "
             "for --cht also writes <output>/coupling_scan.json).",
    )
    p.add_argument(
        "--solid-pattern",
        action="append",
        default=None,
        metavar="REGEX",
        help="Regex for solid zone names (repeatable). Default: solid_region / solid.*",
    )
    p.add_argument(
        "--fluid-pattern",
        action="append",
        default=None,
        metavar="REGEX",
        help="Regex for fluid zone names (repeatable).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    cgns_path = args.cgns_file
    if not os.path.isfile(cgns_path):
        print(f"error: not a file: {cgns_path}", file=sys.stderr)
        return 2

    if args.scan and not args.cht:
        report_path = args.report
        scan_file(
            cgns_path,
            report_path=report_path,
            verbose=not args.quiet,
            solid_patterns=args.solid_pattern,
            fluid_patterns=args.fluid_pattern,
        )
        return 0

    out_dir = args.output_dir
    if out_dir is None:
        stem = os.path.splitext(os.path.basename(cgns_path))[0]
        out_dir = os.path.join(os.path.dirname(os.path.abspath(cgns_path)), stem)

    write_opts = WriteOptions.openfoam_native() if args.openfoam_native else None
    convert_file(
        cgns_path,
        out_dir,
        verbose=not args.quiet,
        write_options=write_opts,
        cht=args.cht or False,
        solid_patterns=args.solid_pattern,
        fluid_patterns=args.fluid_pattern,
    )

    if args.cht and args.report:
        # Extra copy of the coupling report outside the case dir if requested
        from pathlib import Path
        import shutil
        src = Path(out_dir) / "coupling_scan.json"
        if src.is_file():
            Path(args.report).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, args.report)
            if not args.quiet:
                print(f"[cgns2foam] report copied to {args.report}")

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
