"""cgns2foam – CFD CGNS (HDF5) to OpenFOAM project converter.

This package provides a pure-Python converter that reads CGNS files
(stored in the HDF5/CPEX 0001 layout) using only ``h5py`` and ``numpy``
as runtime dependencies, and writes a complete OpenFOAM case directory.

The Python package lives in the ``src/`` directory; import paths are
``src.convert``, ``src.reader`` … and the CLI is invoked as
``python -m src``.

Public entry point: :func:`src.convert.convert_file`.
"""

from .convert import convert_file, scan_file  # noqa: F401
from .cht_direct import convert_cht_direct  # noqa: F401
from .reader import read_cgns      # noqa: F401
from .writer import WriteOptions   # noqa: F401
from .couplings import scan_couplings, format_coupling_summary  # noqa: F401

__all__ = [
    "convert_file",
    "convert_cht_direct",
    "scan_file",
    "read_cgns",
    "WriteOptions",
    "scan_couplings",
    "format_coupling_summary",
]
__version__ = "0.2.0"
