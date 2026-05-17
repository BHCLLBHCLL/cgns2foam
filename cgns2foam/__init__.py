"""cgns2foam – CFD CGNS (HDF5) to OpenFOAM project converter.

This package provides a pure-Python converter that reads CGNS files
(stored in the HDF5/CPEX 0001 layout) using only ``h5py`` and ``numpy``
as runtime dependencies, and writes a complete OpenFOAM case directory.

Public entry point: :func:`cgns2foam.convert.convert_file`.
"""

from .convert import convert_file  # noqa: F401
from .reader import read_cgns      # noqa: F401

__all__ = ["convert_file", "read_cgns"]
__version__ = "0.1.0"
