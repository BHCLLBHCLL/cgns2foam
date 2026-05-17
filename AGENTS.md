# AGENTS.md

## Project Overview

**cgns2foam** — A Python project to convert CGNS (CFD General Notation System) mesh files to OpenFOAM project format.

The repository is in an early stage. It currently contains only test data (CGNS/HDF5 mesh files under `cases/`) stored via Git LFS, and no application source code yet.

## Cursor Cloud specific instructions

- **Git LFS**: All files under `cases/` are tracked by Git LFS. After cloning or pulling, run `git lfs pull` to fetch the actual binary data (CGNS files are ~66 MB each). Without this, the files will be LFS pointer stubs.
- **Python + h5py**: The CGNS files are HDF5 format. Use `h5py` to read and inspect them (`pip install h5py`). The update script handles this dependency.
- **No application code yet**: There is no script, build system, test suite, or linter configured. As code is added, update this file with run/test/lint instructions.
- **Test data**: Three test cases are available under `cases/`: `box_ansa`, `laptop_simplified`, and `tr03`. Each contains a `.cgns` mesh file and a `.zip` archive.
