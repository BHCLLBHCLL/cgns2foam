# AGENTS.md

## Project Overview

**cgns2foam** — Pure-Python converter (h5py + numpy) from CGNS/HDF5 meshes
to OpenFOAM v2412 case directories. Supports ANSA-compatible binary polyMesh
headers, cross-zone BC overlap trimming, coupling scan (`--scan`), and
`chtMultiRegionSimpleFoam` scaffolding (`--cht`).

## Cursor Cloud / agent instructions

- **Git LFS**: Files under `cases/` may be LFS-tracked. Run `git lfs pull`
  after clone if binaries are pointer stubs.
- **Python deps**: `pip install -r requirements.txt` (`h5py`, `numpy`).
- **Run converter** (repo root):
  ```bash
  python -m src path/to/case.cgns [out_dir]
  python -m src --scan path/to/case.cgns --report couplings.json
  python -m src --cht path/to/case.cgns out_cht          # mono + Allrun.pre
  python -m src --cht-direct path/to/case.cgns out_cht   # one-step multi-region
  ```
- **Tests**:
  ```bash
  python -m unittest tests.test_box tests.test_bc_overlap tests.test_couplings -v
  ```
  Note: `.gitignore` has `tests/*`; new test modules may need
  `git add -f tests/test_*.py`.
- **Test data**: Prefer `tests/*.cgns` when present (e.g.
  `laptop_thermal_steady_scaled_v3_orig_fix.cgns`). `cases/` archives are
  large LFS objects.
- **Docs**: Algorithm and CHT scan details in `docs/TECHNICAL.md`
  (§3.5 BC trim, §3.6 couplings / CHT).
