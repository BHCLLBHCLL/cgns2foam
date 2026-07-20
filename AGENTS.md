# AGENTS.md

## Project Overview

**cgns2foam** — Pure-Python converter (h5py + numpy) from CGNS/HDF5 meshes
to OpenFOAM v2412 case directories. Supports ANSA-compatible binary polyMesh
headers, cross-zone BC overlap trimming, coupling scan (`--scan`), and
one-step multi-region `chtMultiRegionSimpleFoam` (`--cht-direct`).

## Cursor Cloud / agent instructions

- **Git LFS**: Files under `cases/` may be LFS-tracked. Run `git lfs pull`
  after clone if binaries are pointer stubs.
- **Python deps**: `pip install -r requirements.txt` (`h5py`, `numpy`).
- **CHT mode** (`--cht-direct`) requires a sidecar
  `<cgns-basename>.json`. Minimal format (see
  `tests/laptop_thermal_steady_scaled_v3_orig_BCs_fix.json`):
  ```json
  {"fluid_regions": ["zone.a", "..."], "solid_regions": ["zone.b", "..."]}
  ```
  All fluid zones go into one region `air` (`constant/air/polyMesh`),
  not a separate `fluid` region.
  Fluid–fluid → `cyclicAMI`; fluid–solid / solid–solid → `mappedWall`.
  Optional keys: `mrf_regions`, `heat_sources` (total W →
  `scalarSemiImplicitSource` + `volumeMode absolute`), `materials`
  (per-region rho/Cp/kappa, fluid mu/Pr/Cp), `external_convection`
  (regex patches → `externalWallHeatFluxTemperature`), `g`,
  `initial_conditions`, `n_procs`, `endTime`/`writeInterval`/`purgeWrite`.
  The old two-stage `--cht` (mono + splitMeshRegions) mode was removed —
  its fluid–solid interfaces were never converted to `mappedWall`.
- **Run converter** (repo root):
  ```bash
  python -m src path/to/case.cgns [out_dir]
  python -m src --scan path/to/case.cgns --report couplings.json
  python -m src --cht-direct path/to/case.cgns out_cht   # needs case.json
  ```
- **Tests**:
  ```bash
  python -m unittest tests.test_box tests.test_bc_overlap tests.test_couplings tests.test_regions_config -v
  ```
  Note: `.gitignore` has `tests/*`; new test modules may need
  `git add -f tests/test_*.py`.
- **Test data**: Prefer `tests/*.cgns` when present (e.g.
  `laptop_thermal_steady_scaled_v3_orig_foam.cgns`). `cases/` archives are
  large LFS objects.
- **Docs**: Algorithm and CHT scan details in `docs/TECHNICAL.md`
  (§3.5 BC trim, §3.6 couplings / CHT).
