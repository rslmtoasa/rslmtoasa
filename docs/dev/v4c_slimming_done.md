# fable_v4c slimming: done note

Branch `cleanup/v4c-slim` (pushed as `cleanup/v4c-slim-rc`), base `origin/fable_v4c` at `9a87c09`. `source/` is unchanged
(0 diff lines against the base). The plan file `fable_v4c_slimming_plan.md` was not in the branch; the approved task text was the plan.
Everything deleted stays at the tags `archive/{fable_v4b,fable_v4b_u,fable_v4b_gbt,main}-2026-11-05`; the paths are in
`v4c_delete_manifest.txt` (411 entries, 342.8 MB).

| Measure | Before | After | Target |
|---|---|---|---|
| Tracked tree | 351.6 MB | 8.8 MB | <= 20 MB |
| `tests/` tracked | 17.4 MB | 2.5 MB | <= 5 MB |
| Registered tests (non-MPI, libXC on) | 285 | 112 (68 quick, 32 standard, 12 validation) | <= 130 |
| `quick` tier, `-j4` | over 10 min | 37 s native, 39 s CI configuration (largest test alone: 7.7 s) | <= 2 min, each <= 10 s |
| `standard` tier, `-j4` | not defined | 543 s native, 603 s CI configuration | <= 20 min |

| Family | Before | After |
|---|---|---|
| unit core | 49 | 36 |
| SCF example (incl. post-processing) | 54 | 30 |
| LDA+U | 29 | 10 |
| backend regression (incl. triads) | 18 | 9 |
| validation | 11 | 8 |
| spin response | 22 | 8 |
| GBT | 62 | 8 |
| XC | 40 | 3 |

Choices: the spin-response oracle `.nml` files are regenerated at test time (generators run in about 3 s; option 1);
the XC, spin-response, LDA+U local-moment and projector unit tests are consolidated into executables with sub-cases
(`tests/unit/subcase_runner.f90`); the LDA+U guards are one test (`tests/ldau/test_guards.py`); the reference binary of
`tools/ldau_regression.py` is replaced by JSON references (`tests/ldau/references/`, five rows registered).

Open items (details in `tests/KNOWN_ISSUES.md`):
- `UnitLdaULocalMoments`, sub-case `spinor_frame_kernel`, misses its bitwise check by 2.2e-16 with `-march=native` on FMA hardware;
  it passes in the CI configuration. The same failure is present on the unmodified `fable_v4c`.
- `LdaURefRsImpurityNsp2` (tier `validation`) differs by 1e-6 against the developer tolerance of 1e-7 on Linux.
- The stored LDA+U references are bound to the CI configuration (`-DENABLE_MARCH_NATIVE=OFF`).
- The merge of `main` into this branch was tried and stopped: 43 conflicting files, 52 hunks in 12 files under `source/`
  (including `self.f90`), merge base `9651a900`.
