# GBT WP0 / G0 report

Date: 2026-08-03.  Scope: freeze q=0/noncollinear baselines and prevent the
legacy finite-q reciprocal GBT path from being interpreted as physics.

## Recorded baseline ledger

Machine-readable copy: `tests/gbt_wp0/baselines.json`.

| Case | Input / output source | Recorded result | Classification |
| --- | --- | --- | --- |
| Collinear bulk | `tests/scf/cases/bulk/bccFe`, `Example_bulk_bccFe_nsp2_block/ref.json` | `etot=-2541.9818992061555 Ry`, `mom(3)=1.0` | Golden ordinary q=0 regression |
| Periodic noncollinear bulk | same input with `nsp=4`, `Example_bulk_bccFe_nsp4_block/ref.json` | `etot=-2541.9818992061555 Ry`, `mom(3)=0.9999999999999996` | Golden periodic-NC q=0 regression |
| Explicit texture / commensurate supercell | `example/bulk/bccFe/q050/super/FM`, four site outputs | per-site `etot=-2541.9814920214894, -2541.9814919811238, -2541.9814919366772, -2541.9814918777261 Ry` | Diagnostic texture baseline; no reciprocal GBT claim |
| Real-space MFT | `example/bulk/bccFe/bccFe_gamma_h/frozen_magnon.dat` | at q=(0.01,0.01,0): `eband=-1.99664723 Ry`, `omega=2.02165721e-5` | Diagnostic MFT baseline; preserve path, not a new physics calibration |
| Old reciprocal GBT/MFT | `example/bulk/bccFe/q050/gbt/KS` | `etot=-2541.9697460049274 Ry`, `Eband=-2.1382627426 Ry` | Diagnostic only; explicitly not golden physics |

The existing q=0 ordinary references remain the only golden values in this
ledger.  Older finite-q GBT values are retained solely to make a change in the
migration observable.

## WP0 implementation and commands

Build:

```bash
cmake -S . -B build_13
cmake --build build_13 -j2
```

Ordinary q=0 regression:

```bash
/Users/andersb/envs/p311/bin/python tests/run_test.py --binary build_13/bin/rslmto.x \
  --cases-json tests/scf/cases.json --case-name Example_bulk_bccFe_nsp2_block \
  --scratch-root /tmp/rslmto-wp0 --compare-ref tests/scf/references
```

Finite-q diagnostic (expected failure):

```bash
/Users/andersb/envs/p311/bin/python tests/run_test.py --binary build_13/bin/rslmto.x \
  --cases-json tests/gbt_wp0/cases.json --case-name finite_q_gbt_soc_guard \
  --scratch-root /tmp/rslmto-wp0
```

Expected result: the q=0 setup proceeds; at finite q the code rebuilds a full
chemical BZ, and unsupported combinations stop before eigensolution. Finite-q
SOC, CCOR, Hubbard-V, HOH/overlap products, and every
non-`ham_only` reciprocal overlap mode stop earlier with an explicit WP0 fatal.
For every standard or generalized eigensolve, H(k), and O(k) if allocated, are
checked against their adjoints with a relative `1e-10` tolerance before LAPACK
is called.

Observed on this workspace:

- `build_13` was reconfigured and rebuilt successfully with GNU Fortran 13.4.0.
- `/Users/andersb/envs/p311/bin/python gbt_model.py` confirms the legacy transform is non-Hermitian away
  from q=0: at q=0.0208 and theta=90 degrees it reports
  `|H-H^H|=4.70e-02` (q=0 reports zero).
- The ordinary q=0 bcc-Fe run completed and matched 8/9 reference samples.
  The sole difference is `totaldos.out` row 1500, column 2:
  `21.60499` vs `21.60495` (`4e-5`). The untouched pre-existing
  `build_v13/bin/rslmto.x` reproduces exactly the same difference, so this is
  a gfortran-13/reference tolerance difference, not a WP0 change.
- `legacy_finite_q_gbt` passed its smoke check and logged full-BZ rebuilding;
  its nonzero-q runs used all 64 points of the requested 4x4x4 mesh.
- `finite_q_gbt_soc_guard` exited nonzero at the first finite-q point with
  `nonzero-q GBT with SOC is unsupported (WP0 guard)`.
- `UnitGammaSupercell` and `UnitDysonEquivalence` pass under `build_13`.

## Remaining risks

- The old `fourier_transform_gbt_array` is intentionally still present and is
  expected to fail the Hermiticity check at finite q; WP0 does not alter GBT
  physics.
- The exact bond gauge, multi-sublattice phase ownership, transformed overlap
  products, SOC, CCOR, and Hubbard-V remain WP1+ work.
- The reciprocal MFT band energy still uses projected-DOS moments; it is a
  diagnostic only until a projection-free occupied-eigenvalue sum is added.
- The stored ordinary DOS reference needs a gfortran-13-aware tolerance at the
  one observed sample before its strict 1e-6 comparison can be called green.
