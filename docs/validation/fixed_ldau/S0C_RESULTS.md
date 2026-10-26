# S0c results

Work branch `fix/ldau-s0c`, base `a67715a`, merged locally as `7d28265c`. Reference binary: clean build of `90f57c3`
(sha256 `8dc179ed…`); candidate: build of the merged tree (sha256 `afa98866…`). All solver launches went through
`tools/run_rslmto_checked.py`, single-threaded.

## T0 regression (`tools/ldau_regression.py`, 4 SCF steps)

| check | measured | expected | command |
|---|---|---|---|
| no-U RS `bulk/bccFe` | 26 files byte-identical | byte-identical | `.venv/bin/python tools/ldau_regression.py --ref REF --cand build/bin/rslmto.x` |
| no-U KS `k_space_scf/bccFe` | 24 files byte-identical | byte-identical | same |
| +U RS `bccFe_hubbard_nsp1` | 26 files, 600 occupation values equal | by value | same |
| +U KS bccFe U=2 J=0.5 | 24 files, 600 values equal | by value | same |
| +U impurity `B2FeCo_nsp1` | 28 files, 600 values equal | by value | same |
| 13 guard fixtures (RA2), candidate | pass | pass | same |
| dm: legacy inline vs dm-backed restart | 26 files byte-identical | identical | same |
| dm: perturbed occupation (negative control) | outputs differ | must differ | same |
| dm: missing sidecar, hash mismatch | rc=1, fatal message | fatal | same |
| auto fixtures (`UnitAutoImpurityLocalMode`) | 1/1 | pass | same |
| CTest label `lda_u` | 11/12; `Val04RsImpurityEnergyDiagnostic` fails by 1.24e-12 Ry | 12/12 | same |

Found and fixed on the way: `radial_reason` printed uninitialised bytes in the HUBBARD_RUN line for explicit modes
(`self.f90`, one line); VAL-04, RA1 and RA2 read the removed inline `ldm_flatten` and an old log token.

## T1 Chebyshev, bcc Fe rc=80, nsp=1, local, J=0, window -1.0 to 1.3 Ry, 2500 points, 60 steps from accepted states

Chebyshev minus block (block lld=21, Chebyshev lld=150, beta 0.005):

| quantity | U=0 | U=2 eV | ratio |
|---|---|---|---|
| d occupation up | 3.31e-3 | -4.98e-3 | 1.5 |
| d occupation down | -4.64e-3 | 6.71e-3 | 1.4 |
| sphere moment (μB) | 8.27e-3 | -7.53e-3 | 0.91 |
| E_F (Ry) | 2.7e-5 | 9.75e-4 | 36 |
| E_Hub (Ry), U=2 | n/a | 4.46e-4 (0.109291 vs 0.108845) | n/a |

The S0c gate (U=2 difference at most 10x the U=0 difference) failed on E_F. T1(a) source audit: deterministic, no stochastic estimate.
T1(c) frozen replay was not run.

## T2 Ni bulk (user deck `RS-utils/tmp/Cheb/fccNi`, Chebyshev, rc=60, lld=150, nsp=1)

| run | steps | final Diff | moment (μB) | E_F (Ry) |
|---|---|---|---|---|
| U=0, restarted from the user's `Ni_out.nml` | 18 | 7.5e-9 | 0.670705 | -0.063509 |
| U=3 eV, J=0.7 eV, `auto` (resolved local) | 76 | 5.0e-9 | 0.786506 | -0.077760 |

U=3: d occupations 4.719841 up, 3.869480 down, E_Hub 0.09116 Ry. Logs carry six decimals; the 20-step spreads are at that resolution.

## Not checked

Stage-end full ctest; T1(c) replay; RD2/RD3(a) on the magnetic host; the acceptance spreads at full precision.

## Notes

- Quick tier 39/40: `Example_bulk_bccFe_nsp2_hubbard` stops at the SOC guard, identically at `90f57c3`.
- The reference RA2 row fails only because the updated test expects the sidecar format.
- Raw artifacts are under the session scratchpad; the JSON sidecar lists hashes and commands.
