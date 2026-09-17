# TG-FZ-R1 — Fixed-z Diagnostic Repair

Status: **PASS — diagnostic fixture repaired and hardened**.

This document records the scoped repair of
`tests/unit/test_dresp03tg_native_fixed_z.f90`. It does not add an energy
integration, a native pair or q-space evaluator, a Fe shell comparison, or any
DRESP-04 work.

## Defects confirmed and repaired

The fixture had two layout defects:

1. `gamma_spin_block` declared an `n x n` result and then wrote the down-spin
   block into indices `n+1:2*n`. The result is now an explicit `2*n x 2*n`
   collinear spin block, with checked global input shape and fixture-site
   ranges.
2. The two-site orbital structure-constant block was lifted to spin-major
   order through a shape-obscuring reuse of a `4*norb x 4*norb` array. The
   fixture now constructs `S_sites_orbital(2*norb,2*norb)` first and then
   explicitly forms
   `S_spin_major(4*norb,4*norb) = diag(S_sites_orbital,S_sites_orbital)`.

The screening shift `D = alpha - qpar` and finite-H endpoint width matrix are
also held in separate arrays. This prevents the finite-H construction from
overwriting the matrix needed by the P/S/path-operator covariance checks.

## Basis layout

For the controlled two-site fixture:

```text
norb             = (lmax+1)^2 = 9 for spd
nsite_fixture    = 2
norb_sites       = 2*norb = 18
nspin_orb_sites  = 4*norb = 36
```

The global spin-major ranges are:

```text
up/site1       1:norb
up/site2       norb+1:2*norb
down/site1     2*norb+1:3*norb
down/site2     3*norb+1:4*norb
```

The orbital two-site matrix uses `(site1,site2)` ordering. The spin lift is
block diagonal in spin, so the fixture contains no `x` or `y` spin blocks.

## Checks and covariance ledger

The test evaluates all values at the three existing fixed complex energies.
It reports the two ordered contractions independently:

```text
A = Im Tr[d_i Gup_ij d_j Gdown_ji]
B = Im Tr[d_i Gdown_ij d_j Gup_ji]
S = (A+B)/2
```

The gamma coefficient-space, active-alpha path-operator, and finite-H routes
all print `ordered_ud`, `ordered_du`, and `symmetrized`. The historical Pauli
helper is printed separately, together with its residual against finite-H
`ud`, `du`, and `S`.

The checked run produced the following compact ledger (printed values are
rounded here):

| `z` | gamma `ud` | alpha `ud` | finite-H `ud` | finite-H `du` | finite-H `S` | historical Pauli |
|---|---:|---:|---:|---:|---:|---:|
| `-0.91+0.83i` | `-5.380807e-6` | `-3.582616e-6` | `-5.380807e-6` | `-5.380807e-6` | `-5.380807e-6` | `-5.380807e-6` |
| `-0.17+0.04i` | `-2.202285e+0` | `-2.125794e+0` | `-2.202285e+0` | `-2.202285e+0` | `-2.202285e+0` | `-2.202285e+0` |
| `0.62+0.31i` | `-3.412346e-4` | `-1.377152e-3` | `-3.412346e-4` | `-3.412346e-4` | `-3.412346e-4` | `-3.412346e-4` |

The two orderings happen to coincide in this symmetric controlled fixture,
but they are evaluated by separate matrix products and are not replaced by a
fitted factor.

The maximum residuals were:

| check | maximum residual |
|---|---:|
| complex `P(z)` real-axis reduction | `0.000000e+0` |
| `d_matrix` versus width-scaled `DeltaP` | `2.827599e-16` |
| direct one-site native `[P_alpha-S_alpha]^-1` solve | `2.714385e-16` |
| fixed-z global native `[P_alpha-S_alpha]^-1` solve | `4.518280e-16` |
| P transformation | `2.842171e-14` |
| S transformation round trip | `2.220446e-15` |
| off-site path-operator covariance | `6.497414e-16` |
| gamma versus finite-H, ordered `ud` | `4.440892e-16` |
| gamma versus finite-H, ordered `du` | `4.440892e-16` |
| alpha versus gamma, ordered `ud` | `7.649156e-2` |
| alpha versus gamma, ordered `du` | `7.649156e-2` |
| alpha versus finite-H, either ordering | `7.649156e-2` |
| Pauli helper versus finite-H `ud`, `du`, or `S` | `4.440892e-16` |
| Pauli helper versus independent explicit Pauli reconstruction | `0.000000e+0` |

The independent Pauli reconstruction uses the frozen historical semantics:

```text
d_i G0_ij d_j G0_ji
- d_i Gx_ij d_j Gx_ji
- d_i Gy_ij d_j Gy_ji
- d_i Gz_ij d_j Gz_ji
```

with `Gx=Gy=0`, `G0=(Gup+Gdown)/2`, and
`Gz=(Gup-Gdown)/2`. The matrix order is preserved, giving the algebraic
collinear reduction `S=(A+B)/2`. The historical helper now agrees with that
explicit reconstruction to roundoff; its residual against the ordered
contractions is reported rather than interpreted as a new physics result.

## Checked-build evidence

Compiler:

```text
/bin/gfortran 13.3.0
```

The checked configuration was GNU Fortran Debug with the relevant compile
flags including:

```text
-g -O0 -fbacktrace -fcheck=all,no-recursion
```

Commands:

```bash
cmake -S . -B build -DRUN_UNIT_TESTS=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j2
ctest -R 'UnitDresp03tgNativeFixedZ' --output-on-failure
ctest -R 'UnitDresp03|UnitLrKlStaticBridge' --output-on-failure
```

The full Debug build completed. The focused test passed under bounds/runtime
checking without array-bound, shape-assignment, or related runtime
diagnostics. The nearby regression selection passed all 8 tests.

## Scope and remaining blocker

The P transformation, S transformation, and off-site path-operator
transformation are now separate diagnostics and close at numerical roundoff.
The complete active-alpha versus gamma/finite-H exchange-vertex contraction
still differs by up to `7.649156e-2`. That residual remains exposed and is not
accepted by a sign, factor, scale, broadening, or normalization adjustment.

The spin-dependent screening/vertex covariance derivation is deferred to
TG-FZ-R2. The historical Pauli relation is numerically well-defined in R1;
its physical/native ordering interpretation remains a later campaign task.

`source/exchange.f90` was not modified. Its SHA-256 remains:

```text
6e9ae7da6af46367aaf18a4fef37bdd4a9a0d6f69fe8e08b589e5c5f278a70e6
```
