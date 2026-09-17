# TG-FZ-R5 — Production-H / native-gamma representation bridge

Status: **BLOCKED at the static gamma-to-production map**.

This audit follows TG-FZ-R4 without changing the accepted production
Hamiltonian or the native exchange route.  It stops before mixed derivatives,
complete curvature, energy integration, native `J`/`J(q)`, Fe-shell checks, or
DRESP-04.

## Scope and invariants

The diagnostic is in
`tests/unit/test_dresp03tg_native_fixed_z.f90`.  It uses the existing
two-site `spd` fixture, the production screened structure-constant block, the
live `B/Q/E_nu` assembly seam, and the three fixed complex energies from R1/R2.
No fitting, sign repair, scale tuning, or change to `source/exchange.f90` was
made.

The reported `max` values below are maximum absolute matrix-entry norms.  The
truncation is also reported with a Frobenius norm.

## Live notation map

The source trace is `source/symbolic_atom.f90:270-311`:

```text
wow       = wsm / ws_r
DELE      = srdel * wow**(1/2 - (l+1))
QI        = qpar  * wow**(1 - 2*(l+1))
X         = 1 - (QI-QM)*(C-ENU)/DELE**2
Y         = (QI-QM)/((C-ENU)*(QI-QM)-DELE**2)
center    = (C-ENU)*X + ENU + VMAD
shifted   = (C-ENU)*X
width     = DELE*X
obar      = Y
```
The notation used by R5 is therefore:

| Live quantity | R5 notation | Meaning |
|---|---|---|
| `potential%c + vmad` | `C` | native orthogonal-gamma center |
| `potential%dele` after `predls` | `Delta^(1/2)` | native endpoint width factor `DELE` |
| native `potential%qpar` | `gamma` | stored native-gamma screening input |
| `potential%center_band - vmad` | `bar C` | transformed center |
| `potential%shifted_band` | `bar C - E_nu` | onsite part used by live `bar h` |
| `potential%width_band` | `bar Delta^(1/2)` | transformed width |
| `potential%obar` | `bar o` | transformed overlap coefficient |
| production `lattice%sbar` | `bar S = S^alpha` | target-screened structure constants |
| `center_band - shifted_band - vmad` | `E_nu` | live onsite center correction |

The fixture reconstructs

\[
\bar h = \bar C-E_\nu + \bar\Delta^{1/2}\bar S\bar\Delta^{1/2}
\]

as `hbar`.  `source/lr_kl_hessian.f90:951-972` forms `B` from the same
bond blocks and target-endpoint `obar`; the onsite path adds the shifted-band
term.  Thus the matrix products, including their site ordering, are tested
directly rather than inferred from variable names.

## Static production-H identities

The live source documents the second-order convention at
`source/lr_kl_hessian.f90:208-230`:

\[
H^{(2)} = E_\nu + B-QB.
\]

The R5 test evaluates both `Q` and `QB` explicitly:

| Check | Result |
|---|---:|
| `||B-hbar||_max` | `1.110223e-16` |
| `||Q-hbar*obar||_max` | `5.551115e-17` |
| `||QB-hbar*obar*hbar||_max` | `1.110223e-16` |
| `||H_live-H2_explicit||_max` | `0.000000e+00` |

This establishes that the R4 `B` object is the live first-order tight-binding
matrix `bar h`, and that the live `Q` is the reciprocal product
`bar h*bar o`, so the live second-order term is `bar h*bar o*bar h` in the
source multiplication order.

The independently reconstructed untruncated object is

\[
H_{\rm orth}^{\rm exact}
 = E_\nu + \bar h( I+\bar o\bar h)^{-1}.
\]

The native coefficient-space object is constructed independently as

\[
H_\gamma = C + \Delta^{1/2}S^\gamma\Delta^{1/2},
\qquad
S^\gamma=(I+S^\alpha D)^{-1}S^\alpha,
\quad D=\alpha-\gamma,
\]

using the R1/R2 `P^gamma`, `qpar`, `dele`, and native target-screening
quantities.  Both matrices are mapped to site-major, up/down spin, orbital
ordering before comparison.

The hard static result is:

| Check | Result |
|---|---:|
| `||H_exact-H_gamma||_max` | `2.427984e-05` |
| audit-only replacement using the `predls` canonical `alpha(0)` literal | `3.375502e-05` |

The residual is far above the roundoff gate.  No similarity, congruence, or
normalization was fitted to remove it.

The source-level reason for treating this as unresolved rather than silently
repairing it is visible in the two producers.  Legacy/default structure
constants use the lattice/native default alpha (`0.3485` for `l=0`), while
`predls` uses `qm_canonical` (`0.348485` for `l=0`) unless an explicit
potential screening array is authoritative.  `predls` also carries the
`wow`-scaled `QI` and `DELE` internally, whereas the certified native path
consumes the stored `qpar` and post-`predls` `dele`.  The R5 fixture does not
invent a conversion between these source-level conventions.  The canonical
alpha substitution was measured only as a diagnostic and does not close the
map.

Consequently the accepted result is **not** that the native gamma Hamiltonian
and production `H2` differ by a proven truncation alone.  The exact TB-to-
native-gamma static bridge has not closed for the live fixture.

## Resolvent bridge

For each existing fixed complex energy, the test evaluates

\[
G_\gamma=(zI-H_\gamma)^{-1},
\qquad
g^\gamma=(P^\gamma-S^\gamma)^{-1},
\]

and applies the live endpoint scaling with the inverse `dele` width factor on
both ends.  The scaling is performed in the same site-major ordering as the
Hamiltonian comparison.

| `z` | `||G_gamma - endpoint_scaled_g_gamma||_max` | `||G_exact-G_gamma||_max` |
|---|---:|---:|
| `-0.91 + 0.83i` | `2.482534e-16` | `7.911172e-06` |
| `-0.17 + 0.04i` | `1.069508e-14` | `7.835880e-05` |
| `0.62 + 0.31i` | `1.047382e-15` | `1.523696e-04` |
| maximum | `1.069508e-14` | `1.523696e-04` |

Thus the coefficient resolvent-to-`H_gamma` endpoint relation closes, while
the exact transformed production resolvent does not coincide with that native
gamma resolvent because the static Hamiltonian gate already fails.

## Production second-order truncation

The internal production truncation diagnostic is

\[
\delta H^{(2)} = H^{(2)}-H_{\rm orth}^{\rm exact}.
\]

It gives:

| Check | Result |
|---|---:|
| `||H2-H_exact||_max` | `4.070182e-01` |
| `||H2-H_exact||_F` | `1.607265e+00` |

This is a valid measure of the finite second-order expansion on the fixture.
Because `H_exact` versus `H_gamma` did not pass, it is recorded as an
internal truncation measure only and is not promoted to a proof that it
explains the native-gamma discrepancy.

## First-derivative bridge and R4 interpretation

With

\[
A=I+\bar o\bar h,
\]

the test implements the noncommuting matrix derivative

\[
T_i^{\rm exact}=E_{\nu,i}+\bar h_iA^{-1}
 -\bar hA^{-1}(\bar o_i\bar h+\bar o\bar h_i)A^{-1}.
\]

The production derivative is independently expanded as

\[
T_i^{(2)}=E_{\nu,i}+B_i-Q_iB-QB_i,
\]

and compared with the live rotation-term API and central finite differences.

| Check | Result |
|---|---:|
| maximum FD residual for `T_exact` | `5.394976e-12` |
| maximum FD residual for `T2` | `5.351788e-12` |
| `||T_live-T2_explicit||_max` | `3.469447e-18` |
| `||T2-T_exact||_max` | `1.165563e-02` |

This proves locally that the R4 live torque is the derivative of the
production `H2`, while the exact transformed derivative contains the
additional `A^{-1}` and `obar_i` terms.  It explains the *production-side*
part of the R4 mismatch.  Since the static gamma bridge is blocked, it is not
claimed as a complete native-gamma/TB explanation.

The transformation is rotation dependent: `bar h`, `bar o`, and their spinor
blocks are built from the local moments, and the exact derivative contains
`bar o_i`.  Therefore a rotation-independent rule such as `T_gamma == T_alpha`
cannot be assumed.  No fitted connection matrix was introduced.

## Complete-curvature gate

Not entered.  The R5 protocol requires the static gamma map and first-
derivative representation map to pass before deriving mixed representation
derivatives.  Accordingly this run reports no `FD mixed residual`,
`TT_gamma`, `C_gamma`, `complete_gamma`, `TT_exact`, `C_exact`, complete
covariance residual, or production-curvature covariance result.

## Required distinction from R1/R2 and R4

The R1/R2 object previously labelled `finite_H` is the coefficient-space
gamma Hamiltonian

\[
C+\Delta^{1/2}S^\gamma\Delta^{1/2},
\]

not the accepted production `H2`.  R4 remains valid: its live torque is the
derivative of `H2` and is not elementwise equal to either tested Turek
coefficient vertex.  R5 establishes the live `B/Q/E_nu` identities and the
exact-derivative diagnostic, but stops at the unresolved static
gamma-to-production representation gate.

## Verification

```text
cmake --build build -j2
ctest --test-dir build -V -R UnitDresp03tgNativeFixedZ --output-on-failure
ctest --test-dir build -R 'UnitDresp03|UnitLrKlStaticBridge' --output-on-failure
```

`source/exchange.f90` was not modified.  Its required SHA-256 is:

```text
6e9ae7da6af46367aaf18a4fef37bdd4a9a0d6f69fe8e08b589e5c5f278a70e6
```
