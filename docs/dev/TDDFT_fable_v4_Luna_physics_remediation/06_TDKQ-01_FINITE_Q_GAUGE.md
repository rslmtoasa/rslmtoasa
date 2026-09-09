# TDKQ-01 — Close the finite-q folded/unfolded endpoint gauge contract

## Prerequisites

TDCOV-01 through TDCOV-03 and TDWARD-01/02 must be green.

## Goal

Prove that the production finite-q TD-DFT path treats `k+q` folding consistently with the LMTO site/basis Fourier gauge, especially for multisite cells.

The code already contains folded/unfolded helper concepts. This task must establish whether production actually applies the required transformation and patch it only according to the derived Fourier convention.

## Closure record

### Fourier derivation (verified against production assembly)

The reciprocal assembler uses the actual directed bond displacement

\[
d_{ij}=R+\tau_j-\tau_i
\]

and the convention

\[
H_{ij}(k)=\sum_R H_{ij}(d_{ij})\exp(+i2\pi k\cdot d_{ij}).
\]

For \(k+q=\bar k+G\), where \(\bar k=[k+q]_{\rm BZ}\), the integer reciprocal
translation gives

\[
H(\bar k+G)=U_G^\dagger H(\bar k)U_G,\qquad
(U_G)_i=\exp(+i2\pi G\cdot\tau_i).
\]

Consequently, the eigensolver's folded target spinor is converted to the
explicitly unfolded endpoint by

\[
u_U=U_G^\dagger u_F,qquad
(u_U)_i=\exp(-i2\pi G\cdot\tau_i)(u_F)_i.
\]

The pair endpoint phase is independently verified as
\(\exp(+i2\pi(k+q)\cdot d_{ij})\) in `lmto_endpoint_phases`; its sign is not
changed for this remediation.

### Production representation and audit result

Production selects Route U. `kpoint_workset%shifted()` folds each arbitrary
`k+q` point and copies weights only from a complete-BZ workset. The arbitrary
endpoint eigensolve therefore returns \(u_F\), and
`unfold_kq_eigenvectors_for_response` applies the site-block
\(\exp(-i2\pi G\cdot\tau_i)\) transform exactly once before the arrays reach
the chi0 and Xi transition accumulators. The pair-potential source retains its
direct `+i` endpoint phases and is not gauge-transformed a second time. A
Route-F row-gauged operator remains available only for the equivalence oracle.
The production helper converts the stored Cartesian-in-`alat` primitive-site
coordinates through the inverse real-lattice matrix before applying this phase,
so its \(\tau_i\) are the same direct coordinates used by `ham_vec_type_direct`.

Finite-q response setup explicitly disables reciprocal symmetry reduction and
time reversal, with a rank-zero warning when reduction was requested. The
shifted-workset API rejects a reduced workset with a fatal diagnostic rather
than reusing irreducible weights after the endpoint little group changes. The
workset constructor requires its `complete_bz` declaration, so this policy
cannot be bypassed by omitting a flag at a new call site.

The output metadata records the requested fractional `q`, Cartesian `q`, the
Fourier phase convention, whether any endpoint folded, the first observed
integer `G` shift, and the number of folded endpoints.

---

## Required derivation

From the actual reciprocal Fourier convention used by the LMTO Hamiltonian and pair-potential vertices, derive the transformation under

\[
k+q = [k+q]_{\rm BZ} + G.
\]

For a multisite basis at positions tau_i, determine the exact site-dependent phase, schematically

\[
u_i(k+q)
=
e^{\pm i G\cdot\tau_i}
u_i([k+q]_{\rm BZ}),
\]

with the sign fixed from the code convention.

Do not copy the sign from a helper name or comment without verifying it against the Fourier transform.

---

## Production audit

Trace:
- `kpoint_workset%shifted()`;
- arbitrary k+q eigensolve;
- pair-potential endpoint construction;
- any `lmto_unfold_site_spinors()` helper;
- transition accumulation.

Determine whether:
1. the target eigenvectors are transformed;
2. the vertex is transformed equivalently;
3. or no transform is needed under the exact basis convention.

Only one representation should be chosen in production. Avoid double-applying the gauge.

---

## Decisive two-site BZ-crossing oracle

Construct a deterministic two-site cell with nontrivial basis positions.

Choose q and k points such that at least one target endpoint crosses the chosen reciprocal-cell boundary.

Compute the same transition in two mathematically equivalent ways:

### Route U — explicitly unfolded
Use the unfolded k+q representation consistently.

### Route F — folded plus gauge
Use `[k+q]_BZ` eigenvectors plus the derived reciprocal-translation site gauge.

Compare:
1. transformed eigenvectors/site coefficients as applicable;
2. pair-potential transition amplitudes;
3. chi0 at selected omega;
4. Xi.

Require agreement at tight numerical tolerance.

The test must fail if the required site gauge is intentionally omitted.

---

## Shifted workset / symmetry weights

Audit the validity of copying weights from a reduced irreducible workset to its shifted k+q set.

If finite q requires full-BZ sampling under the current implementation:
- enforce that contract;
- reject invalid reduced shifted worksets;
- clearly log when a user request for symmetry reduction is overridden for TD-DFT.

Do not silently reuse irreducible weights when the little group changes.

---

## q-sign and phase metadata

Ensure the output records:
- input q coordinates;
- Cartesian q;
- Fourier phase convention;
- whether any k+q endpoint was folded;
- reciprocal G shift if useful for debug/validation mode.

Do not change q sign to make a spectrum look physically convenient.

---

## Forbidden shortcuts

- No monoatomic-only test: site gauge is trivial there.
- No Gamma-only test.
- No assumed sign for exp(±i G·tau).
- No copying reduced weights without proof.
- No simultaneous change to mode extraction/convergence tuning.

---

## Acceptance checklist

- [x] Exact folded/unfolded gauge derived from production Fourier convention.
- [x] Production k+q path traced end to end.
- [x] Two-site BZ-crossing fixture added.
- [x] Explicitly unfolded and folded+gauge transition amplitudes agree.
- [x] chi0 agrees between the two representations.
- [x] Xi agrees between the two representations.
- [x] Test fails when required gauge is intentionally omitted.
- [x] Shifted reduced-workset policy is rigorous and explicit.
- [x] q/Fourier metadata is sufficient to debug the path.
- [x] Existing finite-q tests pass.

### Validation record

The deterministic two-site crossing reports transition, chi0, Xi, and omitted-
gauge errors of `3.4694e-18`, `0.0000e+00`, `8.6736e-19`, and `8.1553e-02`,
respectively. The first three satisfy the `1e-10` oracle tolerance; the
omitted-gauge control exceeds `1e-6` as required. The focused finite-q CTest
command was:

```text
ctest --test-dir build --output-on-failure -R 'UnitLmtoPairPotential|UnitKpointWorkset|UnitKpointWorksetReducedReject|UnitTddft(ChiKS|DirectXi|CircularCovariance)|UnitTddftDispatch'
```

It completed with 7/7 tests passed.

## Required one-line commit message

`td-dft: close finite-q endpoint gauge contract`
