# TDCOV-03 — Enforce circular-response covariance

## Prerequisites

TDCOV-01 and TDCOV-02 must be committed and green.

## Goal

Make `plus_minus` and `minus_plus` faithful implementations of the corresponding retarded circular correlators and prove their covariance at Gamma, finite q, negative q, negative omega, and global spin reversal.

Do **not** force both circular sectors to display a positive-frequency pole by redefining q or swapping k/k+q endpoints.

---

## Mathematical contract

Use the code's exact Fourier and operator conventions to derive the implemented relation. The expected structure is

\[
\chi^{+-}(q,\omega)
=
\left[\chi^{-+}(-q,-\omega)\right]^\dagger.
\]

At omega = 0:

\[
\chi^{+-}(q,0)
=
\left[\chi^{-+}(-q,0)\right]^\dagger.
\]

For inversion-symmetric collinear systems this may reduce to same-q conjugacy at omega=0, but only use that simplified oracle when the symmetry assumptions are explicitly satisfied.

The derivation must state:
- operator definitions;
- source/measurement ordering;
- Fourier convention;
- Lehmann denominator convention;
- matrix dagger/transposition in response space.

---

## Required implementation

1. Ensure both circular channels use the identical resolved occupation state from TDCOV-01.
2. Audit channel construction and metadata.
3. Remove hard-coded assumptions that "reverse" always means `minus_plus`.
4. Build filenames and metadata from the actual channel enum/object.
5. If the output layer wants a positive-energy view of opposite chirality, implement it only as an explicitly derived/relabelled output product; do not alter the retarded correlator.

---

## Mandatory covariance fixtures

### Fixture A — Gamma
SOC-off one-site collinear LMTO case:
- compute both channels at omega=0;
- require the appropriate conjugacy identity to tight numerical tolerance.

This is a necessary algebraic test but is not sufficient.

### Fixture B — finite +q and -q
Choose a nonzero q small enough to remain a simple test but large enough to exercise the arbitrary-k+q path.

Compute:
- `chi^{+-}(+q, 0)`;
- `chi^{-+}(-q, 0)`.

Verify the derived matrix identity.

This test must not use the Gamma short-circuit.

### Fixture C — +/- omega
For one finite q, use a small symmetric frequency grid containing negative and positive omega.

Verify the full retarded circular covariance relation within a justified tolerance.

The purpose is to prove that opposite circular sectors naturally carry mirror spectral weight; do not transform one channel to force its pole to the other side.

### Fixture D — global spin reversal
For otherwise identical +z and -z collinear LMTO states:
- the physical circular sectors must exchange appropriately;
- excitation energies must be invariant under global spin reversal;
- pair Goldstone closure from TDCOV-02 must remain +1.

### Fixture E — response rank > 1 where practical
Add at least one multisite response-space covariance test so matrix dagger/transposition is actually exercised and not reduced to scalar conjugation.

---

## Numerical tolerances

For algebraically identical eigenpair sums using the same states, target much tighter agreement than the old 5% dynamic backend gate.

Do not pick a tolerance first and make the test fit it. Record observed roundoff/quadrature differences and choose a tolerance with clear margin.

For deterministic same-backend covariance, expect near machine precision unless there is a documented source of numerical asymmetry.

---

## Forbidden shortcuts

- No q -> -q substitution inside production solely to make a pole positive.
- No k/k+q swap solely to make `minus_plus` look like `plus_minus`.
- No sign flip in the loss writer to enforce positivity for both sectors.
- No assumption that same-q conjugacy holds without checking spatial symmetry.
- No Gamma-only validation.

## Implementation record

The implemented response-space derivation is the following.  For a site
projector \(P_i\), the measured operators are

\[
A_i^+ = P_i O_+,\qquad A_i^- = P_i O_-,\qquad
O_\pm = \sigma_x \pm i\sigma_y,
\]

so \(O_+^\dagger=O_-\).  The transition vertex routine evaluates the left
factor as the measurement matrix element
`<n,k|A_i|m,k+q>` and the right factor as the source matrix element
`<m,k+q|B_j|n,k>`.  The eigenpair backend therefore implements

\[
\chi_{ij}^{AB}(q,\omega)=\sum_{k,n,m}w_k
\frac{f_{n,k}-f_{m,k+q}}
{\omega+\epsilon_{n,k}-\epsilon_{m,k+q}+i\eta}
\langle n,k|A_i|m,k+q\rangle
\langle m,k+q|B_j|n,k\rangle.
\]

The time convention is the retarded one associated with
`delta W(t) = Re[delta W(omega) exp(-i omega t)]`; spatial real-space output
uses `exp(-i q.R)`, and the eigenpair endpoint contract is explicitly `k+q`.
For the SOC-off collinear fixture, complex conjugation together with the
time-reversed endpoint relabelling \(k\mapsto-k\), \(q\mapsto-q\), and
\(\omega\mapsto-\omega\) gives

\[
\chi_{ij}^{+-}(q,\omega)=
\left[\chi_{ji}^{-+}(-q,-\omega)\right]^*,
\]

which is exactly
`chi_plus(q,omega) = transpose(conjg(chi_minus(-q,-omega)))` in the
response-space matrix, i.e. the required dagger relation.  At zero frequency
the same equation applies.  Same-q conjugacy is not used by the regression
unless the fixture's spatial symmetry makes it valid.

The production driver resolves one ordered channel object, obtains its
opposite through `opposite_circular_channel`, and derives the measurement and
source channels, metadata labels, peak labels, and filename tags from those
codes.  Both ordered options are copies of the single resolved
`response_occupation` state.  No q/k+q substitution or positive-frequency
relabeling is performed in production.

The dedicated `UnitTddftCircularCovariance` fixture records these same-backend
residuals:

```text
Gamma dynamic:       0.0000E+00   (tolerance 1E-12)
Gamma static:        0.0000E+00   (tolerance 1E-12)
finite q and omega:  3.5804E-15   (tolerance 1E-11)
rank-2 off diagonal: 2.1875E+01
spin reversal:       0.0000E+00   (tolerance 1E-12)
```

---

## Acceptance checklist

- [x] Circular covariance derived using actual code conventions.
- [x] Primary/reverse channels share identical occupation metadata.
- [x] Gamma conjugacy oracle passes.
- [x] Finite +q/-q oracle passes.
- [x] Symmetric +/-omega covariance oracle passes.
- [x] Global spin reversal exchanges circular sectors correctly.
- [x] At least one rank>1 response-space test exercises matrix dagger semantics.
- [x] Channel metadata and filenames reflect the actual channel object.
- [x] No production redefinition was introduced to force positive-frequency poles.
- [x] Existing relevant tests pass.

## Required one-line commit message

`td-dft: enforce circular-response covariance`
