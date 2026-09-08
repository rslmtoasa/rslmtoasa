# TDCOV-02 — Separate pair moment amplitude from signed Goldstone magnetization

## Prerequisite

TDCOV-01 must already be committed and green. Do not start from a tree where circular channels may use different Fermi levels.

## Goal

Repair the global-spin-reversal failure of the LMTO pair-potential route without introducing an ad hoc sign convention.

The key physical distinction is:

\[
\mathbf m_i = M_i \hat{\mathbf e}_i,\qquad M_i=|\mathbf m_i|>0.
\]

For a transverse orientation derivative at fixed moment magnitude,

\[
\frac{\partial H}{\partial m_{ix}}
=
\frac{1}{M_i}
\frac{\partial H}{\partial e_{ix}},
\qquad
\frac{\partial H}{\partial m_{iy}}
=
\frac{1}{M_i}
\frac{\partial H}{\partial e_{iy}}.
\]

But a Goldstone displacement of a collinear reference state may depend on signed

\[
m_{iz}=\pm M_i.
\]

These two roles must not share one ambiguously named variable.

---

## Reproduce the old failure first

Using a real LMTO-backed one-site fixture with otherwise identical physical parameters:

1. reference moment genuinely `+z`;
2. reference moment genuinely `-z`.

Run the pair-potential static Goldstone route.

Record the pre-patch closest Goldstone eigenvalues. The expected old pathology is approximately:
- one orientation -> `+1`;
- reversed orientation -> `-1`.

Do not simulate reversal by changing only a scalar `signed_moment`; reverse the underlying LMTO orientation itself.

---

## Required implementation

1. Trace `lmto_magnetic_tangent` and document exactly what variable it differentiates.
2. Replace the pair-potential normalization API semantics:
   - remove/rename `signed_moment` where it is acting as a Jacobian denominator;
   - pass a **positive `moment_amplitude`** instead.
3. Enforce `moment_amplitude > tolerance`.
   - Do not hide caller mistakes by applying `abs()` deep inside every consumer.
4. Retain signed `m_z` only where the physics requires a signed Goldstone displacement/vector.
5. Audit all callers of the pair-potential builder for the same semantic confusion.
6. Preserve the current circular operator normalization unless the derivation independently proves it wrong.

---

## Mandatory physical oracles

### Oracle A — true global spin reversal
For the same one-site collinear LMTO system:

\[
+\hat z:\quad \Xi(0,0)m = m
\]

and

\[
-\hat z:\quad \Xi(0,0)m = m.
\]

Require the Goldstone eigenvalue to be `+1` for both, to tight static-solver tolerance.

### Oracle B — finite-difference tangent covariance
For both +z and -z references:
- apply small ±delta rotations about x and y;
- finite-difference the Hamiltonian;
- compare to analytic LMTO orientation tangents;
- verify the orientation-to-moment conversion uses positive amplitude.

Do not use a single one-sided finite difference if a centered difference is feasible.

### Oracle C — real two-sublattice +/-z fixture
Repair the existing multisite test so the underlying LMTO site orientations are genuinely opposite.

Use:
- positive moment amplitudes for pair normalization;
- signed site magnetizations for the Goldstone vector.

The acoustic global-rotation mode must close the Ward/Goldstone identity with eigenvalue +1.

### Oracle D — no regression under +z
The previously good Fe-like +z pair-Ward closure must remain at essentially machine/static-solver precision.

---

## Test cleanup

Any existing unit test that currently expects the pair operator to flip sign merely because a scalar `signed_moment` changes sign while the underlying LMTO orientation remains fixed must be rewritten. That expectation encodes the bug.

Keep arithmetic/unit coverage, but make its variable names reflect real semantics.

## Execution record

The production-backed regression fixture was run with the same LMTO parameters
for both underlying orientations.  The explicit pre-patch signed-denominator
replay produced `(+z,-z) = (+1.00000000,-1.00000000)`.  With the repaired
positive amplitude path, the one-site static Goldstone eigenvalues were
`(+z,-z) = (+1.00000000,+1.00000000)`, and the genuine alternating two-site
fixture returned the acoustic action `(1.00000000,-1.00000000)` for the signed
site vector `(m_z,-m_z)`.

The centered ±x/±y finite-rotation tangent oracle passed for both +z and -z;
the maximum reported Cartesian tangent error was `5.0875e-10` at the test
tolerance `2.0e-9`.  The existing +z finite-q/static pair coverage remained
green with a maximum reported error of `8.8558e-10`.

Validation commands:

```text
cmake --build build --target UnitLmtoPairPotential -j2
./build/bin/UnitLmtoPairPotential
python3 tests/unit/test_tddft_dispatch.py
```

---

## Forbidden shortcuts

- No blanket `abs(signed_moment)` patch without API cleanup and covariance tests.
- No extra minus sign in `Q+`, `Q-`, `Xi`, Dyson, or Goldstone vectors to force +1.
- No material-specific special cases for Ni.
- No tolerance relaxation.
- No modification of Fermi logic from TDCOV-01.
- No change to finite-q endpoint gauge in this task.

---

## Acceptance checklist

- [x] Old +z/-z pair-Ward asymmetry reproduced.
- [x] LMTO tangent variable documented.
- [x] Pair builder consumes positive `moment_amplitude`.
- [x] Signed `m_z` is retained only for signed Goldstone displacement/related physics.
- [x] +z one-site Goldstone eigenvalue is +1.
- [x] -z one-site Goldstone eigenvalue is +1.
- [x] Centered finite-rotation tangent oracle passes for both orientations.
- [x] Actual +/-z two-sublattice Goldstone oracle passes.
- [x] Existing +z static Ward precision is preserved.
- [x] No ad hoc compensating sign introduced elsewhere.
- [x] Existing relevant tests pass.

## Required one-line commit message

`td-dft: separate pair moment amplitude from Goldstone sign`
