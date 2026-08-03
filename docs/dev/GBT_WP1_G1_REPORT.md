# GBT WP1 / G1 report

Date: 2026-08-03. Scope: independent algebraic acceptance oracles only. No
production GBT routing, Hamiltonian assembly, or reciprocal routing is changed.

## Convention

The tests use the blueprint's physical-displacement convention:

\[
\alpha=\mathbf q_{\rm cart}\mathbin\cdot
(\mathbf R+\boldsymbol\tau_b-\boldsymbol\tau_a)_{\rm cart},
\qquad
\mathbf q_{\rm cart}=2\pi A^{-T}\mathbf q_{\rm frac}.
\]

Additional sublattice azimuths are encoded only in the endpoint reference
moments. They are not added to \(\alpha\).

## Tests

`tests/unit/test_gbt_oracles.py` is self-contained and imports no production
GBT or Pauli helper. Its dense reference route builds explicit complex
\(2\times2\) spin matrices and evaluates

\[
W_a(\mathbf m_a^0)\,[S_{ab}\otimes D(\alpha)]\,W_b(\mathbf m_b^0).
\]

It compares this to a separate local component transcription, using 24 seeded
random complex orbital matrices across two and three orbitals with unequal
endpoint species. It also checks the one-orbital Stoner
\(\mathbf k\mp\mathbf q/2\) identity, reverse directed-bond Hermiticity,
and a triclinic two-sublattice fixture that fails conspicuously if the basis
displacement or reference azimuth is added to an already physical \(\alpha\).

Run:

```bash
/Users/andersb/envs/p311/bin/python tests/unit/test_gbt_oracles.py
```

The test is registered as `UnitGbtOracles` when `RUN_UNIT_TESTS=ON`.
On this workspace, `ctest --test-dir build_13 --output-on-failure -R
'UnitGbtOracles|UnitGammaSupercell|UnitDysonEquivalence'` passed all three
tests; the latter two provide the retained G0 unit baseline check.

## Result

The executable test output records the maximum errors and enforces a strict
relative threshold of \(<10^{-12}\). The multi-sublattice fixture separately
requires a wrong double-phase construction to differ from the dense oracle by
more than \(10^{-3}\), so a passing result cannot be obtained by silently
testing an insensitive geometry.

Observed with `/Users/andersb/envs/p311/bin/python`:

| Check | Maximum relative error / separation |
| --- | ---: |
| Cartesian/fractional phase equivalence | \(1.780\times10^{-15}\) |
| Stoner \(\mathbf k\mp\mathbf q/2\) identity | \(5.551\times10^{-17}\) |
| Dense unequal-species pair oracle | \(2.505\times10^{-16}\) |
| Reverse directed-bond Hermiticity | \(3.105\times10^{-16}\) |
| Multi-sublattice correct convention | \(1.380\times10^{-16}\) |
| Deliberate double phase separation | \(2.421\times10^{-2}\) |

The maximum algebraic relative error is \(1.780\times10^{-15}<10^{-12}\).

**G1: PASS.**
