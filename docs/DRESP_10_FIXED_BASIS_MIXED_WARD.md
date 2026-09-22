# DRESP-10F: fixed-basis mixed-representation static Ward gate

DRESP-10F is the fixed-ground-state static response formulation

```text
chi0^(P<-SR) = R_P L_f(H) F_SR
```

The input field is a direct physical scalar-relativistic (SR) source.  The
response is measured in the live Pauli product basis.  These are intentionally
different coordinate spaces; the primary calculation does not project the
field into the Pauli product basis and does not add a moving-basis delta-O or
contact term.

## Primary fixed-basis map

```text
Pauli density U_P
       |
       v
raw m_P = P3(U_P)
       |
       v
Kxc = Bxc / P3  --------------------+
       |                              |
       v                              |
raw SR field F_SR = Kxc * m_P <-------+
       |
       v
delta H = F_SR(source vertex)
       |
       v
delta rho = L_f(H)[delta H]
       |
       v
Pauli measurement R_P (348 branches)
       |
       v
U_P(response)
```

The product basis is the six-branch LMTO basis `00, 10, 01, 11, 20, 02` and
has dimension 348 for the certified bcc-Fe input.  The primary action is
matrix-free.  Its acceptance gates are, in order: corrected six-branch span
audit, direct raw-SR source vertex, exact static Fréchet action, mixed
`chi0^(P<-SR)` action/oracles, and the direct `Kxc=P3` identity.

The radial span diagnostic is independent of the production branch routine.
Its endpoint expansion applies the `Enu` shift exactly once.  The old Pauli
residual `7.801e-1` was therefore a diagnostic defect, not a license to add a
contact term.  Historical values retained for comparison are the superseded
SR-total `7.920e-1` and delta-O `5.990e-1`.

## Separate covariance oracle

```text
global spin rotation
        |
        +--> production H representation tangent
        |          |
        |          v
        |     native observable/basis tangent
        |          |
        |          v
        |     exact rotated P3
        |
        +--> DRESP-09U representation artifact
        +--> DRESP-09Y augmentation artifact
```

This is a separate covariance diagnostic, not an input to the primary fixed
basis action.  It reports the native tangent together with the DRESP-09Y
augmentation contribution and reruns DRESP-09U and DRESP-09Y as sidecar
artifacts.  A closed fixed-basis action and closed covariance oracle are
reported as `PASS-A FIXED_BASIS_MIXED_WARD_CLOSED`; finite LMTO covariance
disagreement after all hard gates is `PASS-B FINITE_LMTO_WARD_INCONSISTENCY`.
Unclosed span, source, Fréchet, or mixed-action gates are `BLOCKED`.

BES/Halle and Goldstone correction are off for this gate.  Dynamics are not
run.  The implementation is intended as a static, auditable bridge between
the fixed-ground-state Fréchet response and the LMTO endpoint representation;
further physics extensions should follow only after the gates and the
sidecar covariance oracle are reviewed.
