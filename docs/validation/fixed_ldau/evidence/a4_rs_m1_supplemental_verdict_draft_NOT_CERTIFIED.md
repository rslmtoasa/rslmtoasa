# RS Gate M1 supplemental verdict draft — NOT CERTIFIED

This is an independent conditional assessment prepared before the pending
user decision. It is not an RS Gate M1 certification and does not supersede
the preregistered v1 result. No new measurement or rerun is proposed.

The governing v1 pair report is a failure because its extra fit-point
criterion did not pass. If, and only if, the user approves the prospective
criterion clarification recorded in
`a4_rs_m1_prereg_governance_proposal.md`, the existing measurements would
support a supplemental RS Gate M1 pass under the protocol's literal carrier
decision and the other registered checks:

* The three RS bitwise identity pairs and the local U=2 eV, J=0 capture with
  the required R05 streak passed. The local capture converged at iteration
  103 with final charge difference `4.8e-9` under `3e-8`.
* The independent LIVE anchor reproduced PL, SUMEC, SUMEV, EV, and all ten
  ATOMSC input-density components exactly; rho maximum difference was
  `7.11e-15` under `1e-10`.
* In the frozen P_OUT pair, `Delta SUMEV (spectral-local)` was
  `-0.2590196642435809 Ry`, versus `sum(t) = -0.26088126108452275 Ry`.
  The literal carrier error `0.001861596840941837 Ry` is below
  `A_RS = 0.0019179425667351765 Ry`, even before the `2e-8 Ry` numerical guard.
* Both local and spectral q2 values are nonnegative: minima
  `0.01075725460260394 Ry^2` and `0.009745244577861018 Ry^2`.
* The fit-point prediction error is `1.0234967739219947e-5 Ry`, exceeding the
  separately locked `2.0002e-8 Ry` criterion. This remains recorded as the v1
  failure and is not relabeled as a pass.

If approval is granted, the only proposed supplemental determination is that
the literal carrier and q2 gates, plus the already-passed identity, local-SCF,
and LIVE-anchor checks, satisfy RS Gate M1. The report must retain the v1
failure and explain that the fit-point comparison was an over-tight extra
criterion whose uncertainty did not include the A2 fit approximation. Until
then, RS Gate M1 remains unresolved/stopped, and RS M2 and KS lambda responses
remain unauthorized.
