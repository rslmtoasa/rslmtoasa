# DRESP-12 — finite-LMTO rigid-response decomposition

Starting HEAD for this closure: `e86d2885b64f8d36761f2d03fc61f6143eebbaaf`.

DRESP-12 is a static, diagnostic-only Gamma-point closure for the accepted
64-k bcc-Fe `ham_only` state. It does not modify the Hamiltonian-side
covariance bridge, response kernel, denominator, SVD, or production dynamics.
Goldstone restoration, BES/Halle paths, fitting, rescaling, and correction
routes are all off.

## Complete finite-LMTO response

The six production endpoint objects are

```text
D_pq = H**p rho H**q,   pq = 00, 10, 01, 11, 20, 02.
```

Their tangent is split exactly as

```text
delta D_pq = delta D_pq^(rho) + delta D_pq^(H)
delta D_pq^(rho) = H**p delta_rho H**q
delta D_pq^(H)   = delta D_pq^(complete) - delta D_pq^(rho).
```

`endpoint_tangent_branches_second_order` is authoritative for the complete
product rule. The independent
`endpoint_fixed_h_branches_second_order` primitive supplies the frozen-H
term. The endpoint-H subtraction is independently checked by calling the
complete routine with `delta_rho = 0`; branch `00` is identically zero on the
endpoint-H side.

The density-side/Kubo response is measured with the same physical Pauli
six-branch radial and angular dual as the certified DRESP-09X/Y path. Thus

```text
delta_m_cov_complete_fixedO = delta_m_cov^rho + delta_m_endpoint-H
delta_m_cov^rho = delta_m_B + delta_m_conn.
```

The observable-frame term remains separate:

```text
delta_m_O = Tr[rho delta_O].
```

The complete rigid response accounting is therefore

```text
delta_m_cov = delta_m_B + delta_m_conn + delta_m_endpoint-H + delta_m_O.
```

The fixed-basis defect and the master accounting residual are different
quantities:

```text
r_fixed = delta_m_B - m_P3
R_account = ||r_fixed + delta_m_conn + delta_m_endpoint-H + delta_m_O|| / ||r_fixed||.
```

`fixed_basis_goldstone_relative = ||r_fixed|| / ||m_P3||` is the actual
fixed-basis Goldstone consistency residual. It must not be confused with
`covariance_accounting_relative_to_fixed_defect`, the endpoint-inclusive
accounting residual. The historical incomplete value is retained as:

```text
old covariance accounting residual = 1.4619
status = SUPERSEDED; endpoint-H contribution omitted
```

## Compact D*m_G sign

DRESP-11 stores

```text
D*m_G = m_G - A*m_G = -r_fixed.
```

After the complete accounting closes, its independent reconstruction is

```text
D*m_G = delta_m_conn + delta_m_endpoint-H + delta_m_O.
```

The report records this sign explicitly and compares the reconstruction with
the frozen DRESP-11 compact vector.

## Interpretation and scope

The physical decomposition is

```text
delta_m_LMTO = delta_m_Kubo(delta_rho)
              + delta_m_endpoint-H       (energy-moment/contact)
              + delta_m_basis/observable (deltaO).
```

The static Frechet derivative is an exact fixed-matrix oracle for the first,
Kubo/bubble term. It is not by itself the complete LMTO physical response.
The Lehmann/GF susceptibility is the production dynamical counterpart of the
Kubo term; endpoint, basis, and radial-observable response terms must be
handled separately in a frequency-dependent theory. DRESP-12 does not
generalize this endpoint-H term to finite frequency.

## Decision gate

```text
PASS-A  COMPLETE_LMTO_RIGID_RESPONSE_DECOMPOSITION_CLOSED
PASS-B  RESIDUAL_BASIS_RESPONSE_REMAINS
BLOCKED ENDPOINT_RESPONSE_REGRESSION
```

PASS-A closes the Ward-focused campaign and sets
`NEXT = LMTO_DYNAMIC_RESPONSE_FORMULATION`. PASS-B sets
`NEXT = LMTO_BASIS_RESPONSE_AUDIT` and reports the remaining radial/profile
residual. An endpoint branch or independent isolation failure is blocked and
does not authorize a further interpretation.

## Required static checks

The integration artifact reports the six complete branches, six frozen-H
branches, endpoint-H subtraction and `delta_rho=0` isolation, branchwise
closure, Frechet linearity, complete fixed-observable closure, DRESP-09Y
upper/lower-small/lower-angular observable pieces, corrected master
accounting, compact `D*m_G`, endpoint-H Bxc/connection diagnostics, and the
frozen DRESP-10F/DRESP-11 regressions. Dynamics is not run.

The artifact is `/tmp/dresp12_fe_4k.dat` with the k-resolved table and the
existing DRESP-09U/Y sidecars.

## Closure result

The accepted Fe run closes the endpoint algebra, frozen-H equivalence,
Frechet linearity, complete fixed-observable covariance identity, and the
DRESP-09Y observable regression. The endpoint-H term reduces the historical
accounting residual from `1.4619` to approximately `9.4754e-1`, but does not
close the remaining radial/profile/vector defect. The resulting classification
is `PASS-B RESIDUAL_BASIS_RESPONSE_REMAINS`, with
`NEXT = LMTO_BASIS_RESPONSE_AUDIT`. The fixed-basis Goldstone residual remains
approximately `0.2363`, and the DRESP-11 orthogonal fraction remains
approximately `0.9066`.
