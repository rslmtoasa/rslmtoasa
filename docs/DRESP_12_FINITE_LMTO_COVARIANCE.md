# DRESP-12 — finite-LMTO rigid-response decomposition

Starting HEAD for this closure: `e4fe233d47418a004641f91cb03c40456bf37def`.

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

## Final circular/Cartesian Pauli seam audit

The production input selects `chi_plus`. Its native convention is retained:
the independent circular paths are `m_plus` (up to down) and `m_minus` (down
to up), with

```text
m_x = m_plus + m_minus
m_y = i*(m_plus - m_minus).
```

The frozen DRESP-09Y authority is evaluated without rewriting its path:
`sr_l0_density_from_second_order_endpoint_branches` is called independently
for channels 1 and 2 with `lower_factor = 0`, followed by the exact
`convert_response` convention. DRESP-12's current
`endpoint_pauli_measurement` is evaluated on the identical
`H(k)`, `rho(k)`, `deltaH_cov(k)`, `delta_rho_cov(k)`, and six-branch endpoint
tangent. The audit compares both paths after making the DRESP-09Y conversion
explicit:

```text
weighted = sqrt(4*pi)*r^2*physical_L0_density.
```

The factor-of-two single-block expression is not assumed equivalent to the
two circular paths. The accepted Fe audit establishes the aggregate
equivalence only after the endpoint Hermitian and radial branch swaps:

```text
endpoint branch swap: 00<->00, 10<->01, 01<->10, 11<->11, 20<->02, 02<->20
Hermitian residuals: 7.15e-17, 1.09e-12, 1.09e-12, 6.98e-13, 3.97e-12, 3.97e-12
radial swap residuals: 1.65e-15, 1.74e-15, 2.81e-15, 0, 0, 0
```

The branch-resolved current DRESP-12 values differ for the swapped 10/01 and
20/02 labels; the aggregate measurement closes because those endpoint and
radial product swaps are included. The independent coefficient-space trace
oracle agrees with both aggregate paths, so this is not a remaining
measurement ambiguity.

The accepted seam values are:

```text
DRESP-09Y endpoint norm                         4.901514498244956e-1
DRESP-12 endpoint norm                          4.901514498244927e-1
12-vs-Y residual                                 5.88e-15
trace oracle vs Y                                2.67e-16
trace oracle vs DRESP-12                         5.84e-15
m_plus norm                                      2.450757249122463e-1
m_minus norm                                     2.450757249122492e-1
x reconstruction residual                        1.00e-16
y reconstruction residual                        5.80e-15
observable plus/minus -> x residual              2.08e-16
embedded DRESP-09Y Pauli_complete_vs_P3         5.30e-15
current DRESP-12 Pauli_complete_vs_P3            3.81e-15
```

The measurement seam is therefore `MEASUREMENT_CONVENTIONS_IDENTICAL`.

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

The corrected accepted-Fe result is `PASS-B RESIDUAL_BASIS_RESPONSE_REMAINS`.
The fixed-basis Goldstone residual remains `2.363016977393652e-1`, while the
embedded circular static sum rule and Cartesian rigid covariance oracle close
at approximately `5.3e-15`. The complete DRESP-12 accounting residual remains
`9.475388813621777e-1`; this is the residual basis-response result, not a
Pauli measurement seam failure. The circular plus/minus single-channel
residuals are diagnostic only; their explicit recombination is the Cartesian
closure target.

## Required static checks

The integration artifact reports the six complete branches, six frozen-H
branches, endpoint-H subtraction and `delta_rho=0` isolation, branchwise
closure, Frechet linearity, complete fixed-observable closure, endpoint
Hermitian swaps, radial spin-direction swaps, independent plus/minus paths,
Cartesian reconstruction, coefficient-trace oracle, DRESP-09Y observable
plus/minus reconstruction, embedded DRESP-09Y P3 closure, current-versus-
authoritative measurement, corrected master accounting, compact `D*m_G`,
endpoint-H Bxc/connection diagnostics, and the frozen DRESP-10F/DRESP-11
regressions. Dynamics is not run.

The artifact is `/tmp/dresp12_fe_4k.dat` with the k-resolved table and the
existing DRESP-09U/Y sidecars.

## Closure result

The accepted Fe run closes the endpoint algebra, frozen-H equivalence,
Frechet linearity, complete fixed-observable covariance identity, the
endpoint measurement seam, the independent coefficient trace, and the
DRESP-09Y observable regression. The endpoint-H term reduces the historical
incomplete accounting residual from `1.4619` to
`0.9475388813621777`, but does not close the remaining radial/profile/vector
defect. The resulting classification is
`PASS-B RESIDUAL_BASIS_RESPONSE_REMAINS`, with
`NEXT = LMTO_BASIS_RESPONSE_AUDIT`. Goldstone correction is `OFF` and
dynamics is `NOT RUN`.
