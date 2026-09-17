# DRESP-03G — finite-H complex-contour Green-function bridge

Status: **finite-H spectral/contour bridge certified**.

This bridge compares the same finite-H force-theorem observable in two
representations:

* A: the certified finite-temperature spectral/eigenstate expression in
  `lr_kl_hessian_mod`;
* B: direct complex-energy resolvents of the live finite-H matrices.

It does not implement the native LMTO/Turek path operator.  In particular it
does not use `P(E)`, `g=[P-S]^{-1}`, `d(E)`, auxiliary LMTO `g`, screening
transformations, `d_matrix(E)`, or `source/exchange.f90`.

## 1. Exact derivation

For the fixed chemical potential used by DRESP-03QM,

```text
Omega(H) = -kT Tr log(1 + exp(-(H-mu)/kT)),
f(z)     = 1/(exp((z-mu)/kT) + 1).
```

For a matrix perturbation `T_i=dH/dx_i`, differentiation of the matrix
logarithm gives

```text
Omega_i = Tr[f(H) T_i].
```

With `G(z)=(zI-H)^(-1)`, the resolvent derivative is

```text
D G(z)[T] = G(z) T G(z),
D f(H)[T] = (1/(2*pi*i)) integral_C f(z) G(z) T G(z) dz.
```

For a finite-q perturbation, the certified convention is

```text
T_q(k):       rows at k+q, columns at k,
T_minus_q(k): rows at k,   columns at k+q,
C_ab(k;q,-q): rows and columns at k.
```

The exact ordered-pair Hessian implemented by the new module is therefore

```text
Omega_TT_ab(q,-q) = (1/Nk) sum_k (1/(2*pi*i)) integral_C dz f(z) *
  1/2 { Tr[G_(k+q)(z) T_a(q;k) G_k(z) T_b(-q;k+q)]
      + Tr[G_(k+q)(z) T_b(q;k) G_k(z) T_a(-q;k+q)] },

Omega_C_ab(q,-q) = (1/Nk) sum_k (1/(2*pi*i)) integral_C dz f(z) *
  Tr[G_k(z) C_ab(k;q,-q)],

Omega_ab = Omega_TT_ab + Omega_C_ab.
```

The first trace is written with the endpoint resolvent on the left because
`T_q` maps source columns to endpoint rows.  The second trace is the explicit
`a <-> b` ordering required by the certified symmetric ordered-pair form.  No
additional factor of two is present.  The `1/2` remains because both orderings
and all `(n,m)` pairs are retained.  The contact term is not symmetrized or
dropped; it is the source-space trace shown above.  The BZ normalization is
exactly `sum_k w_k / sum_k w_k`, matching DRESP-03QM.

Inserting the Lehmann identity only as a proof, not as an implementation,

```text
G_k(z) = sum_n |n k><n k|/(z-e_nk),
```

gives

```text
(1/(2*pi*i)) integral_C f(z)
       /[(z-e_m,k+q)(z-e_n,k)] dz
 = (f(e_n,k)-f(e_m,k+q))/(e_n,k-e_m,k+q)
 = K_nm.
```

For coincident energies the residue is the analytic limit
`K_nn=f'(e_n)=-f_n(1-f_n)/kT`.  This proves the exact connection to the
metallic DRESP-03QM kernel without evaluating a small denominator in the GF
route.

## 2. Legal finite-temperature contour

The Fermi function is meromorphic, with poles

```text
z_l = mu + i*pi*(2*l+1)*kT,       Res[f,z_l] = -kT.
```

A naive zero-temperature occupied contour is therefore not the finite-T
observable.  The implementation constructs a counter-clockwise ellipse that
encloses the real spectra of both `H(k)` and `H(k+q)`.  The spectral bounds are
obtained from independent Hermitian Gershgorin bounds; no eigensolve is used
to make the contour.

The ellipse may enclose Fermi poles.  Their locations are enumerated exactly.
For stable quadrature, the contour integrand uses

```text
f_reg(z) = f(z) + kT * sum_(z_l inside C) 1/(z-z_l),
```

which cancels the enclosed Fermi poles.  The direct-resolvent pole residues
are then added with coefficient `+kT`:

```text
physical response = integral_C f_reg(z) R(z) dz/(2*pi*i)
                   + kT sum_(z_l inside C) R(z_l).
```

Here `R` is either the symmetrized TT trace or the contact trace.  This is an
exact residue identity: the added rational terms have zero complete contour
integral when the physical poles and the enclosed Fermi pole are both inside
the same contour, while the explicit `+kT R(z_l)` converts the meromorphic
contour into the physical-spectrum contour.  The implementation evaluates
the pole terms with the same direct LU solve as the contour nodes.

If `contour_account_fermi_poles=.false.`, the height is restricted below the
first Fermi pole and no residue correction is used.  This is a useful audit
mode but is much less practical at low temperature.  The production default
accounts for poles explicitly.

The controls are deliberately small:

```fortran
finite_h_response_backend = 'spectral' | 'contour' | 'both'
contour_points            = 32
contour_shape             = 'ellipse'
contour_margin            = 0.25
contour_height_fraction   = 0.35
contour_account_fermi_poles = .true.
```

The default response backend remains `spectral`.  `both` writes the spectral
and contour TT/contact/total columns side by side and includes their absolute
total residual.

## 3. Zero-temperature limit

For a finite gapped fixture, choose a contour enclosing all occupied poles
and no unoccupied pole, set `f(z)=1`, and take `T -> 0`.  The bridge evaluates

```text
Omega_TT = (1/(2*pi*i)) integral_C R_TT(z) dz,
Omega_C  = (1/(2*pi*i)) integral_C Tr[G_k(z) C] dz.
```

The residue coefficient is one for an occupied physical pole and zero for an
unoccupied pole.  This is the classical occupied-state contour formulation.
It is a limit of the fixed-`mu` finite-T observable, not the same numerical
algorithm: finite-T uses `f(z)` and Fermi-pole accounting, while the legacy
route uses an occupied projector/step contour.  The zero-T API receives
explicit occupied bounds strictly inside the gap and does not diagonalize to
discover them.

## 4. Direct-resolvent implementation

`source/lr_kl_contour.f90` is a focused matrix-only module.  For every contour
node and every `(k,q)` pair it forms

```text
A_k(z)   = z I - H(k),
A_kq(z)  = z I - H(k+q),
```

and factors both matrices with LAPACK `zgetrf`.  Vertex and contact products
are evaluated by `zgetrs` solves, reusing the LU factors at that energy.  The
code never receives or imports eigenvalues/eigenvectors and contains no
Lehmann reconstruction.  In the production adapter, both endpoint matrices
are assembled directly from the copied live finite-H fixture with the same
Hamiltonian builder used by the accepted finite-H representation.

The T/C matrices are assembled once per `(k,q)` using the existing
`assemble_lmto_finite_q_torques` and
`assemble_lmto_finite_q_mixed_derivative` routines, then reused at every
contour energy.  No T/C assembly occurs inside the energy loop.

## 5. Independent finite-matrix oracle

The noncommuting 4x4 fixture has complex Hermitian off-diagonal structure,
noncommuting `H`, `T_i`, `T_j`, and `C_ij`.  It evaluates:

* A — direct four-point finite differences of the diagonalized grand
  potential;
* B — the existing finite-T spectral DRESP-03QM formula;
* C — direct complex-resolvent contour solves.

Observed residuals from `UnitDresp03gContour`:

| comparison | residual |
|---|---:|
| B minus C, TT/contact/total max | `3.47e-14` |
| A minus C, complete Hessian | `1.04e-08` |
| gapped zero-T contour vs legacy occupied spectral | `1.56e-17` |
| exact-degeneracy unitary-similarity residual | `2.78e-17` |

The finite-difference residual is consistent with the deliberately finite
`2e-4` difference step.  The fixture includes an exact degeneracy away from
the Fermi edge, near-degenerate structure, and levels close to the finite-T
Fermi crossover.  The GF route has no small-denominator branch.

## 6. Contour-node convergence

For the same noncommuting finite-T fixture, with the same ellipse and explicit
Fermi-pole handling:

| contour nodes | TT residual | contact residual | total residual |
|---:|---:|---:|---:|
| 8  | `5.26e-03` | `1.35e-03` | `6.61e-03` |
| 16 | `1.63e-04` | `7.89e-04` | `9.51e-04` |
| 32 | `6.99e-06` | `6.06e-07` | `7.60e-06` |
| 64 | `1.22e-08` | `1.29e-08` | `6.42e-10` |

Convergence is judged on TT and contact independently, not only on their
cancellation.  The gapped zero-T check was converged with 128 nodes.

## 7. Production dispatch and provenance

The user-facing dispatch is in `source/exchange_q.f90`.  The default remains
the certified metallic spectral backend.  `contour` selects B, and `both`
evaluates A/B with the same live fixture, q path, k points, k weights, Fermi
level, temperature, torque convention, and contact matrices.

Contour output records the response backend, temperature, fixed Fermi level,
ellipse type, node count, margin, height fraction, Fermi-pole setting, k mesh,
endpoint mode, and q commensurability.  Timing records separate endpoint
work, T/C assembly, direct Hamiltonian assembly, spectral contraction, total
GF response, direct LU/solve time, contour-node integration, and total
`exchange_q` response time.

## 8. bcc-Fe material gate

The reproducible comparison deck is
`example/exchange_q/bccFe/input_dresp03g_both_24.nml`.  It uses
`native_crosscheck=.false.`, `finite_h_spectral_mode='metallic'`, the same
fixed-`mu` 300 K state, a 24³ mesh, and the six-point path

```text
xi = 0, 1/24, 2/24, 3/24, 1/4, 1/2.
```

The generated `both` table is the material closure record.  It must report
spectral and contour TT, contact, and total separately at every q, including
the explicit q=0 cancellation.  An optional 48³ run is not a gate.

The six-point 24³ run was completed with the same fixed-`mu` 300 K state on
both paths.  The production record used the ellipse with 32 nodes and
explicit Fermi-pole residues; the final column is the absolute total A/B
residual in Ry:

| xi | spectral TT | spectral C | spectral total | contour TT | contour C | contour total | abs residual | relative residual |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | -2.0429057593e-1 | 2.0429057593e-1 | 5.5605e-17 | -2.0444517058e-1 | 2.0444517058e-1 | 4.8577e-17 | 7.0278e-18 | cancellation-dominated |
| 1/24 | -2.0396429877e-1 | 2.0436348981e-1 | 3.9919103e-4 | -2.0403675686e-1 | 2.0443868274e-1 | 4.0192588e-4 | 2.7348e-6 | 6.85096e-3 |
| 2/24 | -2.0310855251e-1 | 2.0457726249e-1 | 1.4687100e-3 | -2.0310097603e-1 | 2.0458002950e-1 | 1.4790535e-3 | 1.0343e-5 | 7.04257e-3 |
| 3/24 | -2.0195246799e-1 | 2.0491732571e-1 | 2.9648577e-3 | -2.0189964370e-1 | 2.0488628168e-1 | 2.9866380e-3 | 2.1780e-5 | 7.34614e-3 |
| 1/4 | -1.9796318930e-1 | 2.0643043354e-1 | 8.4672442e-3 | -1.9786968933e-1 | 2.0639799016e-1 | 8.5283008e-3 | 6.1057e-5 | 7.21092e-3 |
| 1/2 | -1.8649911230e-1 | 2.0857029116e-1 | 2.2071179e-2 | -1.8604554103e-1 | 2.0819152372e-1 | 2.2145983e-2 | 7.4804e-5 | 3.38921e-3 |

The compact path residuals are quadrature errors of the deliberately
low-cost 32-node material record, not a tuned change of sign or prefactor.
At the common q=1/4 point, increasing the direct contour from 64 to 128 to
256 nodes gave total residuals `1.08149e-5`, `2.44279e-6`, and `2.00986e-6`
Ry; TT and contact residuals decreased independently at each refinement.
The 256-node result was `spectral=(TT,C,total)=(-1.9796318968e-1,
2.0643043394e-1,8.4672442597e-3)` and
`contour=(-1.9796445101e-1,2.0642968541e-1,8.4652343987e-3)` Ry.
Its relative total residual was `2.374e-4`; at Γ relative residuals are not
interpretable because the physical total is a cancellation at roundoff.

The previously certified spectral evidence remains the reference context:
24³ and 48³ produce smooth positive small-q curves, with the 48³ first few
`DeltaJ/q²` values approximately `0.02582`, `0.02511`, `0.02440`, and
`0.02441 Ry Angstrom²`; TT and contact remain near `-0.20` and `+0.20 Ry`.

For the k-mesh sensitivity question, the required common-q comparison is
12³ versus 24³ (and 48³ if practical), with spectral and contour results
reported at the same physical q.  Agreement at each mesh supports the
interpretation that the sensitivity is intrinsic BZ/Fermi-surface sampling,
not a Lehmann-specific defect.  A contour-only improvement is claimed only
if the independently converged mesh sequence demonstrates it.

At q=1/4, the 12³ `both` run gave
`spectral=(−1.9300452350e-1,2.0076474357e-1,7.7602200713e-3)` and
`contour=(−1.9292462318e-1,2.0068445393e-1,7.7598307487e-3)` Ry, with
total residual `3.89323e-7 Ry`.  The 24³ 256-node audit gave
`spectral total=8.4672442597e-3 Ry` and `contour total=8.4652343987e-3 Ry`,
with residual `2.00986e-6 Ry`.  The spectral mesh change is about 9.11%; the
contour follows it (about 9.08% between the converged totals).  This supports
Case A: the GF representation does not materially improve k-mesh convergence;
the sensitivity is intrinsic finite-BZ/Fermi-surface sampling.

The material timing records were 24³/six-q/32-node:
Hamiltonian assembly `15.019 s`, T/C assembly `209.262 s`, spectral
contraction `2.291 s`, direct GF total `306.812 s` (LU/solve `300.309 s`,
contour integration `83.411 s`), total exchange-q response `533.993 s`.
The separate 24³/q=1/4/256-node audit took `293.444 s` for GF total inside
`369.320 s` total response.  The contour route is therefore an independent
oracle, not a production speedup claim.

## 9. Relationship to the later native bridge

This module proves only A <-> B for orthogonal finite-H matrices.  It leaves
the native LMTO/Turek representation intentionally separate.  The next bridge
may introduce the path operator and its energy-dependent transformations only
after this finite-H contour result is independently closed.  DRESP-04 is not
started here.

## Final verdict

**PASS — FINITE-H SPECTRAL/CONTOUR GF BRIDGE CERTIFIED**
