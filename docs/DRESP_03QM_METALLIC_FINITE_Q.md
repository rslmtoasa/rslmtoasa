# DRESP-03QM — Metallic finite-q finite-H response

Status: **PASS — METALLIC FORMULATION CERTIFIED, PERFORMANCE OPEN**.

The finite-H route is now a fixed-chemical-potential grand-potential Hessian
with a finite-temperature occupation-difference kernel.  It has independent
gapped, finite-temperature, degenerate-subspace, finite-metal, and
commensurate-endpoint tests.  The native LKAG q oracle remains blocked as
recorded in [`DRESP_03Q_FINITE_Q_LKAG_CLOSURE.md`](DRESP_03Q_FINITE_Q_LKAG_CLOSURE.md).

`source/exchange.f90` was not modified.

## Spectral derivation

The force-theorem object is the grand potential at the accepted SCF chemical
potential `mu`:

```text
Omega(H) = -kT Tr log(1 + exp(-(H-mu)/kT)).
```

For two rotation coordinates `x_i` and `x_j`, write
`T_i = dH/dx_i` and `C_ij = d2H/(dx_i dx_j)`.  First differentiation gives
`Omega_i = Tr[f(H) T_i]`.  Differentiating the eigenprojector, or equivalently
using the first-order eigenvector response, gives

```text
Omega_ij = sum_n f_n <n|C_ij|n>
         + sum_(n,m; m/=n) f_n/(e_n-e_m)
             [ <n|T_i|m><m|T_j|n>
             + <n|T_j|m><m|T_i|n> ].
```

Reindexing the second term and retaining all ordered band pairs gives the
implemented form

```text
Omega_ij(TT) = 1/2 sum_(n,m) K_nm
              [ <n|T_i|m><m|T_j|n>
              + <n|T_j|m><m|T_i|n> ],

K_nm = (f(e_n)-f(e_m))/(e_n-e_m).
```

For finite `q`, `n` belongs to `(k,n)`, `m` belongs to `(k+q,m)`, and the two
vertices are the certified `T_q(k)` and `T_-q(k+q)`.  The Brillouin-zone
average is

```text
sum_k w_k / sum_k w_k
```

with no additional band or spin factor.  The one-half is required because all
ordered `(n,m)` pairs are present.  The contact term remains

```text
Omega_ij(C) = sum_n f(e_n) <n|C_ij|n>,
```

using the same `mu`, `kT`, and k weights.  It is not optional.  Thus the
complete Hessian is `TT + contact`, preserving the certified mixed derivative
and HOH product-rule machinery.

The old occupied-only expression is mathematically equivalent to this formula
for an exactly zero-temperature, discrete, gapped spectrum (and, formally,
for a nondegenerate zero-temperature metallic spectrum away from the Fermi
edge).  It is not the finite-temperature response of the accepted reciprocal
state, and its individual occupied/occupied terms suffer large cancellation
when metallic levels approach one another.  The new production default is
`finite_h_spectral_mode='metallic'`; `legacy_occupied` remains available as a
regression path.

## Coincident energies and intraband response

The code evaluates the smooth divided-difference limit rather than skipping a
small denominator:

```text
lim_(e_m -> e_n) K_nm = f'(e_n)
                       = -f_n(1-f_n)/kT.
```

The same positive `kT = max(T*kB, 1e-10 Ry)` convention used by reciprocal SCF
occupations is used here.  A midpoint Taylor continuation is used only to
avoid loss of digits in a numerically coincident pair; the response term is
retained.  A zero-temperature degeneracy at the Fermi edge is nonanalytic,
not a denominator that can be physically regularized by dropping it.

Because the kernel depends only on the two eigenvalues and is identical in a
degenerate subspace, the complete ordered-pair sum is invariant under an
arbitrary unitary rotation within that subspace.  The unit test measures a
`1.39e-17` residual.

## Thermodynamic ensemble

The declared observable is the grand-potential Hessian at fixed chemical
potential `mu`, with `mu` and the electronic temperature inherited from the
accepted reciprocal SCF state.  This is the force-theorem object used by this
route; no heuristic fixed-number term is added.

If a canonical fixed-electron-number free-energy Hessian is required instead,
the exact Legendre-transform correction would be

```text
F_ij|N = Omega_ij - Omega_i,mu Omega_mu,j / Omega_mu,mu.
```

That is a different observable and is intentionally outside this DRESP-03QM
closure.  The contact term and TT term in this implementation are evaluated
in the same fixed-`mu` ensemble.

## q=0 and Goldstone behavior

The one-site production path no longer replaces Gamma by a hard-coded zero.
Gamma is evaluated with the same metallic kernel and contact term as finite q.
For the final 12³ and 24³ bcc-Fe runs the raw Gamma totals were `2.81e-17` Ry
and `2.67e-18` Ry, respectively, while TT and contact remained visible at about
`-0.20` and `+0.20` Ry.  This is the numerical Goldstone cancellation.  The
multi-sublattice route continues to return the complete Gamma matrix.

## Commensurate endpoint reuse

For a full uniform mesh, the detector constructs a coordinate-index table from
the stored folded k points.  It accepts a nonzero q only when each component
is an integer mesh translation and every folded `k+q` lands on an existing
point with the same weight.  It does not round an arbitrary q.

For accepted reciprocal Hamiltonians, arbitrary-k assembly folds its input by
the production convention `k -> k-floor(k+1/2)`.  Consequently the accepted
eigenvectors at the mapped folded point are exactly the eigenvectors used by
the explicit endpoint solver; no additional basis-position gauge transform is
needed.  The finite-q torque still receives the unfurled endpoint coordinate,
so its certified endpoint phase is retained.  At runtime the first nonzero q
also checks

```text
T_-q(k+q) = T_q(k)^dagger
```

against an independent assembly.  The unit test covers q translations
`(0,0,1/N)`, `(0,0,2/N)`, and `(1/N,1/N,0)`, weight preservation, exact map
closure, and rejection of an off-mesh q.

The output header records `mesh_reuse` versus
`explicit_diagonalization`, the commensurability flag, and the coordinate
residual for every q.

## Independent validation

The `UnitExchangeQ` target reports:

| check | residual |
|---|---:|
| gapped metallic reduction | `0.0` |
| finite-temperature two-level spectral oracle | `0.0` |
| finite-temperature finite-difference grand potential | `1.53e-8` |
| degenerate-subspace rotation | `1.39e-17` |
| commensurate endpoint map/reuse oracle | `0.0` |

The finite-temperature oracle uses a two-level analytic spectrum and an
independent grand-potential evaluation, not the production kernel.  The
fixture also checks the intraband value at coincident energy.

## bcc-Fe material evidence

All runs use `native_crosscheck=.false.` and the accepted scalar-relativistic
first-order Hamiltonian.  The extra reproducible decks are:

```text
example/exchange_q/bccFe/input_commensurate_12.nml
example/exchange_q/bccFe/input_commensurate_24.nml
example/exchange_q/bccFe/input_legacy_24.nml
```

The finite-temperature metallic result is:

| mesh | q | TT (Ry) | contact (Ry) | total (Ry) | total/q² (Ry Å²) |
|---|---:|---:|---:|---:|---:|
| 12³ | 1/12 | -0.1973170169 | +0.1989376413 | 1.6206243e-3 | 2.4196444e-2 |
| 12³ | 2/12 | -0.1986364080 | +0.1997098660 | 1.0734580e-3 | 4.0067685e-3 |
| 24³ | 1/24 | -0.2039643013 | +0.2043634923 | 3.9919102e-4 | 2.3840203e-2 |
| 24³ | 2/24 | -0.2031085551 | +0.2045772650 | 1.4687099e-3 | 2.1928313e-2 |

The 24³ result is smooth and positive at the first two commensurate points,
but the 12³ to 24³ change shows that this remains an exchange-curvature
diagnostic, not a claimed converged magnon stiffness.

On the same 24³ state and q path, the legacy occupied-only values are:

| q | legacy TT (Ry) | legacy contact (Ry) | legacy total (Ry) | metallic-minus-legacy total |
|---:|---:|---:|---:|---:|
| 1/24 | -0.2038059543 | +0.2038789845 | 4.0674491e-4 | -7.5539e-6 |
| 2/24 | -0.2025504290 | +0.2040930981 | 1.5426690e-3 | -7.3959e-5 |

The legacy values reproduce the starting-point baseline's `4.0674491e-4` and
`1.5426690e-3` scale, while the finite-temperature response changes the
individual terms and residual.  This confirms that the old expression was a
useful zero-temperature regression oracle but not the accepted finite-T
metallic formulation.

## Performance and scaling

For `N_k` k points, `N_q` q points, `N_b` bands, and `N_s` sites, the spectral
contraction scales as `O(N_q N_k N_s² N_b²)` and the explicit torque/contact
assembly is linear in `N_q N_k` times the LMTO bond work (with site-pair
contact factors).  The old endpoint path added `O(N_q N_k)` dense eigensolves;
an exact commensurate path now reuses accepted endpoint eigenpairs and removes
that term.  Storage for the current local implementation is
`O(N_k N_s² N_b²)` for vertex/contact stacks.

Measured material timings on this workspace were:

| mesh | Nq | wall time | endpoint | assembly | contraction |
|---|---:|---:|---:|---:|---:|
| 12³ | 4 | 21.43 s | 0.056 s | 19.442 s | 0.196 s |
| 24³ | 3 | 136.92 s | 0.490 s | 109.686 s | 1.181 s |

The reported wall time includes the accepted-state SCF/DOS setup; component
times are the `exchange_q` path itself and include the one-time exact endpoint
identity audit.  The steady-state commensurate endpoint reuse is therefore
slightly cheaper than the reported endpoint column.  Relative to the clean
starting commit's 24³ legacy deck (89.26 s full-deck wall time, 67.60 s in the
un-instrumented exchange interval), the final metallic deck is not a like-for-
like wall-time comparison: it evaluates Gamma instead of hard-coding it,
uses all finite-temperature bands, and includes the audit.  The measured
component improvement is nevertheless clear: caching the eigenbasis vertices
reduced metallic contraction from 10.613 s to 1.181 s on 24³.  Finite-q
torque/contact assembly remains the dominant cost.  A full 48³ production run
was not required to establish this scaling and was not used as a closure gate.

## Final verdict

**PASS — METALLIC FORMULATION CERTIFIED, PERFORMANCE OPEN.**

The physics is certified for the declared fixed-`mu` finite-H observable and
the gapped reduction, degeneracies, independent metallic fixture, Goldstone
limit, contact cancellation, and commensurate endpoint reuse all close.  The
remaining performance work is local torque/contact assembly optimization; the
native LKAG q comparison remains intentionally blocked and DRESP-04 was not
started.
