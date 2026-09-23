# Native transverse rotation dynamics

## Kernel and conventions

For each magnetic site the response keeps both real transverse rotation
coordinates, `theta_ix` and `theta_iy`, with combined index `A=(i,a)`. The
kernel dimension is `2*nsite`. Vertices and contacts are derivatives of the
accepted second-order Hamiltonian `H=B-QB+E_nu`, assembled by
`assemble_lmto_finite_q_torque` and
`assemble_lmto_finite_q_mixed_derivative`.

With the stored vertex orientation
`T_A(m,n)=<m,k+q|T_A(q)|n,k>`, the retarded electronic term is

```text
Pi_AB(q,w) = 1/Wk sum_knm wk (f_nk-f_m,k+q)
              T_A(m,n) T_B(n,m)
              / (w + e_nk - e_m,k+q + i eta)
K_AB(q,w) = C_AB(q) + Pi_AB(q,w)
```

The finite-frequency accumulator preserves `Pi_xy` and `Pi_yx` independently.
`C_AB` is the occupied expectation of the mixed derivative and has no
frequency dependence. The production service pretransforms both torque
vertices into endpoint eigenbases. A separate literal orbital-space/band-loop
oracle evaluates the formula without calling that accumulator.

For the exact static branch, the Fermi divided difference supplies the
coincident-energy limit. In the live endpoint convention, the force-theorem
Hessian is the same-q Cartesian index-symmetric part

```text
H_AB(q) = 1/2 [K_AB(q,0) + K_BA(q,0)]
```

with no conjugation in this reduction. This identity closes against
`force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch` at Gamma and a
generic non-self-inverse q. It is separate from the finite-frequency
reality relation. With the stored real Cartesian coordinates that relation is

```text
K_AB(q,w) = conjg(K_AB(-q,-w))
```

Conjugating the stored band sum, applying `T_A(q)^dagger=T_A(-q)`, and
relabeling `k -> k+q`, `n <-> m` gives the same-index equation above. A `BA`
form results only if the external coordinate labels are transposed along with
the endpoint ordering; that is not how this module stores `K_AB`. The
regression checks the live relation at a non-self-inverse q that translates
the sampled k mesh exactly; the direct band oracle separately uses a generic
q that is not mesh-commensurate.

## Berry and circular coordinates

The q=0 Berry value is evaluated directly from the accepted eigenvectors as
`-i Tr[rho [Gx,Gy]]`, using global spin-1/2 generators. It independently
matches `M_band/2` for g=2. The one-site unitary is

```text
U = 1/sqrt(2) [[1, -i], [1, +i]]
K_pm = U K_xy U^dagger
G_pm = U K_xy^-1 U^dagger
```

For the fresh 12^3 Fe state, the direct Berry commutator is
`1.140190722385`, `M_band=2.280381444770`, and the direct commutator versus
`M_band/2` residual is zero at printed precision. The low-frequency slopes
are `b_plus=-1.140190726270` and `b_minus=+1.140190726267` (relative
normalization residual `3.41e-9`). Thus the small-q positive-frequency
candidate is in the `+` channel for this 12^3 state; both slopes are reported
so the channel is not assumed. The 10^3 sensitivity state changes the
finite-q curvature sign and selects `-`, evidence that Fe material
convergence is still open.

The 12^3 q=0 static Goldstone residual is `9.44e-16`; the circular
off-diagonal residual is `6.48e-20`. No mode is zeroed or subtracted.

## Pole search and bounded Fe results

For each finite q the driver uses a 61-point coarse scan and a 41-point local
refinement, then compares the positive-frequency `Re K=0` crossing, the
minimum of `|K|`, and the maximum of `-Im G`, with `G=K^-1`. It scans
`eta = 1e-4, 2.5e-5, 6.25e-6 Ry`. A row is marked `UNRESOLVED` if those
frequency diagnostics do not agree within regulator/refinement resolution;
no pole energy is promoted to `E(q)` in that case. `FWHM=-1 meV` means the
scan did not provide a valid half-maximum width.

The first fresh primitive bcc-Fe state used a 12 x 12 x 12 full k mesh,
second-order `ham_only`, auto-found EF, 300 K, SOC off, CCOR off, and zero
constraining field. The accepted-state values were:

```text
EF                         -0.08562518103 Ry
Electron-number residual   -1.5241e-11 electrons
M_band                      2.28038144477
SCF magnetic moment         2.28038156289
Berry commutator            1.140190722385
```

The reduced direct q coordinates and actual Cartesian magnitudes are below.
The last three columns are the final-eta (`6.25e-6 Ry`) diagnostics. These
crossings and spectral extrema are shown as diagnostics only because all
three 12^3 rows are unresolved.

| q_fraction `(xi,0,0)` | `|q|` (1/Angstrom) | finite-H adiabatic (meV) | Turek adiabatic (meV) | ReK crossing (meV) | min `|K|` (meV) | max `-Im G` frequency (meV) | max `-Im G` (1/Ry) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1/24 | 0.1294003 | 18.05018 | 13.47215 | 31.80704 | 31.81344 | 30.08363 | -333.6643 |
| 1/16 | 0.1941005 | 15.59925 | 15.27153 | 23.16341 | 23.13889 | 24.69881 | -427.9907 |
| 1/12 | 0.2588007 | 6.15170 | 11.03610 | 5.94256 | 5.93190 | 6.89756 | -1018.161 |

At these q values, the refined frequency resolutions were respectively
`5.527766e-6`, `4.777182e-6`, and `3.379743e-6 Ry` (0.0752, 0.0650,
and 0.0460 meV). ReK crossing shifts over the eta ladder were 0.00160,
0.00130, and 0.01865 meV. The kernel crossings are therefore eta-stable at
this scan resolution, but the `-Im G` maxima are displaced and negative in
the selected 12^3 channels; no FWHM is reported and none of the three is
accepted as a controlled magnon energy.

At the same q=`1/24`, a separate fresh 10^3 state gave finite-H and native
Turek adiabatic estimates of `-13.8291` and `-33.2688 meV`, and a resolved
`-`-channel crossing near `6.8184 meV`. This large change, including the sign
change in the static estimates and a circular-channel switch, confirms that
the 12^3 values are not material-converged. The 10^3 run is only a bounded
mesh sensitivity check.

No `E(q)` dispersion or stiffness `D` is fitted from these unresolved,
mesh-sensitive results. The finite-H versus native Turek finite-q values are
diagnostic comparisons, not a normalization identity.

## Gate status and scope

The production direct-loop oracle closes at `2.78e-17`. The exact-static
Hessian reduction residuals are `4.55e-16` at q=0 and `3.97e-16` at the
generic finite q; the independent q-even mixed x/y contact finite-difference
residual is at most `1.73e-15`. The independent Berry oracle residual is
`2.39e-15`, and the non-self-inverse q/omega covariance residual is
`1.67e-16` (mesh-translation residual `5.55e-17`). The fresh 12^3 response also
closes its q=0 Berry/Goldstone gates and has a positive-frequency causal
kernel crossing. The bounded Fe campaign does not yield three controlled
magnon energies, so the result is:

```text
Implementation              = PASS-B
Fe material convergence     = OPEN
Goldstone correction        = OFF
Strict-ASA L0 dynamics      = NOT RUN
```

The runtime tables and SCF products were written under `/tmp`; none are
committed.
