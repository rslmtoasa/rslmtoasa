# RS-LMTO Sternheimer and Liouville-Lanczos TDDFT blueprint

**Task:** TDDFT-15

**Status:** design-only; no Sternheimer or Liouville-Lanczos solver is
implemented by this document.

## 1. Decision summary

The future fourth response engine should be a self-consistent
Sternheimer/Liouville-Lanczos engine for the orthogonal, first-order RS-LMTO
Hamiltonian. It should share the response operators, XC/kernel provenance,
metallic occupations, output conventions, and validation fixtures with the
existing Dyson/GF subsystem, while retaining a distinct numerical path.

There are two deliberately separate operating modes:

1. **Bare Sternheimer mode** solves the independent-particle response with the
   induced Hxc field disabled. It returns the same `chi_KS` response object as
   the existing backends and is the backend-equivalence oracle.
2. **Self-consistent Sternheimer/LL mode** includes the induced Hxc field in
   the Hamiltonian action and returns the physical response to an external
   source. It must not be passed through the existing Dyson kernel a second
   time, because that would double count the interaction.

The current `tddft_chi0_backend` interface is intentionally a bare-response
interface. The future physical-response engine should therefore be a sibling
`tddft_response_engine`/`tddft_response_backend` contract, or an explicitly
tagged extension with `contains_kernel=.true.`. It should not pretend that a
self-consistent Sternheimer result is a `chi0` result merely to fit the
existing factory.

The first implementation boundary remains narrow:

- bulk `calctype='B'`;
- the current collinear, SOC-free response boundary (`nsp=1` and
  `magnetic_representation='periodic_nc'` in the TDDFT driver);
- first-order orthogonal `reciprocal_mode='ham_only'` with no metric response;
- no HOH/second-order, GBT, CCOR, orbital-polarization, constraint, SOC,
  noncollinear, or Hubbard response;
- site-projected charge/spin response, with the transverse circular channel
  as the first physical target; and
- the same resolved Fermi level, electronic temperature, q mesh, and eta
  convention as the current reciprocal response route.

Unsupported combinations must fail at the response boundary. In particular,
the presence of a Hamiltonian branch that can be diagonalized is not evidence
that its first-order response operator or metric derivative is available.

## 2. Current RS-LMTO representation to be preserved

### 2.1 Orthogonal electronic problem

For this design, one k-point has an Euclidean-orthonormal coefficient vector

```text
H(k) c(n,k) = epsilon(n,k) c(n,k),
                         c(m,k)^H c(n,k) = delta(m,n).
```

The production first-order reciprocal assembler forms the same matrix used by
the existing band and transition paths,

```text
H(k) = sum_R ee(R) exp(+i k.R),
```

where the phase convention is the one owned by
`reciprocal_fourier.f90`; the Sternheimer code must call that assembler rather
than reimplementing reciprocal phases. The shifted endpoint is the exact
`k+q` point, folded only by reciprocal lattice vectors. It must not be
replaced by the nearest point of the ground-state mesh.

In the spin-major block layout, the magnetic decomposition used by
`hamiltonian_build.f90` is

```text
H(up,up)     = H0 + Hz
H(down,down) = H0 - Hz
H(up,down)   = Hx - i Hy
H(down,up)   = Hx + i Hy.
```

This is the notation used below. The design is not based on extracting a
response vertex from `hamiltonian%hxc`: that array is an assembled magnetic
Hamiltonian quantity, not the functional derivative of the response field.

The current implementation points to reuse are:

| Responsibility | Existing source | Sternheimer use |
| --- | --- | --- |
| response operators and circular factors | [`response_basis.f90`](../source/response_basis.f90), [`response_vertices.f90`](../source/response_vertices.f90) | measurement operators, external source derivatives, and ordered contractions |
| frequency/sign convention | [`tddft_conventions.f90`](../source/tddft_conventions.f90) | resolvent sign and positive eta |
| exact k and shifted-k endpoint assembly | [`reciprocal_fourier.f90`](../source/reciprocal_fourier.f90) | `H(k)`, `H(k+q)`, and q provenance |
| response occupations | [`tddft_chi0.f90`](../source/tddft_chi0.f90), [`reciprocal_occupations.f90`](../source/reciprocal_occupations.f90) | finite-temperature weights and Fermi-surface limits |
| kernel data | [`xc_response_kernel.f90`](../source/xc_response_kernel.f90), [`self_xc_response.f90`](../source/self_xc_response.f90) | induced Hxc field action, keeping XC provenance separate from source/measurement |
| common bare-response result | [`tddft_chi0.f90`](../source/tddft_chi0.f90), [`tddft_backend.f90`](../source/tddft_backend.f90) | bare-mode output and backend equivalence |
| Dyson reference and loss matrix | [`tddft_dyson.f90`](../source/tddft_dyson.f90) | independent physical-response oracle and common spectral products |
| real-space H/GF actions | [`recursion_transport.f90`](../source/recursion_transport.f90), [`green.f90`](../source/green.f90), [`green_lanczos.f90`](../source/green_lanczos.f90) | later preconditioner/action reuse only after residual tests |

### 2.2 Response basis and kernel separation

Let (A) label the response coordinate (site plus charge, Cartesian spin,
or circular component). Define:

- (Gamma_A): the measurement operator used to project an induced density;
- (Lambda_A = partial H/partial b_A): the Hamiltonian generator for a
  field or Hxc-field coordinate; and
- (K_{AB} = partial b_A^{\rm Hxc}/partial d_B): the interaction kernel in
  the response basis.

They are different objects even when their local matrix shapes look similar.
The total first-order Hamiltonian is

```text
delta H_tot = delta H_ext + sum_A Lambda_A delta b_Hxc,A,
delta b_Hxc,A = sum_B K_AB delta d_B.
```

The project’s circular convention is retained exactly:

```text
O+ = sigma_x + i sigma_y = 2 sigma_+
O- = sigma_x - i sigma_y = 2 sigma_-
dH/db+ = O-/2
dH/db- = O+/2.
```

Thus a `PLUS` measurement and `MINUS` source column use the ordered factors
already pinned by [`TDDFT_CONVENTIONS.md`](TDDFT_CONVENTIONS.md). The
Sternheimer engine must use `external_source_operator` for the source and
`response_operator`/the site projector for measurement. It must not put
`dH_xc/dm` on both wavefunction ends and then apply `K` again.

## 3. Response equations in RS-LMTO notation

### 3.1 Authoritative density-matrix equation

For a response with spatial phase (q) and time dependence
(\exp(-i\omega t)), let (delta\rho_{k,q}(\omega)) map the k endpoint to
the k+q endpoint. Define the independent-particle Liouvillian action

```text
L0,q X(k) = H(k+q) X(k) - X(k) H(k).
```

The sign authority for every implementation is

```text
(omega + i eta - L0,q) delta_rho_q(omega)
    = [delta H_tot,q(omega), rho0].
```

In the ground-state eigenbasis, the matrix element is

```text
delta_rho_mn(k,q; omega)
 = (f_n(k) - f_m(k+q)) delta_H_mn(k,q; omega)
   / (omega + epsilon_n(k) - epsilon_m(k+q) + i eta).
```

This is exactly the retarded denominator used by the existing eigenpair
backend. It also makes clear that a metal needs occupation differences, not a
binary occupied/unoccupied assumption.

The response coordinate is projected as

```text
delta d_A(q,omega) = sum_k w_k P_A[delta_rho_k,q(omega)],
```

where `P_A` is the ordered site/orbital/spin measurement contraction. For a
source column (B), use

```text
delta H_ext,q = Lambda_ext,B(q) delta b_ext,B,
delta H_tot,q = delta H_ext,q + sum_A Lambda_A(q) sum_C K_AC delta d_C.
```

The same equations define both the bare mode (`K=0`) and the self-consistent
mode (`K` active). The self-consistent solution is the fixed point of the
Hamiltonian perturbation, not an empirical rescaling of its output.

### 3.2 Resonant first-order equation

For a ground-state state (c_{n,k}), define the resonant response ket
(x^+_{n,k}(q,\omega)) at k+q by

```text
[ H(k+q) - epsilon_n(k) - omega - i eta ] x+_n,k
    = - Q+_n,k(q) delta_H_tot,q(omega) c_n,k.
```

Equivalently, the left-hand operator is
(H(k+q)-[\epsilon_{n,k}+\omega+i\eta]). If `Q+` is the identity, projection
onto a state (c_{m,k+q}) gives

```text
<c_m,k+q | x+_n,k>
 = <c_m,k+q | delta_H_tot,q | c_n,k>
   / (omega + epsilon_n,k - epsilon_m,k+q + i eta).
```

That equation is the direct Sternheimer realization of the current response
denominator. For an insulating static response, `Q+` may be the conduction
projector after a covariant occupied-space gauge has been fixed. For a metal,
the hard projector must not remove partially occupied Fermi-surface states;
the occupation-aware rules in Section 4 apply.

### 3.3 Antiresonant first-order equation

The antiresonant companion is a separate equation at ((-q,-\omega)):

```text
[ H(k-q) - epsilon_n(k) + omega - i eta ] x-_n,k
    = - Q-_n,k(-q) delta_H_tot,-q(-omega) c_n,k.
```

Its denominator is

```text
-omega + epsilon_n(k) - epsilon_m(k-q) + i eta.
```

The response density at ((q,\omega)) contains both the resonant ket and
the antiresonant bra contribution. In compact wavefunction notation:

```text
delta_rho_q(omega) = sum_nk f_nk [
       |x+_n,k(q,omega)> <c_n,k|
     + |c_n,k> <x-_n,k(-q,-omega)|
   ] + delta_rho_FS.
```

The second term is an independently solved bra contribution. A code may solve
the ((-q,-\omega)) ket and map it to this bra, or solve the equivalent row
system directly. It must record which choice was made. In particular, an
unhalved circular operator must not be replaced by an assumed Hermitian
conjugate: the ordered `PLUS/MINUS` source and measurement vertices remain
authoritative.

The density-matrix equation in Section 3.1 is the implementation oracle for
this two-equation bookkeeping. Re-indexing the two terms must reproduce the
occupation-difference numerator and the retarded/advanced identities tested
by [`test_response_conventions.f90`](../tests/unit/test_response_conventions.f90).

### 3.4 Self-consistency loop

For each q and source column, the first production Sternheimer algorithm is:

```text
1. Set delta H_tot = delta H_ext.
2. Solve resonant and antiresonant equations for all active k/state blocks.
3. Form delta rho, then project delta d_A with the response projectors.
4. Apply the supplied Kxc/Hartree kernel to delta d.
5. Build delta H_Hxc from Lambda_A and the induced field coordinates.
6. Mix the new total perturbation and repeat until both field and density
   residuals converge.
7. Project the converged density onto all requested measurement coordinates.
```

The iteration should expose separate residuals for:

- the resonant and antiresonant linear solves;
- the induced density/field fixed point;
- charge conservation when a scalar source is present; and
- the response-basis Ward residual at static Gamma.

The mixing parameter is a numerical convergence control, not a physical
kernel parameter. A failed fixed point must be reported as failed; it must not
be hidden by returning the last unconverged spectrum.

### 3.5 Static limit and Ward identity

The static limit is not a finite-eta dynamic sample. For q=0, use the stable
divided-difference density-matrix limit

```text
delta_rho_mn(0) = [(f_m - f_n)/(epsilon_m - epsilon_n)] delta_H_mn,
```

with the equal-energy limit (f'(\epsilon_n)). This must agree with the
existing static eigenpair implementation before a Sternheimer static path is
trusted.

For the SOC-free transverse ground state, the self-consistent static response
must retain the Goldstone/Ward relation in the selected response basis. The
engine should report the raw residual of

```text
D(0,0) m_GS = [I - chi_KS(0,0) Kxc] m_GS = 0
```

and compare it with the common Dyson path. No empirical lambda, energy shift,
or silently projected Goldstone vector belongs in the Sternheimer solver.

## 4. Metallic Fe/Ni terms

Fe and Ni are not optional edge cases. The implementation must make the
following terms explicit in the design and in its provenance record.

### 4.1 Fractional occupations and exact Fermi level

Use the same Fermi-Dirac function, temperature in kelvin, energy unit in
Rydberg, and resolved response Fermi level as the current reciprocal path.
The response Fermi level is resolved on the complete response mesh after the
ground-state eigenpairs are available. It must not silently reuse a chemical
potential from a different mesh.

For finite q, the factor

```text
f_n(k) - f_m(k+q)
```

must be retained for every transition represented by the Sternheimer density.
It supplies the intraband/small-q contribution and the low-energy Stoner
continuum. Replacing it with a binary occupied/unoccupied filter removes the
metallic response that is important for damping.

### 4.2 Fermi-surface term at q=0

At q=0 in the static or adiabatic limit, the diagonal density response has a
term that is invisible to a hard occupied-space projector:

```text
delta_rho_nn = f'(epsilon_n) [ delta_H_nn - delta_mu ].
```

For a fixed-electron-number calculation, determine the scalar chemical
potential variation from

```text
delta_mu = sum_nk w_k f'(epsilon_nk) delta_H_nn(k)
           / sum_nk w_k f'(epsilon_nk),
```

for the scalar/charge component when the denominator is numerically
nonzero. In a grand-canonical response, `delta_mu=0` is a deliberate ensemble
choice and must be recorded. For a pure transverse perturbation the charge
constraint may vanish by symmetry, but the implementation should still test
it rather than assume it.

At zero electronic temperature, (f') is a Fermi-surface distribution. The
finite-k calculation therefore needs a controlled smearing/temperature and a
k-mesh convergence study. A dynamic q=0 calculation at nonzero omega uses
the commutator equation; the static occupation derivative must not be
injected into it unless the selected ensemble and adiabatic limit require it.

### 4.3 Degeneracies and partially occupied subspaces

- For equal or nearly equal eigenvalues, use the divided-difference limit or
  diagonalize the perturbation within the degenerate block. Never divide by a
  raw energy difference below the configured threshold.
- Do not use a zero-temperature occupied projector for partially occupied
  states. A finite-temperature covariant density-matrix treatment or an
  explicitly occupation-weighted active-space formulation is required.
- If an active-state window is used for a preconditioned implementation, it
  must contain all occupied/partially occupied states and all response states
  needed through the requested (omega_{max}). Increase the window and
  prove convergence; do not interpret the omitted high-energy tail as a
  physical cutoff.
- At small nonzero q, retain the exact shifted k endpoint. The difference of
  occupations across the Fermi surface is the desired intraband physics, not
  a numerical nuisance.

### 4.4 Metal-specific acceptance evidence

The future test surface must include a small spin-split metallic fixture with
states crossing the Fermi level and demonstrate all of the following:

1. the q=0 static result changes by the analytic (f') term when the
   Fermi-surface contribution is disabled;
2. fixed-N and fixed-μ results differ only by the declared chemical-potential
   constraint;
3. the q-to-zero sequence converges to the static limit under k mesh and
   smearing refinement; and
4. the finite-q low-energy loss spectrum retains the intraband/Stoner
   continuum and agrees with the explicit transition oracle.

## 5. Hxc mapping and response provenance

### 5.1 Required induced-field provider

The current [`xc_response_kernel.f90`](../source/xc_response_kernel.f90)
stores a site-scalar transverse ALSDA value and validated derivative slots,
but it does not yet expose a complete electronic-space perturbation action.
The future engine needs a small adapter with the following conceptual
operations:

```fortran
initialize(ground_state, response_basis, kernel_provider)
apply_hxc(response_coordinates, delta_d, delta_hxc)
describe(provenance)
```

`apply_hxc` must construct the first-order Hamiltonian perturbation in the
same coefficient representation as `H(k)`. It may initially support only the
site-transverse ALSDA kernel. It must preserve:

- XC functional and SCF ordering;
- radial versus response-projector populations;
- sign of `B_xc` and the circular one-half factor;
- local site/orbital basis and q phase;
- whether Hartree is active; and
- whether any field is external, XC, constraining, or Hubbard-derived.

For the initial transverse SOC-free path, `Kxc` is static and local in the
selected response basis. If a richer site-orbital kernel is later needed for
the Ward residual, the response coordinate and `Lambda_A` operators must be
refined together. A site-scalar kernel must not be declared successful merely
because the Sternheimer solve converged.

### 5.2 Direct physical response versus common Dyson

The direct Sternheimer fixed point is algebraically equivalent to

```text
chi = chi_KS + chi_KS Kxc chi.
```

This gives two independent comparison modes:

- bare Sternheimer versus `tddft_chi0` eigenpair/K-GF/native-RS results;
- self-consistent Sternheimer/LL versus `tddft_dyson` applied once to the
  same bare response and kernel.

The result metadata must identify `bare` or `self_consistent`, kernel label,
source/measurement basis, q/omega, eta, temperature, Fermi level, solver
residuals, iteration counts, and whether the output already contains Hxc.

## 6. Proposed implementation architecture

### 6.1 Future types and boundaries

Names are suggestions to fit the existing class-like Fortran style:

```fortran
type :: tddft_sternheimer_options
   real(rp) :: eta
   real(rp) :: solver_tolerance
   real(rp) :: self_consistency_tolerance
   integer :: max_solver_iterations
   integer :: max_self_consistency_iterations
   logical :: fixed_electron_number
   logical :: include_hxc
   character(len=24) :: linear_solver
end type

type :: tddft_sternheimer_engine
   ! owns/borrows validated ground-state and response-basis services
contains
   initialize(...)
   evaluate_bare_column(q, omega, source, result)
   evaluate_self_consistent_column(q, omega, source, result)
   evaluate_frequency_batch(q, omega(:), source, result)
   apply_liouvillian(q, vector, output)
   describe(provenance)
end type
```

The production object should own reusable matrix/RHS scratch and immutable
ground-state fingerprints, but it must not own stale eigenvectors or a
mutable global k mesh. Initialization should record the Hamiltonian
`operator_generation`, `reciprocal_mode`, Hamiltonian order, Fermi-level
source, and response-kernel generation. A changed ground-state operator must
invalidate the engine.

The LL engine should consume the same `apply_liouvillian` action and
source/measurement adapters:

```fortran
type :: tddft_liouville_lanczos_engine
contains
   initialize(sternheimer_state, ll_options)
   build_chain(q, source_block, chain)
   evaluate_spectrum(chain, omega(:), response)
end type
```

The interface must support source blocks. One chain per scalar source is
acceptable for the prototype; block chains are the intended route for both
circular channels or a multi-site response.

### 6.2 Ground-state and k-point data flow

```text
validated ground state
  -> H(k), H(k+q), c/occupation data or H actions
  -> response source Lambda_ext and measurement Gamma
  -> bare or self-consistent first-order solve
  -> delta rho / delta d / induced Hxc
  -> common response result and provenance
  -> common loss/mode/sum-rule analysis
```

For the first prototype, retaining the full small-fixture eigenvectors is
acceptable for an oracle and for dense `zgetrf/zgetrs` solves. The production
design should be able to replace the H application and preconditioner without
changing the response equations or result semantics.

### 6.3 Linear solver and preconditioner reuse

The current repository provides:

- LAPACK Hermitian eigensolvers in the reciprocal execution path;
- `zgesv` in the response-space Dyson solve;
- LU factorizations in the existing Green/Lanczos implementation; and
- real-space Hamiltonian vector actions in the recursion module.

It does not currently provide a general high-level GMRES/FGMRES response
solver. The implementation sequence should therefore be:

1. **Prototype:** use a factorized dense shifted system on a small
   orthogonal fixture, preferably `zgetrf` once and `zgetrs` for a block of
   right-hand sides when shifts can be shared. Since each state has a distinct
   (epsilon_{n,k}) shift, do not claim one factorization serves all states
   unless the implementation actually groups equal shifts.
2. **First scalable path:** add a typed complex shifted-system solver with
   block RHS support and explicit residual checks. A restarted GMRES/FGMRES
   path is appropriate for matrix-free H actions; its convergence and
   restart policy must be recorded.
3. **Preconditioner:** begin with site/orbital block-diagonal or local
   Hamiltonian blocks plus the complex shift. Reuse `reciprocal` assembler
   caches and BLAS-3 contractions.
4. **Later RS reuse:** use `recursion%ham_vec_matmul` or a Green/recursion
   resolvent as an approximate preconditioner only after comparing the
   preconditioned residual with a dense solve. The existing native Green
   storage is a spectral/real-space provider, not automatically an inverse for
   the shifted k-space Sternheimer matrix.

No legacy atomic LMTO routine needs to change for this design. If a future
implementation discovers a basis-response term, it must first show that the
fixed orthogonal `ham_only` representation is insufficient and add a separate
typed derivative path with its own finite-difference test.

## 7. Liouville-Lanczos formulation

### 7.1 Liouvillian action

Let (X_q) be a k-resolved density-matrix perturbation. The self-consistent
linearized action is

```text
L_q[X](k)
 = H(k+q) X(k) - X(k) H(k)
   + [ delta H_Hxc[X](k,q), rho0(k) ].
```

The source is (S_B=[\Lambda_{ext,B},\rho_0]), and the response is

```text
delta_rho_B(q,omega)
 = (omega + i eta - L_q)^(-1) S_B(q).
```

The measured susceptibility is the contraction of this density with
(Gamma_A). This formulation includes the resonant and antiresonant sectors
in one operator equation and is the preferred sign/occupation oracle for the
LL implementation.

An equivalent doubled wavefunction representation may be used for efficiency,
but it must retain both sectors and the induced-field coupling between them.
Dropping the antiresonant block is not a controlled low-frequency
approximation for this design.

### 7.2 Non-Hermitian algorithm policy

The Liouvillian is generally non-Hermitian in the ordinary Frobenius inner
product, especially after the occupation commutator and self-consistent Hxc
action are included. The initial LL implementation must therefore assume a
two-sided bi-Lanczos recurrence:

```text
L |v_j>  -> |v_(j+1)>,
<w_j| L -> <w_(j+1)|,
<w_i|v_j> = delta_ij,
```

with explicit monitoring of biorthogonality, look-ahead/breakdown events, and
loss of numerical stability. It must not call the Hermitian `zheev` path just
because the ground-state H is Hermitian.

A pseudo-Hermitian or symplectic reduction may be added later for the
collinear SOC-free case if a response metric proving the reduction is derived
and tested. It is an optimization, not the initial correctness assumption.

The continued fraction/post-processing step has the schematic form

```text
chi_AB(q,omega) = <w_A | (omega + i eta - T_LL)^(-1) | v_B>,
```

where (T_{LL}) is the non-Hermitian tridiagonal or block-tridiagonal
representation and the source/measurement starts obey the same ordered
operator convention as the direct response. Frequency enters only in this
small post-processing solve when the kernel is adiabatic and the chain is
converged.

### 7.3 LL limitations that must remain visible

- A frequency-dependent XC or Hubbard kernel cannot be represented by one
  frequency-independent LL chain without an additional approximation. Such a
  kernel requires frequency-by-frequency Sternheimer/Dyson solves or a
  separately derived rational representation.
- Truncation and loss of biorthogonality can create noncausal poles or violate
  spectral symmetries. Reject/report unstable chains; do not smooth them into
  an apparently physical spectrum.
- The chain is q- and source-dependent. Reusing a chain across q points or
  source operators requires a proven symmetry and matching response vertex.

## 8. Cost model and method choice

Let:

- (N_k): complete response k points;
- (N_q): q points;
- (N_\omega): frequencies;
- (N_s): source columns (one circular column, two chiral columns, or a
  larger site/component block);
- (N_a): active ground-state states or equivalent response blocks;
- (N_i): iterations of one shifted linear solve;
- (N_{sc}): self-consistency iterations per ((q,\omega,source));
- (N_L): LL iterations per ((q,source)); and
- (C_H,C_P,C_K): costs of one H action, projection, and kernel action.

Ignoring small response-space post-processing, a frequency-by-frequency
Sternheimer estimate is

```text
T_S ~ N_q N_omega N_s N_k N_sc
      [ N_a N_i (C_H + C_P) + C_K ].
```

The square-bracket factors are implementation-dependent: a block solver may
share H actions across states and source columns, while a nested self-
consistency loop may apply the kernel once per iteration. The expression is a
planning model, not a benchmark claim.

For an adiabatic LL chain, the corresponding estimate is

```text
T_LL ~ N_q N_s N_k N_L [ N_a (C_H + C_P) + C_K ]
       + N_q N_s N_L N_omega C_contfrac,
```

where `C_contfrac` is a small dense/block continued-fraction operation. The
second term is normally much cheaper than a new electronic solve. LL pays a
larger one-time chain cost and may need two-sided vectors, reorthogonalization,
or block-width work.

Memory has a different tradeoff:

| Method | Electronic work | Frequency scaling | Main memory pressure |
| --- | --- | --- | --- |
| Sternheimer | shifted solves and Hxc fixed points | approximately linear in `N_omega` | active response/RHS blocks and solver workspace |
| LL | repeated Liouvillian actions and one chain | electronic work nearly independent of `N_omega` | two-sided Krylov vectors; optional reorthogonalization history |
| existing Dyson/GF | bare response plus response-space solves | bare-GF dependent; Dyson scales with frequency | GF/response batches and q reuse |

Choose **Sternheimer** when only a few frequencies are needed, the kernel may
be frequency-dependent, solver conditioning needs frequency-specific control,
or a transparent direct comparison is the priority. Choose **LL** when a broad
frequency spectrum is needed for many q/source points, the kernel is adiabatic,
and a stable Liouvillian action is already available. The decision must be
made from measured `N_i`, `N_sc`, `N_L`, memory, and residual data on Fe/Ni-
shaped fixtures, not from the asymptotic formula alone.

The current performance evidence in
[`TDDFT-14_PERFORMANCE_MPI_REPORT.md`](dev/TDDFT-14_PERFORMANCE_MPI_REPORT.md)
should be extended with these counters. Existing q-level MPI ownership is a
natural first distribution for Sternheimer; LL should preserve q/source
parallelism while keeping a chain's electronic vectors local to its owner.

## 9. Shared validation plan

The future engine reuses the existing tests and adds only solver-specific
checks. The minimum matrix is:

| Stage | Evidence | Required comparison |
| --- | --- | --- |
| algebra | `UnitResponseConventions`, `UnitResponseVertices` | Pauli/circular/source factors and ordered vertices |
| analytic bare response | two-level fixture in `UnitTddftChiKS` | Sternheimer pole, sign, advanced relation, spectral weight |
| Hamiltonian tangent | finite-difference rotation fixture | analytic transverse perturbation reproduces the rotated `H` |
| static limit | `UnitTddftChiKS`, `UnitTddftWardConventions` | divided difference and (f') equal-energy limit; no eta substitution |
| bare backend equivalence | `UnitTddftBackendEquivalence`, `UnitTddftGreenChiKS`, `UnitTddftRealspaceGF` | Sternheimer bare response versus eigenpair/K-GF/native-RS on controlled points |
| self-consistency | `UnitTddftDysonModes`, `UnitTddftWard`, `UnitTddftGoldstone` | direct Sternheimer versus one common Dyson enhancement |
| metal behavior | new finite-T crossing fixture plus Fe/Ni smoke tests | Fermi-surface term, fixed-N policy, small-q continuum and convergence |
| LL correctness | new LL-vs-Sternheimer test | converged chain reproduces complex response over selected frequencies |
| causality | common response/loss checks | retarded/advanced symmetry, positive selected-channel loss, no unstable poles |
| production gate | Fe/Ni validation documents | raw Ward convergence, q² stiffness, damping, eta and k/q convergence |

The bare path must be tested with Hxc disabled and with the same source and
measurement operators. The interacting comparison must use the same kernel
record, Fermi level, temperature, q, omega, and eta. Comparing only a scalar
trace is insufficient; compare the full complex response matrix and the loss
matrix eigenmodes where the response dimension permits it.

Every solver result should report:

```text
backend and mode, q/omega, eta, temperature, Fermi level,
Hamiltonian operator generation and branch, source/measurement convention,
active-state/window policy, linear-solver residual and iterations,
self-consistency residual and iterations, LL chain/breakdown diagnostics,
kernel provenance, static/Fermi-surface policy, and convergence status.
```

No material release claim is permitted until the existing three-backend
release gate is already satisfied and the new engine agrees within the
established numerical uncertainty.

## 10. Hubbard-response path — design only

Current TDDFT Hubbard response is unsupported and must remain rejected. This
section records the future interface only; it is not permission to enable
`+U` in the current response driver.

The future Hubbard contribution must be differentiated from the same
ground-state functional used to construct the Hamiltonian. In a localized
subspace with occupation matrix (N^I), the response needs:

```text
delta N^I = projection of delta rho onto the same localized subspace,
delta V_U = derivative of V_U[N^I] applied to delta N^I,
delta H_tot = delta H_ext + delta H_Hxc + delta V_U.
```

The localized projector, spin/orbital ordering, occupation weights, and any
metric must be part of the response state. A scalar Hubbard kernel or an
empirical U rescaling is not a substitute for this derivative.

The antiresonant branch needs special care: under the time-reversal mapping
used in modern noncollinear TDDFPT, the magnetic Hubbard component changes
sign. The implementation should first use the direct ((-q,-\omega)) equation
to establish the sign, then add an optimized time-reversed form only after a
matrix-level equivalence test.

The LL path can include a static/adiabatic Hubbard derivative in its
frequency-independent Liouvillian. A genuinely dynamical (U(\omega)) cannot
be inserted into one fixed LL chain; that case belongs to the
frequency-by-frequency Sternheimer/Dyson path. This distinction is important
because the modern LL Hubbard work demonstrates the value of self-consistent
response, but does not make the current RS-LMTO Hubbard response implemented.

## 11. Minimal prototype and go/no-go criteria

### 11.1 Fixture sequence

The smallest useful prototype has three layers:

1. **Analytic one-orbital spin split:**
   (H=e_0 I+b\sigma_z), one occupied-up and one empty-down state. Test the
   `PLUS`/`MINUS` pole, source one-half, retarded eta, static response, and
   Goldstone normalization against the existing two-level fixtures.
2. **Small periodic metallic toy:** a two-site or few-site first-order
   `ham_only` Hamiltonian with a band crossing the Fermi level, a complete
   2×2×2 or similarly small k mesh, one nonzero q, and finite electronic
   temperature. Compare all matrix elements against the explicit transition
   builder, including the q=0 Fermi-surface branch.
3. **Material smoke test:** bcc Fe, then fcc Ni, with the current supported
   response boundary and a deliberately small frequency/q set. This is a
   diagnostic smoke test, not a production-validation claim.

The initial code milestone should implement bare Sternheimer first, then the
self-consistent fixed-point loop, then an LL chain using the exact same toy
state and Hxc action. Native RS recursion/Green preconditioning is out of
scope until these three layers are green.

### 11.2 Go criteria

Proceed from prototype to production engineering only when all of the
following are demonstrated:

- residuals of the direct linear solves and self-consistency loop meet their
  declared tolerances on the analytic and metallic fixtures;
- bare Sternheimer agrees elementwise with the explicit transition oracle,
  including finite temperature, q shift, circular channel reversal, and the
  static divided-difference limit;
- the metallic q=0 (f') term, fixed-N chemical-potential policy, and small-q
  occupation difference are independently visible and converge;
- self-consistent Sternheimer agrees with the common Dyson result without
  applying either path twice;
- a converged LL chain reproduces the complex self-consistent response at
  selected frequencies and reports stable biorthogonality/breakdown data;
- Ward/Goldstone residuals shrink under systematic numerical refinement and
  no empirical shift/rescaling is needed; and
- unsupported Hamiltonian/response combinations fail explicitly.

### 11.3 No-go criteria

Stop and redesign if any of the following occurs:

- the implementation needs a generalized-overlap metric derivative while
  claiming `ham_only` support;
- a hard occupied projector is the only way to make a metal solve converge;
- direct and antiresonant equations cannot reproduce the density-matrix
  occupation-difference oracle;
- LL loses biorthogonality or produces noncausal poles without a controlled
  recovery policy;
- preconditioning improves runtime only by changing the response, static
  limit, or spectral weight;
- a large Ward correction is required and does not decrease with refinement;
  or
- Hubbard, SOC, noncollinear, HOH, GBT, or CCOR terms are accepted without
  their first-order Hamiltonian derivative and tests.

## 12. Risks and mitigations

| Risk | Consequence | Mitigation/evidence |
| --- | --- | --- |
| non-Hermitian Liouvillian | unstable or noncausal LL spectrum | two-sided bi-Lanczos, biorthogonality checks, explicit breakdown report, pseudo-Hermitian reduction only after proof |
| metallic degeneracy and Fermi surface | wrong q→0 limit or missing damping | finite-T occupation differences, (f') term, degenerate-block treatment, k/smearing convergence |
| shifted-system conditioning | excessive Sternheimer iterations near Stoner poles | complex-shift block solver, residual-based stopping, q/omega-specific preconditioner diagnostics |
| self-consistency oscillation | false convergence or source-dependent spectra | field and density residuals, bounded mixing/DIIS later, hard failure on max iterations |
| memory for LL vectors | out-of-memory or hidden recomputation cost | block width control, q/source ownership, optional two/three-vector storage, measured reorthogonalization mode |
| stale ground-state data | Ward failure and irreproducible response | operator-generation/kernel-generation fingerprints and immutable response state |
| source/measurement double counting | factors-of-two and wrong Goldstone mode | shared convention module, ordered circular fixtures, separate `Gamma`/`Lambda`/`K` types |
| recursion/GF preconditioner mismatch | apparent speedup with changed physics | compare residuals and spectra against dense solves before enabling |
| static versus dynamic kernel | LL cannot represent requested physics | reject frequency-dependent kernel in LL and route to Sternheimer/Dyson |
| premature Hubbard extension | unsupported claims and broken antiresonant sign | keep `+U` guard; implement projector/occupation derivative as a separate milestone |

## 13. Implementation order and completion boundary

The future implementation should proceed in this order:

1. Add/extend the response-state and induced-Hxc action types without changing
   legacy atomic LMTO code.
2. Implement the analytic bare Sternheimer fixture and compare against the
   existing eigenpair response.
3. Add exact q/k endpoint handling and finite-temperature metallic terms,
   including the q=0 static branch.
4. Add the self-consistent Hxc fixed point and compare with common Dyson.
5. Add the matrix-free H action and a validated block iterative solver.
6. Add two-sided LL for adiabatic kernels and compare with converged
   frequency-by-frequency Sternheimer.
7. Only after correctness and material gates, evaluate recursion/GF
   preconditioning and production MPI/memory regimes.

This task is complete at the design boundary. No production solver, input
keyword, Hubbard path, SOC path, or empirical Goldstone adjustment is added by
TDDFT-15.

## 14. References

1. S. Y. Savrasov, “Linear Response Calculations of Spin Fluctuations,”
   *Physical Review Letters* **81**, 2570–2573 (1998),
   [doi:10.1103/PhysRevLett.81.2570](https://doi.org/10.1103/PhysRevLett.81.2570).
   This is the primary LMTO/muffin-tin-orbital variational Sternheimer
   reference and derives the self-consistent first-order wavefunction and
   density response used in this design.
2. L. Binci, N. Marzari, and I. Timrov, “Magnons from time-dependent
   density-functional perturbation theory and nonempirical Hubbard
   functionals,” *npj Computational Materials* **11**, 100 (2025),
   [doi:10.1038/s41524-025-01570-0](https://doi.org/10.1038/s41524-025-01570-0).
   This is the primary modern reference for resonant/antiresonant TDDFPT,
   metallic occupation handling, Liouville-Lanczos frequency amortization,
   time-reversed antiresonant Hubbard signs, and self-consistent Hubbard
   response.
3. T. Gorni, I. Timrov, and S. Baroni, “Spin dynamics from time-dependent
   density functional perturbation theory,” *European Physical Journal B*
   **91**, 249 (2018),
   [doi:10.1140/epjb/e2018-90465-7](https://doi.org/10.1140/epjb/e2018-90465-7).
   This is the LL/TDDFPT background used for the future chain formulation.
4. T. Gorni, O. Baseggio, P. Delugas, I. Timrov, and S. Baroni,
   “turboMagnon — A code for the simulation of spin-wave spectra using the
   Liouville-Lanczos approach to time-dependent density-functional
   perturbation theory,” *Computer Physics Communications* **280**, 108500
   (2022),
   [doi:10.1016/j.cpc.2022.108500](https://doi.org/10.1016/j.cpc.2022.108500).
   This informs the practical non-Hermitian/pseudo-Hermitian chain and
   continued-fraction considerations.
