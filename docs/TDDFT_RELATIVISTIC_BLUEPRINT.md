# RS-LMTO relativistic and noncollinear TD-DFT blueprint

**Task:** TDDFT-16
**Status:** design-only; SOC and noncollinear TD-DFT remain disabled.
**Scope:** formal response structure, relativistic sum rule, validation gates,
and the implementation boundary required before production enablement.

This document is the implementation contract for a future SOC/noncollinear
response engine.  It is deliberately more restrictive than the algebraic
spinor helpers already present in the repository: being able to contract
arbitrary spinors is not evidence that the corresponding relativistic kernel,
torque terms, static sum rule, or production validation exists.

## 1. Decision summary

Do not enable SOC or general noncollinear TD-DFT in the current production
route.  The current `post_processing='susceptibility'` path is the
collinear, scalar-relativistic, first-order, orthogonal (`ham_only`) path.
The production entry point now rejects the following before it builds the
response eigenpairs or accepts a response kernel:

| Ground-state branch | Current control/source evidence | Response status |
| --- | --- | --- |
| Collinear fully relativistic | `control%nsp=2`, rejected by `control%has_soc()` | Unsupported |
| Noncollinear scalar-relativistic | `control%nsp=3`, rejected by `control%is_noncollinear()` | Unsupported |
| Noncollinear fully relativistic | `control%nsp=4`, rejected by `control%is_noncollinear()` | Unsupported |
| Collinear input carrying nonzero `xi_p`/`xi_d` | `has_soc` potential audit in `calculation.f90` | Unsupported |
| Local-axis response representation | `hamiltonian%local_axis`, rejected at the response boundary | Unsupported |
| `magnetic_representation` other than `periodic_nc` | response-boundary guard | Unsupported |
| GBT, explicit texture, HOH, CCOR, orbital polarization, Hubbard, constraints | response-boundary guards | Unsupported |

The implementation sequence below must not relax these guards until all of
the following are available together:

1. a full spinor `chi_KS` response and its source/measurement convention;
2. a kernel in the same charge-spin basis, including all allowed couplings;
3. a relativistic static sum-rule residual and torque/anisotropy source;
4. independent static anisotropy and dynamic gap validation; and
5. a provenance record that identifies the relativistic response branch.

The existing `tddft_four_component_mod` is a useful algebraic fixture and a
future building block.  Its arbitrary-spinor bare response and local-frame
kernel helpers do not authorize production SOC/noncollinear use.  The
production driver remains the capability boundary.

## 2. Current code boundary and evidence

### 2.1 Meaning of `nsp`

The control module defines the global electronic mode, independently of the
number of local XC density channels:

```text
nsp = 1  collinear scalar-relativistic
nsp = 2  collinear fully relativistic (l.s)
nsp = 3  noncollinear scalar-relativistic
nsp = 4  noncollinear fully relativistic (l.s)
```

This distinction matters because the local atomic/XC path may still carry two
spin density channels when `nsp=1`.  A response guard must use the global mode,
not the local radial channel count.

The relevant source is [`source/control.f90`](../source/control.f90), where
`is_collinear`, `is_noncollinear`, and `has_soc` are executable capability
queries.  The Hamiltonian stores SOC in `lsham`; the onsite construction uses
the orbital angular-momentum matrices and spin blocks in
[`source/hamiltonian_build.f90`](../source/hamiltonian_build.f90).

### 2.2 Current response entry point

The current call sequence is:

```text
post_processing = susceptibility
        |
        +-- parse &tddft
        +-- prepare the normal ground-state stack
        +-- reject nsp != 1, SOC, local-axis and unsupported Hamiltonian branches
        +-- build the first-order ham_only reciprocal Hamiltonian
        +-- solve k and exact k+q spinor eigenpairs
        +-- construct the selected bare chi_KS backend
        +-- construct the current collinear kernel/Xi route
        +-- solve Dyson and write spectra
```

The TDDFT-specific guards are in
[`source/calculation.f90`](../source/calculation.f90), immediately after
generic stack preparation and before the response-mesh eigenpair solve.  They
are intentionally fatal rather than warnings or fallbacks.  In particular,
the presence of a diagonalizable spinor Hamiltonian does not imply that its
response equations are implemented.

The four-component test
[`tests/unit/test_tddft_four_component.f90`](../tests/unit/test_tddft_four_component.f90)
contains synthetic spin-mixed fixtures.  Those tests validate matrix
contraction and local-frame algebra only.  They are not a relativistic
material validation and must not be interpreted as production enablement.
The source-contract checks in
[`tests/unit/test_tddft_dispatch.py`](../tests/unit/test_tddft_dispatch.py)
verify that the SOC/noncollinear guards precede the response eigenpair work.

### 2.3 Existing object capabilities versus production capability

| Object/path | What it can currently do | What it cannot claim |
| --- | --- | --- |
| `response_vertices_mod` | Contract charge and Cartesian/circular Pauli operators with arbitrary complex spinors | A complete relativistic kernel or torque response |
| `tddft_four_component_mod` | Build a synthetic site-major `(n,mx,my,mz)` bare tensor and rotate a local ALSDA fixture | Enabled SOC/noncollinear production TD-DFT |
| `xc_response_kernel_mod` | Store collinear circular and optional full ALSDA derivative data | SOC-induced anisotropy, orbital response, or a relativistic Ward repair |
| `tddft_goldstone_mod` | Diagnose/reconstruct the collinear no-SOC static identity | Enforce a zero mode when SOC or an external field breaks spin rotation |
| reciprocal `G(k,z)` and real-space GF paths | Supply one-electron propagator infrastructure | A validated full spinor GF bubble plus relativistic kernel |

No collinear-only scalar `K_perp` may be passed through Dyson for a SOC or
noncollinear state.  The current `K_perp` is tied to the unhalved circular
convention and the collinear local exchange splitting; its use in a future
relativistic calculation would discard the torque and channel-coupling terms.

## 3. Hamiltonian and response conventions

### 3.1 Electronic basis and Hamiltonian

The normal electronic coefficient vector is spin-major:

```text
(orbital_1 up, ..., orbital_N up, orbital_1 down, ..., orbital_N down).
```

For a scalar-relativistic spin-dependent Hamiltonian, the code’s block
assembly is

```text
H(up,up)     = H0 + Hz
H(down,down) = H0 - Hz
H(up,down)   = Hx - i Hy
H(down,up)   = Hx + i Hy.
```

Equivalently, before SOC is added,

```text
H(k) = H0(k) sigma_0 + Hx(k) sigma_x + Hy(k) sigma_y + Hz(k) sigma_z.
```

The fully relativistic branch adds an orbital-spin matrix `H_SOC(k)` built
from `L . sigma` in the LMTO orbital basis.  It is not generally diagonal in
the global spin labels.  The future engine must treat the complete complex
Hermitian `H(k)` as the operator; it must not split the problem into up/down
bands or assume a conserved `S_z`.

The first relativistic implementation should use the same first-order
Hamiltonian assembler and exact shifted endpoint as the current response
path:

```text
H(k) c_n(k)       = e_n(k) c_n(k)
H(k+q) c_m(k+q)   = e_m(k+q) c_m(k+q)
```

The shifted point is folded by reciprocal lattice vectors only.  It must not
be replaced with the nearest normal-mesh point.  The first implementation
must retain `reciprocal_mode='ham_only'`; a generalized-overlap relativistic
response needs a separate metric derivation and remains rejected.

### 3.2 Density and measurement projectors

Let `P_a` be the site/orbital projector for response site `a`, and define

```text
Gamma_(a,A) = P_a tensor sigma_A,
A in {0,x,y,z}.
```

`A=0` is charge and `A=x,y,z` are the code’s spin-population components.  In
the code’s density convention,

```text
n  = rho_upup + rho_downdown
mx = rho_updown + rho_downup
my = Im(rho_downup - rho_updown)
mz = rho_upup - rho_downdown.
```

The projected observable is

```text
d_(a,A) = Tr[Gamma_(a,A) rho].
```

These are electron-number/spin-population variables.  No `mu_B` factor is
present in the current response operators or XC kernel.  A future physical
magnetization API may add units, but the conversion must be explicit and
must not be inserted into only the SOC path.

### 3.3 External source versus self-consistent kernel

For a magnetic source coordinate `b_(b,B)`, define an electronic-space source
operator

```text
Lambda_(b,B) = dH / db_(b,B).
```

The measurement `Gamma` and source `Lambda` are distinct interfaces.  The
interaction kernel is a third object,

```text
K_(a,A),(b,B) = d b_Hxc(a,A) / d d_(b,B).
```

For the current circular convention,

```text
O+ = sigma_x + i sigma_y = 2 sigma_+
O- = sigma_x - i sigma_y = 2 sigma_-
m+ = mx + i my
m- = mx - i my
dH = (b+ O- + b- O+) / 2.
```

Thus `external_source_operator(PLUS)` is the opposite circular operator with
the explicit one-half source factor.  A bare response matrix may store the
ordered unhalved `O+`/`O-` vertices, but its metadata must state whether the
source half is included.  The code must never put `dH_xc/dm` on both response
vertices and then apply `K` again in Dyson.

## 4. Full spinor Kohn-Sham response

### 4.1 Authoritative Lehmann expression

For a site-major response index `I=(a,A)` and `J=(b,B)`, the first
relativistic bare response is

```text
chi_KS[I,J](q,omega) =
  sum_(k,n,m) w_k
  (f_nk - f_m,k+q)
  <n,k | Gamma_I | m,k+q>
  <m,k+q | Lambda_J | n,k>
  ------------------------------------------------------------
  omega + e_nk - e_m,k+q + i eta.
```

The retarded denominator and finite-temperature occupations are inherited
from the current convention contract.  The two endpoint spinors remain
separate; the right factor is an ordered matrix element and is not generated
by blindly conjugating the left factor.

If the implementation chooses to store both response vertices as unhalved
Pauli operators, then the source-coordinate conversion is applied exactly
once at the source boundary.  A test must compare the stored circular block
to the Cartesian tensor after applying the documented linear transformation.

### 4.2 Response dimension and coupling

For `N_site` response sites, the active full spin-charge tensor has dimension
`4*N_site` with site-major ordering

```text
(site_1, charge), (site_1, mx), (site_1, my), (site_1, mz), ...
```

The response is a block matrix

```text
        charge   mx       my       mz
charge  chi00    chi0x    chi0y    chi0z
mx      chix0    chixx    chixy    chixz
my      chiy0    chiyx    chiyy    chiyz
      mz      chiz0    chizx    chizy    chizz
```

where every entry is itself a site-site matrix.

For a collinear SOC-free state, symmetry can reduce this tensor to the
familiar transverse and charge/longitudinal sectors.  For SOC or a general
noncollinear state, the implementation must retain:

- charge-spin mixing allowed by the Hamiltonian and the XC Hessian;
- all Cartesian spin components, or a mathematically equivalent full complex
  circular representation;
- `xy` and `yx` terms, including ellipticity/chirality;
- site and orbital projectors required to resolve local anisotropy; and
- the full complex response, including non-Hermitian dissipative terms.

Circular channels are a basis change, not a decoupling rule.  In particular,
the code must not assume that `chi^{+-}` alone contains the SOC response.

### 4.3 Kernel and Dyson equation

The physical response uses one kernel in the same response basis:

```text
chi = chi_KS + chi_KS K chi
    = (I - chi_KS K)^(-1) chi_KS.
```

The future relativistic kernel contains:

1. charge-charge Hartree and any consistently projected electrostatic terms;
2. the complete adiabatic XC Hessian in `(n,mx,my,mz)` coordinates;
3. local-frame rotations when the equilibrium moment is not along global z;
4. any external-field or orbital-field source terms in `Lambda`; and
5. no empirical scalar shift or fitted Goldstone correction.

For a spin-rotation-invariant local XC functional, SOC does not turn the XC
Hessian into an arbitrary anisotropic scalar.  Rather, SOC enters the
one-electron Hamiltonian, the spinor eigenvectors, and the non-XC torque
source.  A local-frame ALSDA Hessian may still be useful, but its derivative
provenance must be independent of the SOC capability flag and its use must
be validated in the complete response.

The first implementation should use stable linear solves and report the
condition of `I-chi_KS K`.  Near-singularity caused by a physical resonance
must be distinguished from a failed factorization or an invalid kernel.

## 5. Relativistic static sum rule and torque source

### 5.1 Jülich identity

The relativistic extension in dos Santos Dias et al. starts from a
non-spin-diagonal Kohn-Sham Hamiltonian and separates the exchange-correlation
splitting from all other contributions.  In their half-normalized circular
notation the static relation is

```text
m_z(r) - integral dr' chi_KS,+-(r,r';0) K_perp(r') m_z(r')
    = delta_m_z(r),
```

with `K_perp = 2 B_xc / m_z` in that paper’s convention.  The right-hand
side contains the contribution from SOC and external fields (and any other
terms included in the effective splitting) after the XC part has been
removed.  It is therefore not zero in general.  SOC removes the ordinary
global spin-rotation Goldstone condition and produces an anisotropy gap even
at zero applied field.

The paper also reports weak coupling between longitudinal and transverse
magnetic response and shows that both resonance energy and lifetime are
affected by SOC.  These effects are part of the target response, not
post-processing corrections.

### 5.2 Translation to the RS-LMTO convention

The code measures unhalved `m+/-` with `O+/- = sigma_x +/- i sigma_y` and
stores the circular site kernel for the ordered unhalved response as

```text
K_code_perp = B_xc / (2 M_site).
```

The factor is not interchangeable with the paper’s `K_perp`: the paper uses
half-normalized circular spin operators while the code’s response vertices
are unhalved and the external source contributes the explicit one-half.
For the code’s current ordered circular matrix, the static diagnostic must
therefore be written as

```text
delta_m_code = m_code
             - chi_KS,code(0) * b_xc,code
             - (non-XC torque/source contribution),
```

or, for the scalar site diagnostic when the source convention is fixed,

```text
m_code - chi_KS,code(0) * (K_code_perp M_site) = delta_m_code.
```

This is the only acceptable way to compare the existing site-projected
quantities with the relativistic literature.  A factor-of-two or factor-of-
four conversion must be recorded in metadata and tested, never hidden in a
new kernel constant.

The future full-tensor form should be expressed using infinitesimal spin
rotation generators.  For a rotation axis `u_a`, the rigid spin variation is

```text
delta m_a = theta_a (u_a cross m_a).
```

The exchange-correlation field rotates in the same local frame.  The
variation of `H_SOC`, external Zeeman/orbital terms, and any other non-XC
Hamiltonian contribution produces a residual torque source.  The static
identity must compare the response to that complete source with the directly
computed torque/field variation.  It must not manufacture a zero eigenvalue
by projecting the residual away.

### 5.3 Ward policy

The existing no-SOC Goldstone policies (`diagnose`, controlled sum-rule or
projection paths) are not valid relativistic policies.  For SOC/noncollinear
response:

- no zero-mode enforcement is applied;
- no empirical frequency shift is applied;
- no global scalar lambda is fitted to a target gap;
- any static reconstruction must solve the relativistic identity with the
  non-XC residual explicitly present; and
- the raw residual, torque source, kernel provenance, and any optional
  diagnostic reconstruction are written separately.

A large residual is a failed derivation/convergence result, not evidence that
the gap should be tuned.

## 6. Required ground-state quantities

The production relativistic engine must receive or reconstruct the following
from the same ground-state Hamiltonian and XC calculation:

| Quantity | Required use | Current status |
| --- | --- | --- |
| Complete complex `H(k)` including `lsham` | Spinor eigensystem and KS response | Ground-state capability exists; response path rejects it |
| Exact `H(k+q)` | Endpoint-resolved transitions | Existing arbitrary-k service exists |
| Spinor eigenvectors and occupations | Full `4x4` response contraction | Existing eigenpair storage can supply arbitrary spinors |
| Site/orbital projectors `P_a` | Charge/spin measurement | Existing site projectors; orbital-resolved production is future |
| `m_a=(mx,my,mz)` and signed conventions | Rotation generators and sum rule | Ground-state moments exist; response provenance needs explicit packaging |
| Local `B_xc` or its LMTO-consistent source representation | XC part of static identity | Current scalar/circular provider is collinear-only |
| `dB_Hxc/d(n,m)` | Full Dyson kernel | Full ALSDA fixture exists; relativistic validation absent |
| SOC torque / non-XC field variation | `delta_m` right-hand side | Not implemented |
| External spin and orbital field vertices | Field-dependent gap and g-factor checks | Spin source helper exists; orbital response source is future |
| Metric and `delta S` if overlap is used | Generalized response | Not implemented; `ham_only` remains mandatory |

For noncollinear states, the local magnetization direction is not merely an
argument to rotate a scalar kernel.  It determines the local frame for the
XC Hessian and the rigid-rotation generators, while the spinor Hamiltonian
determines the non-XC torque.

## 7. Validation hierarchy before enablement

### 7.1 Algebraic and one-electron fixtures

1. **Pauli decomposition:** construct random Hermitian `H0,Hx,Hy,Hz` and
   verify the spin-major block assembly and inverse decomposition.
2. **SOC block:** construct an orbital `L.sigma` fixture and verify Hermiticity,
   spin off-diagonal blocks, and the finite-difference change under a spin
   rotation.
3. **Full tensor:** compare the `4*N_site` spinor Lehmann sum with a direct
   density-matrix finite difference for a two-level/two-orbital fixture.
4. **Circular transformation:** compare Cartesian and circular responses with
   the explicit source half and ordered vertices.
5. **Charge coupling:** use a spin-mixed fixture where `chi_(0,x)`,
   `chi_(x,0)`, and longitudinal/transverse off-diagonal terms are nonzero;
   assert they are not silently discarded.
6. **Causality and symmetry:** test retarded/advanced identities and the
   expected complex response symmetries at `q` and `-q`.
7. **Metric guard:** verify that a generalized-overlap input is rejected until
   the metric response is implemented; no `zS-H` result may enter by fallback.

### 7.2 Static anisotropy and torque checks

The SOC gap must be validated against an independent static calculation.  The
minimum fixture is a small magnetic site with a controlled SOC strength and
two orthogonal orientations.

For each orientation:

1. converge the same ground-state XC/Hamiltonian branch;
2. compute the total-energy difference or a documented magnetic-anisotropy
   energy (MAE) with identical numerical settings;
3. compute the static torque as the finite derivative of energy with respect
   to the rotation angle;
4. compute the relativistic static response residual from the full tensor; and
5. compare the dynamic resonance/gap with the static curvature, allowing for
   the documented response convention and dynamical renormalization.

The check must vary SOC strength, magnetization direction, and `eta`.  It must
show that the gap tends to zero as SOC tends to zero, while the zero-SOC
response recovers the ordinary rigid-rotation identity.  This is a limiting
test, not permission to shift a finite-SOC spectrum to zero.

### 7.3 External-field checks

With SOC enabled, apply a small field parallel and perpendicular to the
equilibrium moment and verify:

- the field contribution appears in the static residual and provenance;
- the gap/branch changes with the correct sign and source factor;
- spin-only and spin-plus-orbital field couplings are distinguishable;
- the response changes smoothly through a canted equilibrium direction; and
- the extracted g-factor is a response result, not a fixed `g=2` conversion.

### 7.4 Dissipation and mode analysis

The KS spectrum remains the Stoner/electron-hole baseline.  The enhanced
response must retain the full complex loss matrix and report:

- mode eigenvectors and site/orbital character;
- transverse/longitudinal/charge mixing;
- chirality and ellipticity;
- resonance energy and linewidth where a quasiparticle peak exists; and
- overdamped/non-quasiparticle classification when it does not.

The linewidth must be tested against `eta` and energy-contour convergence so
numerical broadening is not reported as a physical SOC lifetime.

## 8. Software implementation sequence

### Phase R0 — capability boundary and provenance

Keep the current hard guards and expose one response-boundary capability
record containing:

```text
global mode: nsp, collinear/noncollinear, SOC
magnetic representation and local-axis state
Hamiltonian metric/order: ham_only or generalized
Hxc branches: HOH, CCOR, orbital polarization, Hubbard, constraints
response basis and source convention
```

The record should be emitted on rejection and with every future relativistic
response product.  It must distinguish `nsp=2/4` SOC from `nsp=3` scalar
noncollinearity and also catch a collinear input with nonzero `xi_p/xi_d`.

### Phase R1 — full spinor bare response

Add a backend mode that returns the full site-major `4*N_site` `chi_KS` with
ordered measurement/source vertices.  Reuse the current arbitrary-k endpoint
service and finite-temperature occupation code.  Implement the direct
eigenpair reference first; only then adapt K-space and native real-space GF
products.

The backend metadata must identify:

```text
spinor Hamiltonian branch, SOC strength/source, q and -q endpoints,
response basis, source normalization, eta role, metric, occupations,
orbital/site projection, and convergence state.
```

### Phase R2 — full kernel and field vertices

Implement a response-kernel provider that is separate from the measured
`chi_KS` vertices.  The provider must supply the full XC Hessian, local-frame
rotation, Hartree charge block, and external/SOC field source records.  The
existing collinear `K_perp` accessor must be rejected by type/capability when
the relativistic branch is active.

If orbital moment response is included, add explicit orbital measurement/source
coordinates or a separate torque operator.  Do not claim that a spin-only
`(0,x,y,z)` tensor contains orbital Zeeman response automatically.

### Phase R3 — relativistic sum-rule diagnostic

Implement the direct static residual with `eta=0` and the appropriate real
divided-difference/static GF limit.  Compute the non-XC torque/source term
directly from the Hamiltonian decomposition and external-field provenance.
Compare:

```text
full static response applied to the complete infinitesimal rotation source
versus the directly rotated density/moment.
```

Report separate XC, SOC, external-field, and numerical residual components.
No projection or reconstruction is enabled by default.

### Phase R4 — fixture and material release gates

Run the algebraic fixtures, static anisotropy test, and field tests before
material calculations.  Then add at least one SOC-bearing 3d adatom/film
fixture with a known anisotropy trend and one canted noncollinear fixture.
Only after those pass may a material response route be exposed behind an
explicit opt-in capability flag.

### Phase R5 — backend equivalence and production opt-in

Require agreement between the direct eigenpair, K-space GF, and native
real-space GF bare tensors at selected `(q,omega)` points.  Require the same
static residual and gap trends across backends.  Keep a runtime option for
diagnostic raw products and make the first production default fail closed when
the relativistic provenance is incomplete.

## 9. Proposed API boundary

Names are illustrative and should follow the existing class-like Fortran
style.  The important contract is the separation of capabilities and data,
not a new inheritance hierarchy.

```text
type :: tddft_relativistic_capabilities
   logical :: enabled = .false.
   logical :: full_spinor_chi0 = .false.
   logical :: full_charge_spin_kernel = .false.
   logical :: soc_torque = .false.
   logical :: external_spin_source = .false.
   logical :: orbital_source = .false.
   logical :: metric_response = .false.
   character(...) :: reason
end type

type :: tddft_relativistic_metadata
   character(...) :: global_mode
   character(...) :: hamiltonian_branch
   character(...) :: response_basis
   character(...) :: source_normalization
   character(...) :: kernel_provenance
   real :: soc_scale
   real :: torque_residual
   real :: static_sum_rule_residual
   logical :: gap_validation_passed
end type
```

The backend contract should accept batches of q/frequency points and return
the common response result plus the metadata.  It should preserve separate
`k`/`k+q` and `R`/`-R` endpoint provenance.  The physical response engine must
declare whether it already contains the kernel; a self-consistent
Sternheimer-style result must not be sent through Dyson a second time.

## 10. Support matrix after the future implementation

| Feature | Current branch | Earliest acceptable enablement condition |
| --- | --- | --- |
| `nsp=1`, SOC-free, collinear, `ham_only` | Supported transverse baseline | Existing core release gates |
| `nsp=2`, collinear SOC | Rejected | R1–R4 plus SOC gap/torque validation |
| `nsp=3`, noncollinear scalar-relativistic | Rejected | Full spinor tensor, local frames, noncollinear kernel and fixtures |
| `nsp=4`, noncollinear SOC | Rejected | Both preceding capabilities and SOC/noncollinear equivalence |
| `xi_p/xi_d != 0` with `nsp=1` | Rejected | Constructed-Hamiltonian SOC detection and R1–R4 |
| Generalized overlap | Rejected | Metric-consistent vertices, occupations, static identity and tests |
| GBT / explicit texture | Rejected | Endpoint/tangent/source derivation for the selected representation |
| HOH / CCOR / orbital polarization | Rejected | Complete differentiated Hamiltonian response |
| Hubbard U/V | Rejected | Self-consistent Hubbard response kernel and validation |
| External constraints/common moment | Rejected | Differentiated external-field source and sum-rule treatment |
| Orbital magnetic response | Not represented | Explicit orbital basis/source and independent validation |

## 11. Risks and non-goals

- **False Goldstone correction:** SOC gaps must not be removed by the existing
  no-SOC projection or by an empirical shift.
- **Hidden factor errors:** circular operator/source factors differ between the
  code and literature conventions; every conversion needs a fixture.
- **Incomplete torque:** rotating only the XC field while leaving `L.sigma`,
  Zeeman, or constraint terms unchanged produces a false static residual.
- **Charge leakage:** dropping charge components because they vanish in the
  collinear limit loses SOC-induced coupling.
- **Local-frame ambiguity:** a site-local frame must have a deterministic
  orientation and its transformation must be applied to both kernel and
  response projectors.
- **Metric contamination:** an overlap-capable eigensolver is not enough; the
  density metric and first-order metric response must be derived.
- **Damping conflation:** `eta`, contour errors, and physical SOC/electron-hole
  damping must be reported separately.
- **Orbital omission:** a spin-density response is not automatically a
  spin-plus-orbital response.
- **Overbroad scope:** this task does not implement the relativistic solver,
  SOC kernel, orbital response, generalized-overlap response, Hubbard response,
  or material release validation.

## 12. References

Primary relativistic reference:

- M. dos Santos Dias, B. Schweflinghaus, S. Blügel, and S. Lounis, “Relativistic
  dynamical spin excitations of magnetic adatoms,” *Phys. Rev. B* **91**,
  075405 (2015), [arXiv:1501.05509](https://arxiv.org/abs/1501.05509),
  [doi:10.1103/PhysRevB.91.075405](https://doi.org/10.1103/PhysRevB.91.075405).
  The paper derives the non-spin-diagonal static sum rule and identifies the
  SOC/external-field contribution that replaces the zero-SOC Goldstone
  condition.

Related implementation references from the master blueprint:

- S. Lounis, A. T. Costa, R. B. Muniz, and D. L. Mills, *Phys. Rev. Lett.*
  **105**, 187205 (2010), and *Phys. Rev. B* **83**, 035109 (2011), for the
  real-space response and static sum-rule architecture.
- P. Buczek, A. Ernst, and L. M. Sandratskii, *Phys. Rev. B* **84**, 174418
  (2011), for the generalized susceptibility, Dyson layer, loss matrix, and
  Stoner/collective separation.
- L. M. Sandratskii and P. Buczek, *Phys. Rev. B* **85**, 020406(R) (2012),
  and M. M. Odashima et al., *Phys. Rev. B* **87**, 174420 (2013), for AF/ferri
  channel and chirality validation targets.

## TDDFT-16 acceptance evidence

- [x] Current unsupported SOC and noncollinear cases fail explicitly at the
  production response boundary; source-contract tests verify guard ordering.
- [x] Required full spinor response structure is documented in Sections 3–4.
- [x] Relativistic sum-rule implications are mapped to the RS-LMTO circular
  normalization in Section 5.
- [x] SOC gap, anisotropy, torque, and external-field validation plans exist
  in Section 7.
- [x] No collinear-only kernel is silently reused; the current scalar
  `K_perp` remains outside the relativistic capability boundary.
- [x] Feature remains guarded; this document adds no relativistic production
  route.
