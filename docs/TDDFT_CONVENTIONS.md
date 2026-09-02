# RS-LMTO transverse TDDFT conventions

This document is the executable convention contract for the initial linear-
response TDDFT implementation.  It applies to the collinear, SOC-free,
orthogonal (`ham_only`) response path.  Generalized-overlap/HOH, GBT, CCOR,
Hubbard-corrected response, noncollinear response, and SOC response remain
explicitly unsupported until their operators and metrics are derived.

The corresponding compact regression is `tests/unit/test_response_conventions.f90`.
It is registered as `UnitResponseConventions` and is intended to fail when a
circular factor, frequency sign, channel assignment, source vertex, or
spectral-weight normalization changes.

For multi-sublattice transverse runs, `&tddft circular_channel='both'`
(the default) evaluates the ordered `+-` and `-+` responses independently.
The primary `+-` files retain the historical names; reverse-channel products
carry the `_minus_plus_` suffix and every response metadata record identifies
its ordered channel.

## 1. Electronic basis and Pauli algebra

The spin ordering of an electronic block is

```text
(orbital 1 up, ..., orbital N up, orbital 1 down, ..., orbital N down).
```

For one orbital the code uses

```text
sigma_0 = [[1, 0], [0, 1]]
sigma_x = [[0, 1], [1, 0]]
sigma_y = [[0,-i], [i, 0]]
sigma_z = [[1, 0], [0,-1]] .
```

`response_basis_mod` repeats these blocks for the supplied orbital count.
Site response channels in `response_vertices_mod` apply the same local block
to a site projector; each local site block is spin-major even when the full
unit-cell vector is assembled from several sites.

The measured circular variables are defined without a one-half:

```text
O_+ = sigma_x + i sigma_y = 2 sigma_+
O_- = sigma_x - i sigma_y = 2 sigma_-
m_+ = m_x + i m_y = <O_+>
m_- = m_x - i m_y = <O_->
```

The conventional ladder matrices are `sigma_+ = O_+/2` and
`sigma_- = O_-/2`.  In code, `response_operator(PLUS/MINUS)` returns `O_+/-`,
while `ladder_operator(PLUS/MINUS)` returns the conventional half-normalized
ladder matrix.  The factors are centralized as
`tddft_circular_operator_factor = 2` and
`tddft_circular_source_factor = 1/2`.

For a Hermitian 2x2 density matrix, `math_mod%rho2nm` implements

```text
n  = rho_upup + rho_downdown
mx = rho_updown + rho_downup
my = Im(rho_downup - rho_updown)
mz = rho_upup - rho_downdown
```

Consequently `m_+ = 2 rho_downup` and `m_- = 2 rho_updown`.  This is also the
trace convention used by the response transition vertices:
`<a|O|b>`, with the density expectation represented as `Tr(rho O)`.

## 2. Hamiltonian and external source vertex

The spin-dependent Hamiltonian assembled by `hamiltonian_build.f90` is

```text
H_spin = H_0 sigma_0 + H_x sigma_x + H_y sigma_y + H_z sigma_z
```

and therefore has matrix blocks

```text
H(up,up)     = H_0 + H_z
H(down,down) = H_0 - H_z
H(up,down)   = H_x - i H_y
H(down,up)   = H_x + i H_y .
```

The same mapping is used for the magnetic part stored in the LMTO `hxc`
array.  `hxc` is an assembled Hamiltonian block, potentially with LMTO
matrix/hopping structure; it is not, by itself, a local response kernel.

Internal magnetic coefficients have energy units of Rydberg.  For an external
magnetic perturbation written in those energy coordinates,

```text
delta H_ext = delta b_x sigma_x + delta b_y sigma_y + delta b_z sigma_z
            = (delta b_+ O_- + delta b_- O_+)/2
delta b_+ = delta b_x + i delta b_y
delta b_- = delta b_x - i delta b_y .
```

Thus the actual source derivatives are

```text
dH/db_x = sigma_x              dH/db_y = sigma_y
dH/db_+ = O_-/2                dH/db_- = O_+/2 .
```

`response_basis_mod%external_source_operator` returns these derivatives.
This source operator is distinct from the measured transition operator in a
circular channel: the `chi^{+-}` right vertex is the unhalved `O_-`, and the
factor one-half belongs to the circular Hamiltonian/source coupling.  A
finite-difference rotation of `H_0 + b sigma_z` about the y axis gives
`dH/dtheta = b sigma_x`; `UnitResponseConventions` checks this directly.

## 3. Magnetization sign, units, and XC provenance

The internal response population is the electron spin population

```text
M_site = n_up - n_down
```

so a site with more occupied up states has positive `M_site`.  The projected
DOS path in `bands_mod%calculate_magnetic_moments` computes the same Cartesian
spin populations and stores the norm only as a separate display quantity.
TDDFT response vectors use the signed `P_site sigma_z` population; they must
not replace it with `abs(M_site)` or a normalized moment direction.

The numerical identification of one spin electron with one Bohr magneton is
an output convention (`g=2`), not an internal response normalization.  No
`mu_B` factor is present in `rho2nm`, `VXC0SP`, the Hamiltonian coefficients,
the response vertices, or the circular kernel.  `B_xc` below is an energy
coefficient in Rydberg, not a field in Tesla.  Conversion to Tesla or SI
magnetization belongs at an explicitly converted API boundary.

The production radial XC call is the existing SCF path: `VXC0SP` supplies
`rho_down` and `rho_up` to `XCPOT`, which returns `Vxc_down` and `Vxc_up` in
the corresponding SCF ordering.  The response provider records

```text
Vxc_scalar = (Vxc_up + Vxc_down)/2
B_xc_energy = (Vxc_up - Vxc_down)/2 .
```

The constraining/FSM field is added after the XC call and is not XC response
data.  The legacy site-scalar radial projection uses

```text
K_perp,circular(site) = integral[B_xc_energy(r) m(r) dr]
                         / (2 M_site^2)
```

when the radial spin population and the response-projector population are
both available.  This is the kernel normalization for unhalved `O_+/-`
vertices.  The Cartesian x/y kernel is exactly twice this circular scalar.

## 4. Source/measurement vertices versus the interaction kernel

These are separate objects:

| object | code representation | meaning |
| --- | --- | --- |
| measurement/transition vertex | `response_channel` and `response_operator` | `P_site O_A`, used in `<n|A|m>` |
| external source vertex | `external_source_operator` | `dH_ext/db_A`; circular coordinates carry the explicit `1/2` |
| interaction kernel | `xc_response_kernel_provider` | derivative of the self-consistent Hxc field in the same response basis |

The eigenpair `chi_KS` builder contracts the ordered factors

```text
<n|A|m> <m|B|n>
```

and does not infer the second factor by conjugation.  For `chi^{+-}`, the
left/measurement channel is `O_+` and the right/source channel is `O_-`.
The common Dyson layer then uses the independently supplied kernel `K`:

```text
chi = chi_KS + chi_KS K chi
Xi  = chi_KS K
D   = I - chi_KS K .
```

The same XC derivative must not be put on both ends of `chi_KS` and then also
inserted as `K` in the Dyson equation.  A kernel-weighted transition factor
may be used as a computational composite for forming `chi_KS K`, but its
provenance must remain distinguishable from both the source vertex and the
measured observable.

## 5. Frequency, Green function, and Fourier conventions

The time convention is

```text
delta W(t) = Re[delta W(omega) exp(-i omega t)]
chi^R_AB(t) = -i theta(t) <[A(t), B(0)]> .
```

The code’s energy variable `omega` is in Rydberg.  The dynamic retarded
denominator is

```text
omega + epsilon_n - epsilon_m + i eta,    eta > 0,
```

and is implemented by `tddft_retarded_denominator`.  The advanced denominator
uses `-i eta`.  A physical angular frequency must first be converted to an
energy (`hbar omega_SI`) and then to Rydberg.

For a periodic real-space quantity the spatial convention is

```text
f(q) = sum_R exp(-i q.R) f(R) .
```

This is the convention required by the native real-space response path; the
electronic Green function does not need to be transformed to k space at
runtime.

## 6. Circular channels and the analytic two-level fixture

The response labels are operator ordered, `chi^{AB} = chi_{O_A O_B}`.  In a
collinear split one-orbital fixture with

```text
H = e0 sigma_0 + b sigma_z
epsilon_up   = e0 + b
epsilon_down = e0 - b
Delta = epsilon_down - epsilon_up > 0,
```

the positive excitation is an occupied-up to empty-down transition.  With
`M = f_up - f_down`, the production convention gives

```text
chi^{+-}_R(omega) = 4 M / (omega - Delta + i eta)
chi^{-+}_R(omega) = -4 M / (omega + Delta + i eta) .
```

Therefore:

```text
chi^{+-}_R(omega) = [chi^{-+}_R(-omega)]^*
chi^{+-}_A(omega) = [chi^{+-}_R(omega)]^*
```

for this real collinear fixture.  `chi^{+-}` contains the up-to-down positive
transition; `chi^{-+}` contains the reversed down-to-up transition.  They are
not interchangeable, especially in ferrimagnets and half-metals.

At positive excitation energy `Im chi^{+-}_R < 0`, so the scalar KS/Stoner
spectral weight is

```text
S^{+-}(omega) = -Im chi^{+-}_R(omega) / pi .
```

The integrated commutator weight for the unhalved operators is

```text
integral d omega S^{+-}(omega) = 4 M
```

when the pole is integrated over the full real axis.  The finite-grid fixture
uses a sufficiently broad interval and checks this normalization.  A factor
of two or a reversed frequency sign changes this test.

## 7. Static Ward/Goldstone identity

The static response uses the real Fermi divided difference

```text
(f_n - f_m)/(epsilon_n - epsilon_m)
```

with the analytic Fermi derivative for equal energies.  It has no dynamical
`eta`; the dynamic broadening is not allowed to define the static Ward limit.

For a collinear SOC-free ground state, in the same response basis,

```text
B_xc = chi_KS(q=0, omega=0)^(-1) M_GS = K M_GS
D M_GS = [I - chi_KS K] M_GS = 0 .
```

On the two-level fixture, `chi^{+-}(0) = -4M/Delta` and
`K_perp,circular = b/(2M) = -Delta/(4M)`, hence
`chi^{+-}(0) K_perp,circular = 1`.  This is an identity check, not an
empirical energy shift.  Any controlled Ward restoration must report its
magnitude and must not hide a failure of the underlying response
discretization.

## 8. Capability boundary

The conventions above do not silently extend to generalized overlap/HOH,
GBT, CCOR, Hubbard response, SOC, or general noncollinear response.  Those
paths must fail early until metric-consistent source vertices, kernels, and
symmetry identities have separate derivations and tests.
