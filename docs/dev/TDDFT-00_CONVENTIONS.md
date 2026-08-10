# TDDFT-00: Response conventions and XC-kernel interface

This is the executable convention contract for LR-TDDFT.  It covers response
operators and XC provenance only; it intentionally does not implement
`chi_KS`, a TDDFT Dyson solve, Goldstone correction, SOC response, GBT changes,
or longitudinal response.

## Hamiltonian and spin basis

The normal electronic basis is spin-major:

```
(orbital_1 up, ..., orbital_N up, orbital_1 down, ..., orbital_N down).
```

`basis_mod` sets `spin_off = norb`; `hamiltonian_build.f90` writes

```
H(up,up)     = H0 + Hz
H(down,down) = H0 - Hz
H(up,down)   = Hx - i Hy
H(down,up)   = Hx + i Hy.
```

Thus the *assembled electronic block* has the algebraic form
`H = H0 sigma_0 + Hx sigma_x + Hy sigma_y + Hz sigma_z`.  Its magnetic part
`hamiltonian%hxc` is filled by copying the full `Hvec.sigma` block.  It has
LMTO representation-dependent matrix and hopping structure; it is not an
identified local XC field and is never used as the TDDFT kernel input.

`response_basis_mod` builds those matrices in precisely this ordering.  Its
`response_operator(PLUS)` is `sigma_x + i sigma_y`, because the response
variable is defined as `m_plus = m_x + i*m_y`.  The conventional ladder
operator `sigma_plus = (sigma_x + i*sigma_y)/2` is separately available through
`ladder_operator(PLUS)`.  The same statements hold with `-` signs for MINUS.
This distinction is pinned by `UnitResponseConventions`.

## Density, moment, and units

`math_mod%rho2nm` defines the electronic density components directly from the
2x2 density matrix:

```
n  = rho_upup + rho_downdown
mx = rho_updown + rho_downup
my = Im(rho_downup - rho_updown)
mz = rho_upup - rho_downdown.
```

They are spin-population/electron-number components.  `bands%calculate_magnetic_moments`
integrates the same spin-resolved DOS difference and stores `potential%mx,my,mz,mtot`;
no `mu_B` factor is applied in that electronic response path.  Existing user
documentation calls the displayed number a spin magnetic moment in Bohr
magnetons, which is the usual numerical identification for an electron spin
population with `g=2`, but that conversion is not an internal response-kernel
normalization.

The internal energy unit is Rydberg: `math_mod%ry2ev` and `ry2joule` are used by
the code, and the libXC wrapper explicitly converts Hartree outputs by two.
`ry2tesla` exists only as a conversion constant.  The SCF XC potential and the
Hamiltonian coefficients above are energy-valued Rydbergs, not Tesla.  There is
no `mu_B` in `VXC0SP` or the Hamiltonian assembly.  Therefore no `mu_B` may be
introduced into a TDDFT kernel unless its input/output response variables have
first been converted from this spin-population convention to a physical
magnetization convention.

## Ground-state XC provenance

The production radial SCF path is

```
self%VXC0SP -> xc_obj%XCPOT(rho_down, rho_up, rho_total, ...) -> VXC1/VXC2.
```

`VXC0SP` passes `RHO1 = rho_down` and `RHO2 = rho_up`, and calls XCPOT with
output arguments in the order `VXC2, VXC1`.  Consequently `VXC1` is the up
potential and `VXC2` is the down potential.  It adds `+B_fsm/-B_fsm` only after
the XC call; this constraining field is not XC.

For a radial sample, `xc_response_kernel_mod%evaluate_ground_state_xc_sample`
calls that exact `XCPOT` route and stores

```
Vxc_scalar = (Vxc_up + Vxc_down)/2
Bxc_energy = (Vxc_up - Vxc_down)/2.
```

`Bxc_energy` is deliberately named as an energy coefficient, not a magnetic
field in Tesla.  `xc_response_kernel_provider` retains this provenance plus a
site spin population and future slots for `dVxc/dn`, `dVxc/dm`, `dBxc/dn`,
`dBxc/dm`, and `K_perp`.  It does not calculate `K_perp` from a Hamiltonian
block or an untraced site-moment ratio.  `self%VXC0SP` now accumulates its
existing radial SCF quadrature, using the returned XC potentials *before*
the constraining field is added.  The site projection is

```
K_perp(site) = integral Bxc_energy(r) m(r) dr / (2 M_site^2),
```

which is the **legacy site-scalar projection**. It is retained only as an
auditable comparison path: projection and multiplication by spatially/orbitally
varying `Bxc` do not commute, so it is not a Ward-consistent material kernel.
The frozen pair-potential replacement is documented in
`TDDFT_WARD_REPAIR_BLUEPRINT.md`. The response variables are `m+/- = mx +/- i my`, while the
Hamiltonian couples as `(delta B+ m- + delta B- m+)/2`; this supplies the
explicit factor one half.  No `mu_B` factor is introduced.  The provider retains the radial spin
population separately from `potential%mtot`, the latter supplying `M_site`.
SCF MPI
ownership is reduced so every rank receives the same complete provider.

The stored name is `k_perp_circular`: it is only for unhalved circular
`sigma_x +/- i sigma_y` vertices.  `circular_transverse_kernel` returns that
quantity, while `cartesian_transverse_kernel` returns exactly twice it for the
Cartesian `sigma_x`/`sigma_y` block.  A caller must select one of these typed
accessors; it may not reuse the circular scalar in Cartesian response.

`self%refresh_xc_response_kernel()` lets a caller that retains an SCF `self`
object refresh this provider through the existing atomic SCF/VXC path.  It has
no Hamiltonian-derived fallback.

## Retarded convention

No existing source module defines a dynamical TDDFT Fourier transform.  The
frozen TDDFT convention is therefore

```
delta W(t) = Re[delta W(omega) exp(-i omega t)]
chi^R(t)  = -i theta(t) <[A(t), B(0)]>.
```

Future `chi_KS` code must use energy-valued `omega` in Rydberg and the retarded
denominator `omega + epsilon_n - epsilon_m + i*eta`, with `eta > 0`.  Physical
angular frequency requires `omega_energy = hbar*omega_SI` converted to Rydberg.
This contract does not infer a time convention from static SCF code.

## Static Ward limit and response provenance

The Gamma Ward diagnostic does **not** call the dynamical response at
`omega=0` with its finite broadening.  Its dedicated eigenpair solver evaluates

```
(f(e_n) - f(e_m)) / (e_n - e_m)
```

as a real divided difference.  For equal or numerically near-degenerate
energies it uses the analytic derivative of the identical finite-temperature
Fermi function, `f'(e) = -f(e)*(1-f(e))/kT`.  Thus the static sign is negative,
there is no `i*eta`, and the static Xi imaginary norm is expected to be
numerical noise.  The dynamic `eta` remains a spectrum-resolution control and
cannot enter `Xi_static`.

Unless explicitly supplied in `&tddft`, response temperature is inherited from
the reciprocal ground-state object.  Once the complete response-mesh
eigenpairs exist, `auto_find_fermi=.true.` resolves the inherited chemical
potential on that same mesh against the ground-state target electron count.
This prevents a stale SCF-mesh Fermi seed from changing the response particle
number on a refined mesh.  An explicit `&tddft fermi_level` is never silently
replaced: the driver records both values and terminates when its recomputed
count differs by more than `1e-8 * max(1, N)`.  The transverse response vector
is the signed occupied `P_site sigma_z` population, not a site-moment
magnitude.

## Explicitly unresolved

- Derivatives of the legacy LDA/GGA and optional libXC functional paths are not
  yet evaluated.  The provider exposes their typed slots but marks each absent
  until a future evaluator populates it.
- `channel='full'` is capability-gated in production.  It is rejected unless
  the selected XC route supplies `dVxc/dn`, `dVxc/dm`, `dBxc/dn`,
  `dBxc/dm`, the local frame, and circular transverse value, and explicitly
  marks its derivative evaluator validated.  The present SCF XC routes do not
  make that claim.
- `chi0_backend='green'` currently means an eigenpair-resolvent reference, not
  an adapter to a native RS Green-function provider.  Longitudinal response is
  unavailable pending a WR-04 real-static-limit calibration; it is not LLB
  ready.
- The legacy SCF radial path currently calls `XCPOT` directly.  `XCPOT_hybrid`
  (the libXC wrapper) is not called by `VXC0SP`; a future libXC-response change
  must first reconcile that ground-state path rather than assume the wrapper
  was active.
