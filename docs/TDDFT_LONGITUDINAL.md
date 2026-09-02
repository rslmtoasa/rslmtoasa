# Coupled charge–longitudinal TDDFT response

TDDFT-13 adds the collinear, scalar-relativistic charge/longitudinal-spin
sector.  For `N` response sites, both the bare Kohn–Sham susceptibility and
the Dyson kernel use the site-major basis

```text
(charge_1, m_z_1, charge_2, m_z_2, ..., charge_N, m_z_N)
```

The response dimension is therefore `2*N`; `longitudinal_index(site, 0)` and
`longitudinal_index(site, 3)` are the single source of truth for the charge
and `m_z` indices.  The existing eigenpair transition engine and K-space
Lehmann Green-function adapter are reused with this channel list, so the
charge–spin off-diagonal elements are accumulated by the same vertex and
occupation machinery as the transverse response.

## Kernel and units

The production enhancement is

```text
chi = chi_KS + chi_KS (f_H + f_xc) chi
```

`f_H` is copied from the projected `charge%amad` matrix in the existing
Rydberg/LMTO response metric and is inserted only at charge–charge indices.
It is not copied into either spin leg.  The local XC block is

```text
             [ d v_xc / d n   d v_xc / d m ]
f_xc(site) = [ d B_xc / d n   d B_xc / d m ]
```

with `n = n_up + n_down` and `m = n_up - n_down`.  During the response SCF
refresh these derivatives are obtained by symmetric finite differences of the
same `XCPOT_hybrid` route that produced the ground-state potentials.  The
radial values are projected with the existing ASA quadrature and normalized
against the response-site populations, keeping radial and response
projectors distinct.  This local derivative route is currently accepted only
for LDA/ALDA-capable selectors; legacy and native GGA routes fail explicitly
because their nonlocal radial `f_xc` terms are not yet derived.

## Backend scope and guards

Eigenpairs and the K-space Lehmann backend are supported.  The native
real-space Green-function backend remains explicitly rejected for this
sector until its multi-component real-space/Fourier source and static limit
are derived; it must not silently fall back to a different backend.  The
existing production boundary also rejects SOC, noncollinear states,
generalized-overlap Hamiltonians, HOH/GBT/CCOR/Hubbard corrections and
constrained external fields.

The charge–longitudinal block is not a rigid-spin rotation.  Consequently no
Goldstone/Ward constraint or transverse zero-mode correction is applied to it.
The old TDDFT-08 finite-field file controls remain readable for compatibility
with the static prototype unit API, but are not required or consumed by the
production TDDFT-13 route.

## Outputs and interpretation

With `output_xi` or `output_chi`, each q point writes
`*_qNNNNNN_longitudinal_dyson.dat`.  The file contains the ordinary
frequency-dependent `chi_KS`, `Xi`, enhanced `chi`, and loss products, plus
metadata identifying the `(charge,m_z)` layout, Hartree/XC provenance, and
the numerical role of `eta`.  `eta` is a numerical broadening, not a physical
linewidth.

These files are longitudinal susceptibilities.  Poles or loss peaks may be
examined as collective-mode candidates only after the usual determinant,
spectral-weight, and convergence checks.  No susceptibility fit is labeled
as a damping rate, `T_parallel`, or an LLB `alpha_parallel`; a later LLB
study must specify its dynamical equation, normalization, and dissipative
mapping independently.

## Validation status

The unit fixture verifies the site-major `(0,z)` layout, charge-only Hartree
coupling, explicit XC derivative capability, and the shared eigenpair
adapter.  The existing TDDFT-08 static-calibration and transverse regression
tests remain active.  Production Fe/Ni and AF-FeSe qualitative runs are a
separate material-validation gate and are not claimed by the analytic unit
fixture alone.
