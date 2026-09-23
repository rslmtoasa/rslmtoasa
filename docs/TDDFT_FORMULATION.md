# TDDFT response formulation and production path

This document records the accepted representation and the native local-rotation
dynamics path. The physical observable is `E(q)`, obtained from a pole of the
frequency-dependent rotation-coordinate response `K^R(q,omega)^-1`. This is
not the full spatial spin susceptibility.

## Implementation status at d0a223e

```text
DRESP projection semantics       = CLOSED
second-order H adapter           = CLOSED
second-order torque/contact      = CLOSED
q=0 rotation Hessian             = CLOSED
native Turek static reference    = AVAILABLE
native rotation dynamics         = READY FOR IMPLEMENTATION
```

Finite-H versus native Turek finite-q curvature is a diagnostic comparison,
not a normalization identity. Native Turek remains the independent
adiabatic/MFT reference. The dynamic implementation uses the certified
second-order local-rotation Hamiltonian derivatives from `lr_kl_hessian`.

This is the baseline status at d0a223e. The completed implementation and first
bounded Fe campaign are recorded in
[`NATIVE_ROTATION_DYNAMICS.md`](NATIVE_ROTATION_DYNAMICS.md).

## 1. What DRESP-12 established

The accepted bcc-Fe finite-LMTO state closes the strict ASA spherical
(`L=0`) rigid accounting:

```text
ASA L0 rigid accounting residual ~= 1.45e-14
full nonspherical master norm  ~= 1.52084e-1
nonspherical L4 fraction       = 1.0
classification                 = ASA_L0_RIGID_RESPONSE_CLOSED
model boundary                 = NONSPHERICAL_RESPONSE_ON_SPHERICAL_ASA_GROUND_STATE
```

The `L=4` remainder is a model boundary, not a basis-response error,
Goldstone failure, or numerical defect. DRESP-12 closes the Ward-focused
campaign:

```text
Ward-focused campaign = CLOSED
Goldstone correction  = OFF
NEXT                  = NATIVE_ROTATION_DYNAMICS
```

The detailed static evidence remains in
[`DRESP_12_FINITE_LMTO_COVARIANCE.md`](DRESP_12_FINITE_LMTO_COVARIANCE.md).

## 2. Ground-state conjugate density and XC field

The live `VXC0SP` snapshot captures `n_up`, `n_down`, `vxc_up`, and `vxc_down`
from the same XC evaluation. The field used by the transverse strict-ASA
kernel is

```text
m_xc(r) = n_up(r) - n_down(r)
B_xc(r) = 0.5 * (vxc_up(r) - vxc_down(r))
```

For the accepted Fe state the provenance is legacy RS-LMTO Barth-Hedin
(`TXC=1`, spin-polarized LDA), and `constraining_field_ry = 0`. Thus the
stored `B_xc` is not the external/constraining field. The production
convention is explicit: `B_xc` is the XC field only and the ALSDA kernel is a
functional derivative of XC, never a ratio involving the external field.

`m_xc` includes the captured spherical SCF valence plus frozen-core density;
`P3` is the accepted reciprocal valence large-component Pauli polynomial.
Their profile audit is made in the common weighted radial coordinate
`4*pi*r^2*m_xc`, without fitting either profile to the other. The DRESP-12
artifact reports the weighted relative profile difference, integrated
values, signed integrated difference, sign convention, core bookkeeping,
and the physical scalar-relativistic spin comparison.

Accepted-Fe audit:

```text
XC functional / provenance       = Barth-Hedin / legacy RS-LMTO, TXC=1
integrated n_up - n_down         = 2.0000074695332746
integrated m_xc (common metric)  = 1.9989134964272697
integrated P3                    = 1.9984971713843327
integrated physical SR spin      = 1.9983526906817501
weighted profile difference      = 1.6000607534843933e-2
integrated m_xc - P3             = 4.1632504293698247e-4
integrated core m_xc             = -1.8124023342277402e-5
constraining field               = 0 Ry
```

The sign check is `PASS_UP_MINUS_DOWN`. The small core value is retained in
the bookkeeping rather than silently folded into the valence `P3` comparison.

For this LDA state the strict-ASA candidate is

```text
K_xc^ASA(r) = B_xc(r) / m_xc(r)
```

`Bxc/P3` remains historical representation-diagnostic machinery from the
DRESP campaign; it is not silently promoted to the physical kernel. If an
accepted state uses GGA, the complete transverse kernel is not this pointwise
ratio: gradient-dependent and noncollinear functional derivatives require a
separate derivation.

## 3. Native finite-H local-rotation dynamics

The dynamical coordinates are the two local transverse rotations
`theta_ix, theta_iy` for every magnetic site. Their retarded response uses the
accepted second-order `H=B-QB+E_nu` Hamiltonian and the certified finite-q
first and mixed derivatives from `lr_kl_hessian`. The vertex orientation is
`<m,k+q|T_A(q)|n,k>`; the finite-frequency bubble retains `Pi_AB` and `Pi_BA`
independently. The occupied mixed derivative is the frequency-independent
contact term.

The exact zero-frequency branch uses the finite-temperature Fermi divided
difference, including coincident energies. Its Cartesian static reduction is
the same-q index-symmetric part and closes against the force-theorem Hessian.
Finite-H versus native Turek finite-q curvature is a diagnostic comparison,
not a normalization identity. Native Turek remains the independent
adiabatic/MFT reference. No fitted Stoner parameter, scalar exchange
splitting, or projected spin operator replaces the accepted Hamiltonian
vertices.

This rigid local-rotation representation excludes arbitrary intra-atomic
radial/angular transverse deformation and longitudinal dynamics. It does not
claim equivalence to full spatial TDDFT. The global q=0 rotation needs no
Goldstone correction; finite-frequency electron-hole damping remains present.
## 4. Strict-ASA `L=0` physical-field TDDFT

Degrees of freedom are spherical but radially resolved transverse
magnetizations, with

```text
chi_KS^L0(r,r',q,omega)
```

and the XC-conjugate spherical kernel above. This route permits radial
intra-atomic freedom, different `s/p/d` radial amplitudes, and richer
Stoner/electron-hole structure than rigid site rotations.

Its frozen finite-LMTO fixed-basis response has a finite spherical Goldstone
defect before covariant/contact accounting. The relevant future quantity is
`||P0 r_fixed|| / ||m_G||`; the historical full-space `0.2363` is not the
`L=0` defect. Before any BES/Halle decision, the future `D_L0(omega)` scan at
`q=0` must locate the actual spurious pole/gap in Ry and meV.

The direct accepted-state mode geometry is diagnostic only:

```text
cosine(P0 Dm_G, m_G)          = 0.9908090998547531
parallel coefficient lambda0  = 0.09972593836402069
orthogonal fraction           = 0.13526761491581546
compact L0 DmG reconstruction = 1.3409439743737631e-7
```

The direct response-space `raw_residual` is decomposed into the component
outside the finite compact product span and the projected denominator mismatch.
Only the projected raw comparison and compact coefficient residual are DRESP-12
denominator gates. The projection-accounting residual checks the full raw norm
identity using the actual weighted inner product, including the cross term; it
does not assume the compact/raw maps are orthogonal in raw space. The full-space reconstruction residual is
reported only as a diagnostic because it contains the established L4 model
boundary.

## 5. Full-spatial `L>0` ALSDA

Degrees of freedom are intra-atomic spatial transverse responses with
`L=0...Lmax`. The present ASA ground state is stationary only in the
spherical radial-density space. A response that generates finite `L=4`
density therefore has no matching nonspherical ground-state XC field.

```text
STATUS = RESEARCH / NONSPHERICAL_GROUND_STATE_REQUIRED
```

This is not classified as a coding bug, and no nonspherical ground-state
functional is implemented here.

The `norm_by_l` check is retained as a partition/implementation consistency
check under the existing diagonal response metric. It is not an independent
proof of physical angular orthogonality.

## 6. Independent projected GSR route

The existing projected Jülich/GSR route remains separate:

```text
site-projected chi0 + sum-rule interaction
```

It is an independently normalized intermediate model and cross-check. It is
not merged with or substituted for the native KL/Turek route.

## 7. Production-state requirements

The DRESP `4x4x4` fixture remains certification-only. The first physics
state must be a fresh self-consistent bcc-Fe reciprocal state with:

- automatic electron-number Fermi level;
- no pinned historical `EF`;
- zero constraining field;
- second-order k-space Hamiltonian;
- accepted spin-polarized ground state.

Convergence is planned independently in k mesh, temperature/smearing,
frequency broadening `eta`, q, and frequency resolution. The historical
`eta = 0.04 Ry` is not suitable for resolving low-energy Fe magnons.

## 8. Critical production path

The production path is the frequency-dependent native local-rotation
response built from the accepted second-order Hamiltonian derivatives. The
native Turek/LKAG calculation is retained as an independent adiabatic/MFT
reference. No Goldstone patch or kernel rescaling is used. The strict-ASA
`L=0` radial TDDFT route is outside this implementation.

## 9. First bounded dynamical demonstration

The first physical state is a fresh cubic bcc-Fe reciprocal SCF state at 300 K,
with second-order `ham_only`, auto-found Fermi level, SOC and CCOR off, and no
constraining field. The 3x5x7 adapter fixture is certification-only. The first
q set contains Gamma and three small points along one direction. For each
finite q, the dynamic pole is compared with the second-order finite-H and
native Turek adiabatic estimates; neither estimate is substituted for the pole.
The first 12x12x12 campaign did not produce three mutually consistent pole
diagnostics, and its mesh sensitivity check changed the small-q curvature
substantially. No `E(q)` or stiffness fit is reported; Fe material convergence
remains open. The complete kernel, static reduction, Berry, circular-basis,
covariance, and pole contracts are recorded in
[`NATIVE_ROTATION_DYNAMICS.md`](NATIVE_ROTATION_DYNAMICS.md).

## 10. Current scope boundary

No Goldstone correction, strict-ASA radial `L=0` dynamics, Sternheimer
response, or nonspherical ground-state XC is part of this implementation.
