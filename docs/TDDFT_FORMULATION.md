# TDDFT response formulation and production path

This document is the roadmap after the closed DRESP-12 Ward campaign. The
target observable is the transverse susceptibility

```text
chi+-(q,omega)
```

No dynamical susceptibility, Sternheimer/basis response, nonspherical
ground-state XC functional, or long frequency/q campaign is implemented by
this decision.

## Implementation status after the DRESP-12 handoff

```text
DRESP-12 handoff cleanup = RAW_L0_REGRESSION_OPEN
native rotation dynamics = BLOCKED — STATIC_NATIVE_NORMALIZATION_OPEN
Fe magnon poles          = NOT RUN
```

The required static gate is the spectral local-rotation torque/contact kernel
against the certified native Turek curvature on the same accepted state. The
available `exchange_q` spectral path only accepts first-order `ham_only` and
is not the required fresh second-order production route. A first-order 12³
bcc-Fe diagnostic at 300 K found finite-H versus native-Turek differences at
small off-mesh q; it did not close this gate because it is neither the
required Hamiltonian order nor a controlled q/mesh comparison. The requested
raw DRESP-12 reconstruction was measured at `4.2787e-5` on the 4³ smoke
fixture; the existing `1e-6` gate fails there, while the accepted 64-k raw
residual has not been rerun. No tolerance relaxation or dynamical rescaling is
applied. The native-dynamics implementation and pole campaign remain stopped
at their respective open gates.

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

## 3. Native KL/Turek dynamical formulation

Degrees of freedom are sitewise transverse orientations. The certified native
LMTO rotation supplies the vertex

```text
T_i = -i [G_i, H]
```

and the bare response is schematically `Pi_ij(q,omega)` from the accepted
reciprocal eigenstates/Green-function machinery and native site vertices.
This maps onto the existing native LMTO/Turek path-operator and
MFT-LKAG infrastructure documented in
[`DRESP_03TG_NATIVE_TUREK_GF.md`](DRESP_03TG_NATIVE_TUREK_GF.md) and the
native LKAG code path. The global rotation is exact by construction, so this
representation needs no Goldstone correction. Its finite-frequency
electron-hole structure remains, so Landau damping is not excluded.

The route excludes arbitrary intra-atomic radial/angular transverse
deformation and excludes longitudinal dynamics. It is therefore not claimed
to be equivalent to full spatial TDDFT.

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

The critical path is **NATIVE KL/TUREK DYNAMICAL RESPONSE**. It is the
native LMTO representation with a certified global rotation, reuses the
existing Turek/LKAG machinery, needs no arbitrary Goldstone patch, and gives
the direct route to physical magnon dispersion while retaining dynamical
electron-hole damping in the projected rotation channel.

The strict-ASA `L=0` physical-field route is retained as the second
production/research track for radial intra-atomic physics. This is a
formulation choice, not an equivalence claim between the two models.

## 9. First dynamical physics milestone

The first observable is the close-to-Gamma bcc-Fe magnon energy `E(q)`.
The validation ladder is:

1. `q=0` native Goldstone;
2. several small finite-q points;
3. identify the dominant low-energy transverse pole;
4. fit `E(q)=D q^2`;
5. establish k-mesh convergence of `D`;
6. establish eta/frequency-resolution convergence;
7. compare `D` with native/original MFT-LKAG/Turek stiffness and established
   bcc-Fe literature/experiment;
8. extend along symmetry lines only after those checks.

Equality to Bruno-renormalized finite-q exchange is not a requirement. The
long-wavelength stiffness should agree in the appropriate adiabatic limit.

## 10. Immediate next implementation milestone

Implement the native KL/Turek bare dynamical rotation response using the
certified reciprocal eigenstate/GF services and native site vertices. The
first implementation must expose `q=0`, finite-q, pole, damping, and
convergence diagnostics. It must not add a Goldstone correction, Sternheimer
response, or nonspherical ground-state XC.

No dynamics is implemented in the DRESP-12 formulation decision.
