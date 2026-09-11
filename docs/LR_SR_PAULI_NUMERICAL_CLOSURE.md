# LR-02N scalar-relativistic → Pauli numerical closure

## Evidence status

This document records the numerical closure requested by LR-02N. It measures
the difference between the accepted scalar-relativistic radial density and the
large-component Pauli projection for one converged magnetic bcc-Fe ground
state. It does not introduce a response kernel, susceptibility, Dyson
equation, or frequency-dependent formalism.

The evidence was generated on branch fable_v4, starting from exact HEAD
2793c1c47fd6d2669184ac53707d2186d6200919 (the working tree contained only
the focused LR-02N changes). The live run used:

    fixture       tests/scf/cases/bulk/bccFe
    control       nsp=1; scalar-relativistic collinear; no SOC
    SCF           nstep=16, conv_thr=2.0e-6, beta=0.01, Broyden
    k mesh        8 x 8 x 8 full Monkhorst-Pack mesh
    eigenproblem  reciprocal_mode=ham_only, kspace_ham_order=first
    occupations   Fermi-Dirac, T=300 K, EF=-0.069249685116850188 Ry
    k-weight sum  1.0000000000000000
    XC            legacy RS-LMTO Barth-Hedin
    additive terms none: no CCOR and no Hubbard terms

The EF is the accepted final ground-state Fermi level consumed by the
projection. The reciprocal eigenstates are built from the final accepted
potential parameters after SCF; the radial reference is not regenerated from
those parameters.

## Definitions and measures

For each site and local spin channel,

\[
n_\sigma^{SR}(r)
 = \frac{\mathrm{RHO}_\sigma(r)}{4\pi r^2}
\]

is the accepted LR-01 snapshot. The Pauli density is formed from occupied
reciprocal eigenvectors and the first-order large-component LMTO augmentation:

\[
U_{al\sigma}^{n\mathbf{k}}(r)
 = G_{al\sigma}(r)
 + \left(\epsilon_{n\mathbf{k}}-E_{al\sigma}^{\nu}\right)
   \dot G_{al\sigma}(r),
\]

\[
n_\sigma^P(r)
 = \frac{1}{4\pi r^2}
   \sum_{n\mathbf{k}}w_\mathbf{k}f_{n\mathbf{k}}
   \sum_{lm}|c_{alm\sigma}^{n\mathbf{k}}|^2
   |U_{al\sigma}^{n\mathbf{k}}(r)|^2 .
\]

The angular normalization is the certified LMTO spherical normalization: the
m-sum is diagonal for the spherical L=0 density and the radial weighted
quantity is divided by 4πr². No GFAC, lower component, or effective
square-root amplitude enters nP. The frozen occupied core is added separately
using its accepted large component so that the comparison covers the same
occupied core+valence charge as LR-01; this treatment is declared in every
output file.

The reported differences are

\[
\delta n_\sigma=n_\sigma^{SR}-n_\sigma^P,\qquad
m^{SR}=n_\uparrow^{SR}-n_\downarrow^{SR},
\]

\[
m^P=n_\uparrow^P-n_\downarrow^P,\qquad
\delta m=m^{SR}-m^P.
\]

The certified log-mesh measure is

\[
dr_i = a(r_i+b)\,d(\log r)_i,
\qquad
\int q(r)\,dr \equiv
\sum_i s_i\,a(r_i+b)\,q_i,
\]

where s_i is the production composite-Simpson weight. Integrated charge uses
d³r=4πr²dr. Norms use the same physical volume measure:

\[
\|q\|_2^2=\int 4\pi r^2 q(r)^2\,dr.
\]

The magnetically relevant region is objective: the shortest radial prefix
containing 90% of the integrated absolute accepted magnetization
\(\int 4\pi r^2|m^{SR}(r)|dr\).

## Fe pointwise comparison

The emitted file contains every radial point:
lr02n_pauli_projection_Fe_1.dat. Representative rows are:

| radial location | r (bohr) | delta n_up | delta n_down | delta m |
| --- | ---: | ---: | ---: | ---: |
| near core, row 2 | 2.75304689e-6 | 7.69686061e1 | 7.69312582e1 | 3.73478969e-2 |
| valence-region sample, row 248 | 1.89116664e-2 | 2.27179057e1 | 2.27092546e1 | 8.65112265e-3 |
| outer ASA point, row 495 | 2.66220000 | -8.83077856e-5 | -1.96855981e-4 | 1.08548195e-4 |

The large near-core pointwise values are density values; they are multiplied
by the vanishing 4πr² volume factor in the integrated comparison. They are not
replaced by a maximum-relative-error criterion.

## Fe integrated comparison

All quantities below are in electrons, with the spin difference in the native
positive up-minus-down convention (and therefore also in μB for the reported
moment):

| quantity | SR | Pauli | SR − Pauli |
| --- | ---: | ---: | ---: |
| N_up | 14.0519328301 | 14.0729582344 | -2.10254043e-2 |
| N_down | 11.9480745890 | 11.9242315409 | +2.38430481e-2 |
| total charge | 26.0000074191 | 25.9971897753 | +2.81764376e-3 |
| magnetization M | 2.1038582410 | 2.1487266935 | -4.48684524e-2 |

The SR values reproduce the integrated values in the accepted LR-01 snapshot
from the same run. The occupied valence coefficient sum is 7.9410931941,
split as 5.0145806311 (up) and 2.9265125630 (down). The large-component
augmented valence integrals are 5.0892450682 and 2.9405168361. This pair of
diagnostics is intentional: the raw coefficient sum verifies
occupation/eigenvector normalization, while the augmented integral includes
the accepted first-order large-component radial weight. It must not be
silently identified with the raw coefficient sum.

## Volume-weighted norms

| region | relative L2 magnetization | relative L2 total charge |
| --- | ---: | ---: |
| full ASA sphere | 1.43312057e-2 | 9.10157975e-3 |
| 90% absolute-magnetization prefix | 1.40811445e-2 | 9.10159791e-3 |

The objective magnetic prefix ends at radial point 469 of 495,
rmax=1.5826773439 bohr. The table reports the discrepancy; it does not assign
a pass threshold or call it negligible.

## Omitted scalar-relativistic radial contribution

The diagnostic term retained from LR-02R is

\[
\delta D_l^{SR}(r)=
\frac{l(l+1)}{[TMC(r)r]^2}G_1(r)^2+G_2(r)^2.
\]

Its radial integrals and ratios to the corresponding large-component norm were
emitted for every l and spin channel in
lr02n_pauli_projection_Fe_1.dat.summary. The Fe values are:

| l | omitted norm, up | ratio, up | omitted norm, down | ratio, down |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 5.10972962e-5 | 5.10999073e-5 | 5.03820928e-5 | 5.03846313e-5 |
| 1 | 5.72579484e-5 | 5.72612271e-5 | 5.85962308e-5 | 5.85996645e-5 |
| 2 | 1.85547664e-4 | 1.85582099e-4 | 1.48461215e-4 | 1.48483259e-4 |

These are diagnostic radial norm contributions only. They are not asserted to
equal the full density or arbitrary-response difference.

## Reproducibility and consistency checks

The live validation is
tests/validation/val23_lr_pauli_projection.py, registered as
Val23LrPauliProjection when regression tests are enabled. It verifies:

- the exact accepted LR-01 file is the SR reference;
- the logarithmic mesh and certified Simpson/volume measures;
- every emitted delta n, m, delta m, charge, and delta-charge identity;
- integrated SR, Pauli, and delta quantities against the summary;
- SR integrated charge against the LR-01 snapshot;
- raw occupied coefficient weight against the stated Fermi occupations;
- explicit large-component-only/no-GFAC/lower-component exclusion provenance;
- the objective 90% magnetic region and finite radial data.

Run the deterministic repeated execution check with:

    python3 tests/validation/val23_lr_pauli_projection.py \
      --binary build/bin/rslmto.x \
      --scratch-root /tmp/lr02n_validation_repeat \
      --repeat 2

No arbitrary physical discrepancy threshold is used. The numerical tolerances
in the validation script are only quadrature/identity checks. Changing --nk
changes the reciprocal sampling and should be treated as a separate numerical
convergence study, not as a changed definition of the projection.

## Ni disposition

fcc-Ni is deferred. This repository does not contain an existing accepted
LR-01 fcc-Ni radial-ground-state fixture comparable to the bcc-Fe case. Adding
and stabilizing a second physical SCF workflow would make Ni a separate
fixture-validation task rather than a narrow LR-02N closure. Fe is therefore
the mandatory physical result for this evidence.

## Claim level

This computation establishes only the magnitude and radial structure of the
scalar-relativistic-to-Pauli large-component ground-state projection
difference for the chosen converged bcc-Fe state. It does not establish
TD-DFT correctness, an ALSDA/XC kernel, Goldstone compliance, magnon
dispersion, or adequacy for heavier elements.
