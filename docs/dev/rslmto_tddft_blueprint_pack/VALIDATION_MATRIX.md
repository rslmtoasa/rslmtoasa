# RS-LMTO-ASA TD-DFT validation matrix

The goal is not only to get plausible dispersions. Each tier tests a different structural property of the theory and implementation.

| Tier | System / fixture | Main purpose | Required comparisons |
|---|---|---|---|
| 0 | analytic 1- and 2-level models | signs, factors, retarded structure, contour integration | analytic result vs eigenpair vs G(k) vs G(R) where applicable |
| 0 | periodic electronic toy / tiny crystal | one-electron GF consistency | direct inverse vs Lehmann G(k); FT G(R)<->G(k) as a validation-only identity |
| 0 | static transverse fixture | Ward/Goldstone identity | `chi0(0) Bxc` vs ground-state m; `D m ~ 0` |
| 1 | bcc Fe, SOC=0 | simple itinerant FM release gate | all three chi0 backends; Goldstone; q^2; LKAG/frozen-magnon/GBT |
| 1 | fcc Ni, SOC=0 | difficult itinerant FM / previous failure | backend equivalence; near-Gamma q^2; Stoner overlap; broadening/convergence |
| 2 | FeRh FM | multi-sublattice FM | acoustic/optical modes, damping |
| 2 | FeRh AF | AF branch structure/chirality | literature trend + adiabatic reference |
| 2 | CrMnSb | ferrimagnetic chirality | separate `+-/-+` branches and asymmetric damping |
| 3 | Pd / near-magnetic metal | paramagnon optional target | enhanced vs KS response |
| 3 | Fe, Ni longitudinal | `(0,z)` metallic longitudinal response | Buczek et al. 2020 trends |
| 3 | AF FeSe longitudinal | magnetic-second-sound target | energy/FWHM/momentum trends |
| 4 | SOC adatom/film fixture | future relativistic response | relativistic sum rule, gap, linewidth |

## Core numerical convergence axes

For every production benchmark record:

- k mesh (eigenpair and k-GF paths);
- R cutoff / number of lattice shells (r-GF path);
- energy mesh or contour nodes;
- eta / complex offset gamma;
- electronic temperature/smearing used only for numerical integration;
- response-basis resolution (site vs site-orbital if used);
- q path in explicit reciprocal coordinates;
- frequency mesh and peak-fit window;
- Goldstone/Ward residual before any optional correction;
- magnitude and type of any optional numerical correction.

## Release gates

### Gate A — one-electron GF
- [x] `G_k_lehmann(z)` agrees with direct resolvent for representative complex z.
- [ ] native `G_R(z)` Fourier-transforms to the same periodic `G_k(z)` on a controlled fixture within truncation error.
- [x] retarded/advanced and Hermiticity identities hold.

### Gate B — chi0 backend equivalence
- [x] eigenpair and k-GF chi0 agree on small periodic fixtures.
- [ ] k-GF and r-GF chi0(q,w) agree after R convergence.
- [ ] direct and mixed-contour integration agree.

### Gate C — static physics
- [ ] spectral-weight sum rule is satisfied to converged numerical tolerance.
- [ ] `D(0,0)` has the magnetization as its Goldstone null vector in SOC=0.
- [ ] residual decreases systematically under numerical refinement.

### Gate D — Fe/Ni production
- [ ] gapless SOC=0 acoustic mode.
- [ ] q^2 dispersion near Gamma.
- [ ] stiffness consistent with independent adiabatic references.
- [ ] stable linewidth/Stoner-continuum behavior.
- [ ] all three backends agree within established uncertainty.

Only after Gate D may the transverse path be called production validated.
