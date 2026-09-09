# Literature contracts

## REF-BES-2011
P. Buczek, A. Ernst, L. M. Sandratskii, Phys. Rev. B 84, 174418 (2011).

- **BES-01 Dyson:** `chi = chi_KS + chi_KS Kxc chi`.
- **BES-02 ALSDA transverse kernel:** Eq. (22), `Kxc(x) = -mu_B Bxc(x)/m(x)` in their convention.
- **BES-03 Goldstone identities:** Eq. (31), `|Bxc> = chi_KS^{-1}(0)|m_GS>` and `|Bxc> = Kxc|m_GS>`.
- **BES-04 numerical correction:** Eqs. (38–39), modify only the small Goldstone eigenvalue of `D=I-chi_KS Kxc`, reconstruct `D_corr`, derive corrected Kxc.
- **BES-05 spatial basis:** site × real spherical harmonic × radial Chebyshev basis; ASA replaces Voronoi cells with atomic spheres; full spatial kernel retained.
- **BES-06 finite-q mixed representation:** lattice Fourier phase `exp(-i q·R)` in Appendix A.

## REF-LCMM-2011
S. Lounis, A. T. Costa, R. B. Muniz, D. L. Mills, Phys. Rev. B 83, 035109 (2011).

- **LCMM-01 exact static sum rule:** `sum_j int dr' chi0^{ij}(r,r';0) B_eff^j(r') = m_z^i(r)`.
- **LCMM-02 ASA radial reduction:** angular reduction gives the published `4pi` form.
- **LCMM-03 local interaction:** `U^j(r)=B_eff^j(r)/(4pi m_z^j(r))` in their convention.
- **LCMM-04 sum-rule solution:** `Gamma U = m_z`, hence `U = Gamma^{-1} m_z`.
- **LCMM-05:** the construction is intended to avoid arbitrary Goldstone tuning.

## REF-EEB-2026
D. Eilmsteiner, A. Ernst, P. A. Buczek, arXiv:2603.03220v2 (2026).

- **EEB-01:** ASA/FCD explicitly allowed for the closed-lattice systems discussed.
- **EEB-02:** `chi(r,r') = sum_LL' chi_LL'^{ss'}(r_s,r_s') Y_L Y_L'`.
- **EEB-03:** direct Gaussian radial mesh; authors report improved Goldstone-sum-rule fulfillment.
- **EEB-04:** super-index `(density channel, site, angular, radial)`.
- **EEB-05:** local-frame kernel rotated globally with `Kxc = T Kxc_loc T^dagger`.
- **EEB-06:** q=0, omega=0 Dyson denominator contains rigid-rotation zero-energy deformation subspace.

This is a modern preprint feature target, not the sole foundation of basic TDDFT identities.

## REF-SAV-1998
S. Y. Savrasov, Phys. Rev. Lett. 81, 2570 (1998).

Variational Sternheimer TD linear response in a muffin-tin-orbital representation. Use only as a later independent LMTO-native benchmark.
