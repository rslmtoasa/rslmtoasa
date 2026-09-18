# TG-FZ-R4 — Finite-H / Turek Vertex Bridge

Status: **BLOCKED — no unfitted fixed-z identity between the live finite-H
torque and either tested Turek vertex**.

This audit is restricted to the existing two-site spd fixture and the three
fixed complex energies used by TG-FZ-R1/R2. It does not enter contour or
energy integration, native J_ij/J(q), Fe-shell comparisons, or DRESP-04.
The certified R1/R2 results remain unchanged.

The contour machinery in lr_kl_contour.f90 is outside this gate: R4 supplies
the three complex energies directly and calls no contour or integration
routine.

## Live B/Q/E_nu mapping

The source trace gives the following objects.

| object | live construction |
|---|---|
| B | lr_kl_hessian.f90:951-972: Fourier sum of the spinor bond blocks built from lmto_bond_value; these are the ee/h blocks. |
| Q | lr_kl_hessian.f90:967-970: Fourier sum of bond * obar_target; this is the eeo=ee*obarm_target product, not potential%qpar. |
| E_nu | lr_kl_hessian.f90:1109-1120: onsite spinor assembly from enu0/enu1; production build_enim obtains the same coefficient from cx-cex (hamiltonian_build.f90:1158-1190). |
| B_i | lr_kl_hessian.f90:770-802: endpoint moment derivatives of the bond blocks. |
| Q_i | The same routine differentiates both the bond factor and the target obar factor, including the endpoint product rule. |
| E_nu,i | lr_kl_hessian.f90:804-815: onsite derivative of the enu1 spin coefficient. |

The live second-order Hamiltonian is assembled as

    H = B + E_nu - Q B

and the complete fixed-q=0 rotation API added for this audit returns

    T_i = B_i + E_nu,i - Q_i B - Q B_i

(lr_kl_hessian.f90:246-284). The E_nu,i term is retained because it is
present in the actual H(theta) finite-difference oracle; omitting it would
not be the derivative of the live assembled matrix in this fixture.

The production spinor convention is explicit in
hamiltonian_build.f90:1254-1265 and lmto_magnetic_tangent.f90:118-136:

    H_uu = H4 + H3       H_dd = H4 - H3
    H_ud = H1 - i H2     H_du = H1 + i H2

The R4 fixture stores sites as (site1 up, site1 down, site2 up, site2
down). The native R1/R2 path remains in its live global ordering (up site1,
up site2, down site1, down site2) and is not reused to construct the
finite-H vertex. A positive y rotation about the collinear +z moment gives
the displayed ud/du blocks directly; no fitted sign, i, or 1/2 factor is
introduced.

potential%qpar is the orthogonal screening (gamma) input consumed by
symbolic_atom.f90:270-311 to form qi, dele, and the transformed obar=Y. It
is therefore algebraically related to the target obar, but it is not the
matrix Q in H=B-QB+E_nu. The latter is the reciprocal eeo product
constructed from obarm (hamiltonian_build.f90:1121-1156, 1289-1308). The R4
fixture expands the post-predls l-channel center_band, width_band,
shifted_band, and obar arrays into orbital channels so this Q-dependent path
is nonzero.

## Tested vertex chains

The historical chain is retained exactly:

    DeltaP_gamma  ->  d_raw = W_up DeltaP_gamma W_down

The common-alpha chain is the certified R2 transform:

    DeltaP_alpha -> tildeDeltaP_gamma
                  = R_down DeltaP_alpha R_up
                  -> d_tilde_alpha = W_up tildeDeltaP_gamma W_down

The independent finite-H chain is:

    H(theta) -> T_i = B_i + E_nu,i - Q_i B - Q B_i
             -> T_i^ud / T_i^du

The source-level product rule closes internally. The maximum local spin-flip
residual over all three energies is 3.469447e-18:

    T_i^flip = B_i^flip + (-Q_i B)^flip + (-Q B_i)^flip + E_nu,i^flip.

## Finite-H decomposition ledger

The decomposition norms are the maximum of site 1 (ud) and site 2 (du)
blocks. They are independent of z, as expected for a fixed Hamiltonian.

| term | matrix norm |
|---|---:|
| ||B_i^flip||_F | 1.631283e-01 |
| ||Q||_F | 1.809641e+00 |
| ||(-Q_i B)^flip||_F | 1.006105e-02 |
| ||(-Q B_i)^flip||_F | 1.157032e-02 |
| ||E_nu,i^flip||_F | 5.950826e-02 |
| ||T_i^flip||_F | 1.854430e-01 |

## Fixed-z numerical ledger

The vertex residual columns are maximum elementwise absolute residuals over
the site-1 ud and site-2 du blocks. Contractions are the ordered fixed-z
values; the actual-T contraction uses the corresponding spin-flip blocks of
the direct finite-H resolvent.

| z | ||d_raw||_F | ||delta_d_screen||_F | ||d_tilde_alpha||_F | ||T_i^flip-d_raw||_max | ||T_i^flip-d_tilde_alpha||_max | Q-term vs screen residual |
|---|---:|---:|---:|---:|---:|---:|
| -0.91+0.83i | 4.443344e-01 | 1.212928e+00 | 7.886845e-01 | 2.021129e-01 | 3.450894e-01 | 5.429024e-01 |
| -0.17+0.04i | 3.264915e-01 | 6.411666e-03 | 3.205444e-01 | 2.247992e-01 | 2.221090e-01 | 8.418373e-03 |
| 0.62+0.31i | 7.602390e-01 | 6.296638e-01 | 1.334363e+00 | 4.174498e-01 | 6.684216e-01 | 2.780459e-01 |

The maximum residuals are therefore:

    T vs raw d                         4.174498e-01
    T vs transformed-alpha d          6.684216e-01
    Q-term vs R2 delta_d_screen        5.429024e-01
    minimum tested candidate residual  2.021129e-01

Neither candidate closes, and the R2 screening correction is not supplied by
the local (-Q_iB-QB_i) block in this fixture. The exact remaining term is
not guessed or fitted: for either candidate it is the live source expression

    M_i(candidate) = B_i^flip - candidate_i
                    + (-Q_i B)^flip + (-Q B_i)^flip
                    + E_nu,i^flip,

with the product-rule residual above proving that this is precisely
T_i^flip-candidate_i.

The fixed-z contraction ledger is:

| z | raw-gamma | alpha-direct | transformed-gamma | actual-T ud |
|---|---:|---:|---:|---:|
| -0.91+0.83i | -5.380807e-06 | -3.582616e-06 | -3.582616e-06 | 2.466684e-05 |
| -0.17+0.04i | -2.202285e+00 | -2.125794e+00 | -2.125794e+00 | -6.902849e-01 |
| 0.62+0.31i | -3.412346e-04 | -1.377152e-03 | -1.377152e-03 | 4.995552e-05 |

The alpha-direct/transformed-gamma identity still closes to 4.440892e-16;
this is the established R2 covariance chain, not a closure to the finite-H
torque.

## Independent finite-difference oracle

For each of sites 1 and 2, the actual fixture moments were rotated about y,
the complete finite-H matrix was assembled at +epsilon and -epsilon, and the
central difference was compared with the analytic T_i. The two site
residuals are identical to shown precision.

| epsilon | max site residual |
|---:|---:|
| 1e-2 | 1.337936e-06 |
| 1e-3 | 1.337943e-08 |
| 1e-4 | 1.337943e-10 |
| 1e-5 | 1.337958e-12 |
| 1e-6 | 1.337819e-14 |
| 1e-7 | 1.526557e-16 |

The converged oracle gate (epsilon <= 1e-4) is 1.337943e-10 and passes. The
coarse-step values are retained to show the expected O(epsilon^2) truncation
window.

## Contact term

The independently assembled distinct-site mixed derivative is nonzero:

    ||C_12||_F = 7.906970e-04.

The complete fixed-z curvature split, using the full finite-H T_i, T_j, and
C_12, is:

The scalar entries below were regenerated in TG-FZ-R8R with the independent
direct `Tr(H_ij G)` contraction. Any pre-R8R API printout that showed a zero
contact for these distinct-site matrices is superseded; the R4 vertex and
matrix identities above are unchanged.

| z | TT contribution (-Im/pi) | contact contribution (-Im/pi) | complete (-Im/pi) |
|---|---:|---:|---:|
| -0.91+0.83i | -7.718933e-05 | 1.334771e-04 | 5.628774e-05 |
| -0.17+0.04i | 3.295878e-01 | 2.860561e-03 | 3.324483e-01 |
| 0.62+0.31i | -1.639025e-04 | 2.040903e-04 | 4.018783e-05 |

The TT/contact split is internally exact to the reported 0.000000e+00
residual. The contact term cannot be dropped or inferred to vanish from the
fact that the two sites are distinct.

## Disposition

Proved in R4:

* the source-level finite-H B,Q,E_nu product rule and its rotation derivative
  are implemented and independently finite-difference verified;
* the Q-dependent pieces are distinct, nonzero, and do not equal the R2
  coefficient-space screening correction in this fixture;
* the common-alpha direct/transformed-gamma covariance remains exact; and
* the distinct-site C_12 contact contribution is nonzero at fixed z.

Unresolved:

* neither d_raw nor d_tilde_alpha is the actual finite-H torque vertex;
* no source-derived representation map from the native Turek vertex to the
  complete finite-H T_i has been established.

Therefore the verdict is **BLOCKED**, not a fitted PASS-B/PASS-C. The next
step is to derive the missing representation map from the accepted production
Hamiltonian and its endpoint normalization, if that bridge is required. Do
not proceed to energy integration or native J/J(q) on the basis of this audit.

## Verification and frozen invariant

The focused test and nearby DRESP/finite-H selection passed:

    cmake --build build -j2
    ctest -V -R UnitDresp03tgNativeFixedZ --output-on-failure
    ctest -R 'UnitDresp03|UnitLrKlStaticBridge' --output-on-failure

source/exchange.f90 was not modified. Its SHA-256 remains:

    6e9ae7da6af46367aaf18a4fef37bdd4a9a0d6f69fe8e08b589e5c5f278a70e6
