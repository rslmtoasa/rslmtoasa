# A5.3 numeric charge-only U allowance registered before candidate comparison

This additive numerical preregistration binds the measured fine-grid spin-W1
response to the already registered charge-only chain rule. Candidate a/b/c
energy or derivative outcome rows have not been inspected for this reduction.
The original h=0.005/0.0025 W1 result remains separately preserved as unresolved;
the finer-grid result is a supplementary measurement under the already
authorized experiment and unchanged independent-box propagation.

For the fine-grid primary response at h=0.000625 Ry, the registered quantities
are `B=(0.00123034,-0.01042689) Ry/Ry`, with componentwise uncertainty
`(0.00010641,0.00015102)`, and
`R=[[-2.25713,1.70644],[1.69827,-2.33353]] electron/Ry`, with componentwise
uncertainty `[[0.03492,0.03770],[0.03548,0.03854]]`. The determinant of `R^T`
is 2.36908 and its full componentwise uncertainty-box corner range is
[2.07604,2.66213], so the solve remains identifiable throughout that box.
The nominal solution is `g=(0.00626264,0.00904798) Ry/electron`; its
first-order conservative propagated component uncertainty is
`(0.00116391,0.00116646) Ry/electron`.

The measured U-response input is the central h=0.05 eV direct-P derivative;
its uncertainty is the measured h=0.10 versus h=0.05 step difference. Direct
P_OUT q rows were independently checked against per-spin d-shell ql rows.
The mode-specific values and resulting charge-only allowances are:

| radial mode | `dq/dU` (electron/eV) | step uncertainty | nominal `N_ASA` (Ry/eV) | propagated first-order uncertainty | independent-box interval (Ry/eV) |
|---|---:|---:|---:|---:|---:|
| local | `(0.08327349,-0.08073297)` | `(0.00239066,0.00304498)` | `0.00125198` | `0.00023362` | `[0.00102115,0.00155303]` |
| spectral | `(0.10781323,-0.03129201)` | `(0.00635735,0.00392829)` | `0.00095832` | `0.00023734` | `[0.00074110,0.00124622]` |

The interval enclosure propagates the independent B, R and U-response boxes
through `R^T g=B` and `N_ASA=sum_s |g_s dq_s/dU|`; it is conservative and is not
a statistical confidence interval. The raw fine-grid summary and reduction
hashes, along with the U-endpoint records hashes, are recorded in
[`a5_spin_w1_finer_reduction.json`](a5_spin_w1_finer_reduction.json). The
historical coarse W1 box is unchanged and remains unresolved; it is not silently
replaced by this supplementary fine-grid box.

For the candidate decision, the registered base-functional tolerance is the
mode-specific interval above, combined with the separately measured numerical
uncertainty of the production fixed-P Hubbard derivative and the candidate U
finite difference. Candidate b is assessed against the production fixed-P
target plus this A5.3 allowance. Candidate c's predeclared diagnostic miss is
the measured `dt/dU` scale (approximately −0.1486 Ry/eV spectral and
−0.1419 Ry/eV local); its outcome remains blinded until the independent lead
reduction is checked. A-versus-b remains report-only unless separated beyond
the allowance. This registration supplies no adoption verdict by itself.

The J allowance is separate and uses `sum_s |g_s dq_s/dJ|` at J=0.5 and 0.9
eV, with h=0.05 versus h=0.10 response sensitivity and the same propagated
gradient uncertainty. Before inspecting J candidate energy rows, the
spectral and local charge-response measurements give:

| radial mode | J center | `dq/dJ` at h=.05 (electron/eV) | h=.10 minus h=.05 component magnitudes | nominal `N_ASA` (Ry/eV) | independent-box interval (Ry/eV) | fixed-P target `dE_Hub/dJ` (Ry/eV) |
|---|---:|---:|---:|---:|---:|---:|
| spectral | .5 eV | `(-.23319382,+.16540535)` | `(.03875380,.03706411)` | `.00295699` | `[.00204434,.00416579]` | `-.063047077414` |
| spectral | .9 eV | `(-.08876260,+.02910896)` | `(.00919943,.00907050)` | `.00081927` | `[.00057640,.00113991]` | `-.067713496800` |
| local | .5 eV | `(-.13141686,+.06715841)` | `(.00196893,.00359923)` | `.00143066` | `[.00118574,.00174694]` | `-.063520683241` |
| local | .9 eV | `(-.12709680,+.10879738)` | `(.02618505,.02782279)` | `0.00178036` | `[0.00117606,0.00258153]` | `-.064806830829` |

The interval independently propagates all componentwise B/R uncertainties
through `R^T g=B` and the measured h=.10-vs-.05 charge-response boxes; it is
not a statistical confidence interval. Spectral q endpoint records are from
`/private/tmp/a5-m22m23-remaining-ks-grid`, respectively
`u2_j.45/.55/.4/.6` and `u2_j.85/.95/.8/1`; local q endpoint records are from
`/private/tmp/a5-local-followup-ks` at the same centers. The local case
summary SHA256 is `e7463ccb1417e4e9acce0157b6ab1e6bf19d1d993055e0bf9fc72035cf5c9a4d`.
Fixed-P target P_OUT matrices use 50 direct matrix entries from each mode's
center state; local input/target hashes are in `a5_m23_local_j_allowance.json`
and spectral input/target hashes are in the earlier target evidence. These J
allowances were completed before inspecting local J candidate energy rows.
