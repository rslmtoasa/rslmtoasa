# TDVK-02R0 — LMTO product response-basis closure audit

**Verdict: PASS for algebraic representation closure.** The finite LMTO
radial-product space spans the current complete LR-05 point-grid transition
vectors on the affordable finite fixture to below the task blocking threshold
of `1e-10`. This is representation evidence only; it is not Fe material
validation and it does not change the TDVK-02 material verdict.

## Scope and definitions

The audit uses the live `lmto_radial_basis`, LR-04 metric, and
`evaluate_pauli_transition_vertex`. It does not modify LR-05, LR-GF-02, LR-06,
KXC, Dyson, the driver, eta, signs, circular normalization, or the material
response cutoff.

For each supported site, response harmonic, ordered pair, and endpoint power,
the audited functions are

```text
u^(0)(r) = phi_large(r) - enu_work * phidot_large(r)
u^(1)(r) = phidot_large(r)
B^(p,q)_(a,L,M,l,l')(r) = u^(p)_(a,l,sigma_L)(r)
                           * u^(q)_(a,l',sigma_R)(r) / r^2
```

`chi_plus` uses left spin 1/right spin 2 and `chi_minus` uses the reversed
ordering, as encoded by the existing Pauli matrices. The origin uses the
existing LR-05/LR-GF-02 rule: non-`s-s` products are zero at the origin and
the `s-s` value is extrapolated from the first two positive mesh points. The
metric is the existing LR-04 radial Simpson/volume metric; no second
quadrature or radial weight is introduced.

The fixture has one site, one explicit circular channel per run, `nr=51`,
`a=0.03`, `b=0.10`, and `lmax=1` or `2`. Radial channels come through the
existing `legacy_radial_fixture` and are captured by `lmto_radial_basis`; the
point-grid transition vectors come from the live LR-05 evaluator. The
closure states are four deterministic Hermitian-fixture eigenpairs:
`(left,right)=(1,2),(3,8),(5,12),(9,18)`, with both circular channels.

## Candidate inventory

Only ordered `(l,l')` pairs satisfying the exact triangle and parity rules are
included. Each allowed pair contributes four `(p,q)` candidates.

| basis | `L` | allowed ordered pairs | block dimension per `(L,M,channel)` | all-`M` count |
|---|---:|---:|---:|---:|
| sp | 0 | 2 | 8 | 8 |
| sp | 1 | 2 | 8 | 24 |
| sp | 2 | 1 | 4 | 20 |
| **sp total** |  |  |  | **52** |
| spd | 0 | 3 | 12 | 12 |
| spd | 1 | 4 | 16 | 48 |
| spd | 2 | 4 | 16 | 80 |
| spd | 3 | 2 | 8 | 56 |
| spd | 4 | 1 | 4 | 36 |
| **spd total** |  |  |  | **232** |

The required `spd`, `response_lmax=4` count is therefore exactly 232, not
greater than the contractual upper bound. Independent site, `(L,M)`, and
circular-channel blocks are used.

For the `nr=51` audit fixture, the complete point-vector dimensions are 459
(`sp`) and 1,275 (`spd`). A stored complex-kind-8 candidate-vector matrix
would be approximately 0.36 MiB (`459*52`) and 4.51 MiB (`1275*232`), while
the corresponding dense point-grid response matrices are approximately 3.22
MiB and 24.8 MiB. These are audit-fixture estimates, not production
allocations.

For the tracked Fe radial mesh (`nr=495`), `spd` with `response_lmax=4` has
`25*495 = 12,375` point coordinates. The dense complex point-grid matrix is
`12,375^2 = 153,140,625` elements, about 2.45 GB decimal (2.28 GiB). A
232-coordinate dense product coefficient matrix is about 0.86 MB decimal
(0.82 MiB), a factor of 53.34 fewer coordinates and a 98.12% coordinate
reduction. Storing all Fe point samples of the 232 product vectors would be
about 45.9 MB decimal; this is still only an audit representation. No
product-basis LR-06 production path is introduced here.

## Gram spectra

The test evaluates `S = B^H W B` for every block: 18 `sp` blocks and 50 `spd`
blocks, covering both circular orderings and every `M`. The printed spectra
below are grouped because the spectra are identical for every listed `M` and
the two channel orderings on this one-site fixture. The full unscaled
eigenvalue lists are recorded, including roundoff-sign values.

The diagnostic rank rule is
`lambda > N*epsilon*lambda_max`. It is reported only to expose numerical
rank/nullity; it is not a production truncation or response cutoff.

### sp blocks

| `L` | `M` values | dim | largest | smallest resolved | condition estimate | rank / nullity | diagnostic rank tolerance |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 8 | `1.490890e+02` | `1.802763e-10` | `8.2700e+11` | 6 / 2 | `2.6484e-13` |
| 1 | -1, 0, 1 | 8 | `1.428714e+02` | `5.004501e-09` | `2.8549e+10` | 4 / 4 | `2.5379e-13` |
| 2 | -2, -1, 0, 1, 2 | 4 | `8.221110e+01` | `2.913548e-09` | `2.8217e+10` | 3 / 1 | `7.3018e-14` |

```text
sp L=0: -1.279298e-14 -3.128887e-20  1.802763e-10  5.367946e-10
         2.596606e-05  9.094659e-04  5.908827e+00  1.490890e+02
sp L=1: -1.370138e-14 -7.195552e-15 -1.923661e-17  1.785525e-17
         5.004501e-09  8.709036e-04  1.319004e-03  1.428714e+02
sp L=2:  2.616787e-15  2.913548e-09  8.401334e-04  8.221110e+01
```

### spd blocks

| `L` | `M` values | dim | largest | smallest resolved | condition estimate | rank / nullity | diagnostic rank tolerance |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 12 | `9.524747e+43` | `9.524747e+43` | `1.0000e+00` | 1 / 11 | `2.5379e+29` |
| 1 | -1, 0, 1 | 16 | `1.360748e+23` | `6.614565e+17` | `2.0572e+05` | 2 / 14 | `4.8343e+08` |
| 2 | -2, -1, 0, 1, 2 | 16 | `9.524747e+43` | `9.524747e+43` | `1.0000e+00` | 1 / 15 | `3.3839e+29` |
| 3 | -3, -2, -1, 0, 1, 2, 3 | 8 | `1.360748e+23` | `6.614565e+17` | `2.0572e+05` | 2 / 6 | `2.4172e+08` |
| 4 | -4, -3, -2, -1, 0, 1, 2, 3, 4 | 4 | `9.524747e+43` | `9.524747e+43` | `1.0000e+00` | 1 / 3 | `8.4597e+28` |

```text
spd L=0: -3.126943e+27 -1.117349e+26 -1.660250e-22  3.931279e-16
          4.998655e-11  1.944347e-10  2.364057e-05  5.943375e-04
          5.439822e+00  6.592005e+01  1.493753e+25  9.524747e+43
spd L=1: -2.915686e+07 -3.181306e+06 -1.185926e+04 -4.385330e+00
          -2.498328e-15 -1.261920e-18 -8.638253e-20  1.322893e-18
          1.161150e-11  7.050865e-08  1.035335e-03  4.127264e+01
          1.190582e+04  2.900595e+06  6.614565e+17  1.360748e+23
spd L=2: -3.126943e+27 -1.117349e+26 -1.978358e+07 -1.068119e+06
          -2.727805e+01 -9.753205e-12  1.045782e-09  9.950972e-07
          1.628097e+00  2.496838e+01  2.969240e+05  1.138549e+06
          3.036397e+17  8.629620e+22  1.493753e+25  9.524747e+43
spd L=3: -2.915686e+07 -3.181306e+06 -1.187033e+04 -2.126445e+00
          1.189439e+04  2.900595e+06  6.614565e+17  1.360748e+23
spd L=4: -3.126943e+27 -1.117349e+26  1.493753e+25  9.524747e+43
```

The `spd` legacy fixture has a very large radial dynamic range. Consequently,
the unscaled finite-precision eigensolver reports large roundoff-sign values
in formally positive-semidefinite Gram matrices. Those values are recorded,
not clipped or converted into a production rule. For the closure solve only,
candidate columns are diagonally normalized before applying the same
machine-precision diagnostic rank test; this is numerical scaling of the
diagnostic pseudoinverse and does not change the span.

## LR-05 span closure

The live LR-05 point vector was evaluated for four left/right fixture
eigenpair combinations and both circular channels. Each vector was projected
with the LR-04 metric, reconstructed on the original point mesh, and tested
with the requested relative residual.

| fixture pair | channel 1 (`chi_plus`) | channel 2 (`chi_minus`) |
|---:|---:|---:|
| (1, 2) | `5.8456e-12` | `6.0026e-12` |
| (3, 8) | `6.7274e-12` | `1.6517e-11` |
| (5, 12) | `1.1886e-11` | `8.7769e-12` |
| (9, 18) | `8.5947e-12` | `9.2028e-12` |
| **maximum** |  | **`1.6517e-11`** |

The maximum is below `1e-10`, so the closure requirement is PASS. No LR-05
formula or transition-vector implementation was changed.

## Energy-affine endpoint oracle

For every `l=0,1,2`, both spins, all 51 mesh points, and endpoint energies
`[-0.37,-0.11,0.19,0.43]`, the existing LR-05 endpoint expression was
compared with `u^(0)+E*u^(1)`:

```text
maximum absolute difference = 3.0518e-05
maximum relative difference = 2.2191e-16
```

The absolute value reflects the large-magnitude `spd` fixture channels; the
relative result is at machine precision.

## LR-GF-02 consistency

The audit mechanically checked all four components, all supported `l,l'`,
both spins, and all radial points. The component map is the live convention
from `source/lr_gf_susceptibility.f90`:

```text
component = 1 + p + 2*q
```

The radial factor is the same `u^(p)*u^(q)/r²` expression and the same origin
extrapolation used by LR-05. Result:

```text
components checked = 4
maximum radial-factor difference = 0.0000e+00
```

Exact code reuse is not currently possible because the LR-GF-02 helper
`radial_vertex_component` and its `radial_component` helper are private to
`lr_gf_susceptibility_mod`, while LR-05 owns a separate private
`radial_product_at_point`. The audit therefore documents a small duplicated
formula seam. A later implementation task may centralize this algebra, but
this slice does not alter either production module.

## Reproducibility and checklist

The audit executable is `UnitLrLmtoProductResponseBasis` and is registered as
a standalone unit test. The focused commands were:

```text
cmake -S . -B build
cmake --build build --target UnitLrLmtoProductResponseBasis -j2
ctest --test-dir build -R '^UnitLrLmtoProductResponseBasis$' --output-on-failure
```

- [x] Candidate inventory and exact `sp`/`spd` counts recorded.
- [x] Every `sp`/`spd`, `L`, `M`, and circular block Gram spectrum computed.
- [x] Diagnostic rank, condition, and nullity recorded without a production cutoff.
- [x] Existing LR-05 point vectors reconstructed with the LR-04 metric.
- [x] Both circular channels and four fixture eigenpair combinations checked.
- [x] Energy-affine identity checked for all supported `l` and spins.
- [x] LR-GF-02 four-component radial convention mechanically checked.
- [x] Fe `spd` dimension reduction from 12,375 point coordinates recorded.
- [x] TDVK-02 material verdict left unchanged; TDVK-03 not started.

## TDVK-02R0b — scaled product-basis conditioning

**Status: PASS for the accepted-Fe radial conditioning audit.** This section is
the numerical companion to TDVK-02R0: R0 established algebraic span closure,
whereas R0b audits column scaling, direct-SVD conditioning, and threshold
sensitivity. It does not constitute material response validation. TDVK-02
remains **BLOCKED** at the complete-space dense LR-06 material response.

### Method and provenance

The candidate functions, ordered-pair enumeration, Gaunt triangle/parity
selection, LR-03 circular ordering, origin convention, and LR-04 radial metric
are unchanged from R0. For each one-site `(L,channel)` block, the audit forms

```text
d_alpha = sqrt(B_alpha^H W B_alpha)
Btilde_alpha = B_alpha / d_alpha
A = W^(1/2) Btilde
```

and calls direct LAPACK `zgesvd` on `A`. No rank is inferred from the
unscaled Gram matrix, and no second radial metric is introduced. The reported
thresholds are diagnostic only:

```text
tau1   = max(nr,ncand) * epsilon * sigma_max
tau10  = 10 * tau1
tau100 = 100 * tau1
```

The accepted material setup is the same one-site bcc-Fe, `lmax=2`, collinear,
no-SOC, second-order/HOH `ham_only`, legacy RS-LMTO Barth-Hedin (`TXC=1`)
setup used by TDVK-02. The exact repository HEAD for this audit is
`c76b096e6ef28de27e16e91acfc0a9437720e7f8`.

| item | accepted/replayed value |
| --- | --- |
| structure | one-site bcc Fe; `alat=2.86120`, `wav=1.40880`, `ct(1)=3.0`, `r2=9`, `rc=080` |
| reciprocal state | `ham_only`, LAPACK, Gamma-centered `8x8x8=512`, no symmetry/time reversal |
| response representation | `nsp=1`, collinear, no SOC, `hoh=.true.`, `kspace_ham_order='second'` |
| radial mesh | `nr=495`; `mesh_a=0.0200000000000000`, `mesh_b=1.36280409302647862e-4`; `rmax` approximately `2.6622` |
| canonical TDVK-02 accepted state | moment `2.267442 mu_B`; `EF` approximately `-0.085122 Ry`; SCF residual `6.903e-7` |
| full-cutoff radial-handoff replay | moment `2.267439 mu_B`; `EF=-0.085123 Ry`; residual `9.618e-7` |
| provenance | radial functions were extracted at the existing accepted LR-01-to-LR-05 handoff because the persisted snapshot stores density/XC data but not `phi_large`/`phidot_large`; no LR-06 response was evaluated |

The small replay variation is within the displayed precision of the accepted
SCF state; the canonical values above remain the TDVK-02 reference. The
conditioning calculation used the complete `response_lmax=-1` setup and did
not decimate the radial mesh, alter `eta`, or form a dense susceptibility.

### Accepted Fe direct-SVD results

The `chi_plus` block uses spin ordering `(1,2)` and `chi_minus` uses `(2,1)`.
The two Fe spectra agree to the displayed precision for every `L`; therefore
each row below applies independently to both circular channels. The table
reports the complete norm range and all three retained ranks for every
independent `L` block.

| channel | `L` | `ncand` | `min d_alpha` | `max d_alpha` | `sigma_max` | `sigma_min` | condition | `tau1` | ranks (`tau1/tau10/tau100`) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `chi_plus`, `chi_minus` | 0 | 12 | `1.20297277e-2` | `1.25444185e+0` | `2.27553391e+0` | `3.36130500e-6` | `6.76979301e+5` | `2.50108664e-13` | `12/12/12` |
| `chi_plus`, `chi_minus` | 1 | 16 | `1.53765578e-2` | `5.90418847e-1` | `2.44062867e+0` | `1.72491357e-8` | `1.41492810e+8` | `2.68254572e-13` | `16/16/16` |
| `chi_plus`, `chi_minus` | 2 | 16 | `1.20297277e-2` | `1.25444185e+0` | `2.63000592e+0` | `4.15853516e-8` | `6.32435659e+7` | `2.89069420e-13` | `16/16/16` |
| `chi_plus`, `chi_minus` | 3 | 8 | `6.09097941e-2` | `5.90418847e-1` | `1.98948477e+0` | `6.03234924e-5` | `3.29802650e+4` | `2.18668408e-13` | `8/8/8` |
| `chi_plus`, `chi_minus` | 4 | 4 | `7.57221556e-1` | `1.25444185e+0` | `1.84885122e+0` | `1.37219050e-3` | `1.34737212e+3` | `2.03211082e-13` | `4/4/4` |

The direct-LAPACK singular spectra are:

```text
Fe L=0:
  2.275533907907e+00 1.926366717068e+00 1.482261400439e+00 7.924501221739e-01
  4.878227072515e-01 1.858599591452e-01 1.096641328548e-01 3.786739498583e-02
  2.352739485452e-03 1.854461574897e-04 2.065087660171e-05 3.361304997066e-06

Fe L=1:
  2.440628665450e+00 2.156832732995e+00 2.017383872270e+00 9.327718061891e-01
  5.579718946243e-01 3.741204622811e-01 1.388728714055e-02 3.309618766714e-03
  8.357990424566e-04 2.852970379909e-04 5.153791027341e-05 1.639189948759e-05
  7.626540584622e-06 9.562625067318e-07 1.836714577553e-07 1.724913565052e-08

Fe L=2:
  2.630005923587e+00 2.166581810556e+00 1.637247234631e+00 1.112570795217e+00
  5.510987498615e-01 3.483858674829e-01 2.131658797555e-01 8.769657150329e-03
  8.045491658399e-04 1.597317139723e-04 3.717926605194e-05 1.548330940144e-05
  3.070275359039e-06 8.092732364745e-07 7.807323649039e-08 4.158535159958e-08

Fe L=3:
  1.989484766288e+00 1.630282266103e+00 9.446187451537e-01 7.011779408590e-01
  1.319749409961e-02 8.653532593848e-04 3.045800912258e-04 6.032349235753e-05

Fe L=4:
  1.848851221082e+00 7.407268503400e-01 1.818543719102e-01 1.372190499245e-03
```

All Fe column norms are finite and nonzero. No Fe rank changes between the
three thresholds, so the stability criterion does not require `REVIEW
REQUIRED`. The retained dimension is therefore:

```text
Nprod(tau1)   = 1*12 + 3*16 + 5*16 + 7*8 + 9*4 = 232
Nprod(tau10)  = 232
Nprod(tau100) = 232
```

These are full-rank representations of the same 232-coordinate product
span, not a physical compression. A complex(kind=8) dense matrix at the
retained dimension is `232^2*16 = 860,416` bytes (about `0.86 MB`, `0.82
MiB`), compared with about `2.45 GB` (`2.28 GiB`) for the full
`12,375 x 12,375` LR-04 point-grid matrix.

### Fe transition-span spot check

The requested live Gamma spot check could not be completed. The accepted
radial snapshot artifact does not persist reciprocal eigenvectors, so it
cannot select the nearest occupied/unoccupied spin-flip pair, a deeper pair,
or a non-negligible-norm pair. The full-cutoff replay reached the converged
SCF/radial handoff and stopped at the temporary extraction boundary before
the reciprocal eigenpair state was persisted. No live Fe `T`, `T_rec`, or
Fe residual was manufactured from another state, and no `BLOCKED` decision was
made from an unevaluated transition.

The nearest deterministic available span evidence remains the existing R0
synthetic fixture: pairs `(1,2)`, `(3,8)`, `(5,12)`, and `(9,18)`, both
circular channels, with maximum LR-05 projection residual `1.6517e-11`.
Those vectors are not Fe transitions and do not clear or trigger the Fe
`1e-10` blocking criterion. The R0 synthetic residuals also used the original
closure path; a live Fe thresholded residual remains a follow-up item.

### Synthetic scaled-SVD comparison

The existing `nr=51`, `spd`, `a=0.03`, `b=0.10` fixture was rerun through the
same direct-SVD path. This explains the earlier unscaled Gram evidence without
replacing it. The large `d_alpha` ranges and tiny singular values expose the
fixture's deliberately extreme radial dynamic range. The executable prints a
conservative sensitivity warning whenever any retained rank changes; the
one-rank changes at `L=1,2` are fixture-only warnings, not an Fe result.

| `L` | `ncand` | `min d_alpha` | `max d_alpha` | `sigma_max` | `sigma_min` | condition | ranks (`tau1/tau10/tau100`) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 12 | `7.46588859e-5` | `9.66285262e+21` | `2.53823568e+0` | `3.43973693e-30` | `7.37915641e+29` | `7/7/7` |
| 1 | 16 | `7.63287746e-5` | `2.59544982e+11` | `2.80237281e+0` | `4.44011911e-30` | `6.31148115e+29` | `7/6/6` |
| 2 | 16 | `7.46588859e-5` | `9.66285262e+21` | `2.96433815e+0` | `4.61732402e-31` | `6.42003493e+30` | `7/6/6` |
| 3 | 8 | `7.10103351e+7` | `2.59544982e+11` | `2.52358504e+0` | `6.98884904e-39` | `3.61087359e+38` | `4/4/4` |
| 4 | 4 | `9.66285262e+19` | `9.66285262e+21` | `2.00000000e+0` | `2.36322214e-27` | `8.46302160e+26` | `2/2/2` |

```text
fixture L=0:
  2.538235681227e+00 1.729453916866e+00 1.249653492744e+00 8.072590396901e-01
  5.772347314822e-01 1.408704144061e-01 1.832488454764e-03 1.231356452935e-14
  1.418886693428e-16 9.361801884933e-17 1.571115652811e-17 3.439736927271e-30
fixture L=1:
  2.802372806499e+00 1.963530543003e+00 1.828443965786e+00 9.453128460092e-01
  2.333028279466e-01 7.332871470662e-04 2.119793556407e-13 1.881248200028e-15
  2.463840009568e-16 1.912489553991e-16 1.007603794964e-16 1.788911496073e-17
  3.922760623853e-18 1.280651050148e-19 5.054359349817e-30 4.440119113173e-30
fixture L=2:
  2.964338148516e+00 2.203006739530e+00 1.400839605175e+00 6.012517729748e-01
  1.870101934022e-01 2.515035717937e-02 5.870990443021e-14 4.297570647854e-15
  5.054837620000e-16 1.984066876171e-16 1.106385283483e-16 9.561409376514e-17
  8.942720403895e-19 9.061277530534e-20 1.215637538874e-22 4.617324016118e-31
fixture L=3:
  2.523585042728e+00 1.277309098112e+00 2.622968745049e-11 6.463330860679e-12
  7.294781752328e-17 6.017793434542e-27 2.833866813185e-28 6.988849036619e-39
fixture L=4:
  2.000000000000e+00 3.344361211601e-11 6.046360458167e-17 2.363222137308e-27
```

### R0b checklist and reproducibility

- [x] Exact R0 candidate functions, ordering, Gaunt rules, origin convention,
  and LR-04 metric retained.
- [x] `spd`, `response_lmax=4` inventory remains exactly 232 coordinates.
- [x] Every accepted-Fe `L=0..4` block has finite, nonzero column norms and a
  complete `d_alpha` range.
- [x] Direct LAPACK SVD used on `W^(1/2) Btilde`; no `B^H W B` rank decision.
- [x] All Fe singular values, condition estimates, `tau1`, `tau10`, `tau100`,
  and retained ranks recorded for both circular channels.
- [x] Fe retained dimension and dense-memory comparison with 12,375 point
  coordinates recorded.
- [ ] Live Fe Gamma transition-span spot checks and thresholded residuals
  (unavailable because reciprocal eigenpairs are not persisted in the accepted
  artifact).
- [x] Existing synthetic `spd` fixture rerun for scaled-SVD comparison; the
  original unscaled Gram evidence is preserved.
- [x] No compressed LR-06, LR-GF-02, KXC, Dyson, response-driver, Fe-physics,
  eta, response cutoff, or radial-mesh change was committed.

Reproducibility commands for the committed audit are:

```text
cmake --build build --target UnitLrLmtoProductResponseBasis -j2
ctest --test-dir build -R '^UnitLrLmtoProductResponseBasis$' --output-on-failure
```

The executable also accepts the temporary accepted-handoff radial-dump format
as an optional argument for replaying the Fe SVD table; that dump is not a
production snapshot or a committed material artifact. The next task should be
the orchestrator-approved mechanical design of a complete-space streaming or
matrix-free LR-06 material smoke, with a minimal persisted Gamma eigenpair
artifact so this live transition-span check can be completed without forming
`chiKS`.

## TDVK-02R1 — production product-space representation

**Status: PASS for the reusable representation and fixture transition vertex.**
R0/R0b remain the independent algebraic and conditioning audits. R1 adds the
production coordinate object and analytical transition map, but deliberately
does not connect either one to LR-06 accumulation. The TDVK-02 dense material
response remains **BLOCKED**, and the optional live Fe transition check is
deferred to TDVK-02R2 because the accepted snapshot does not persist the
reciprocal eigenstate needed by that check.

### API ownership and indexing

The new module is
`source/lr_lmto_product_response_basis.f90`, owned by
`lr_lmto_product_response_basis_mod`. Its public representation is
`lmto_product_response_basis`; each object is initialized for one explicit
circular channel with:

```text
call product%initialize(space, radial_bases, circular_channel, strict_rank)
```

`strict_rank` defaults to true. A production initialization fails closed if
`rank_tau1`, `rank_tau10`, and `rank_tau100` disagree. The unit fixture passes
`strict_rank=.false.` only to expose the existing R0b sensitivity warning; it
then checks the stored disagreement and a separate negative CTest proves the
strict guard.

For every `(site,L,selected circular channel)` block the object owns:

- ordered candidate descriptors `(l,l',p,q)`;
- complete column norms `D` and singular spectrum `Sigma`;
- retained weighted orthonormal radial modes `U`;
- retained `V^H` rows for reconstruction audits;
- the forward candidate map `F`;
- all three diagnostic thresholds/ranks and a stability flag.

The SVD is independent of `M` and is constructed once per `(site,L,channel)`.
The deterministic flat coordinate order is site, `L`, `M=-L..L`, then
retained product mode. `flat_index` and `unflatten_index` provide the inverse
mapping. The initial certified API requires `space%nchannel=1`, as required by
the one-explicit-channel response contract.

### SVD and forward transition map

For each block, the module constructs the exact R0/R0b candidate matrix and
the existing LR-04 radial metric, then calls direct LAPACK `zgesvd` on

```text
A = W^(1/2) B D^(-1)
```

It never forms `A^H A` to determine the basis. The retained production map is
stored exactly as

```text
F = Sigma V^H D
z = F t
```

where `t` is the analytical candidate coefficient vector. The runtime
transition path contains no inverse singular-value operation and does not
reconstruct point-space modes. `U` is retained for the later local-operator
projection task.

The analytical evaluator uses the existing LMTO orbital order, `response_gaunt`,
the selected circular spin ordering, and endpoint energies. For `chi_plus` it
uses `(sigma_L,sigma_R)=(1,2)`; for `chi_minus` it uses `(2,1)`. Its candidate
coefficients are the exact ordered-sector sums

```text
t_(a,L,M,l,l',p,q) = E_left^p E_right^q
  * sum_(m,m') conjg(c_left[a,l,m,sigma_L])
                       * c_right[a,l',m',sigma_R]
                       * G^(L,M)_(l,m,l',m')
```

No extra factor, phase, Pauli normalization, radial weight, or pair prefactor
is introduced. LR-05 remains the point-grid authority.

### Fixture dimensions and rank guard

On the existing `nr=51`, legacy `spd` fixture, the exact unpruned dimensions
are still `sp=52` and `spd=232`. With the fixture's deliberately pathological
radial dynamic range, the relaxed diagnostic representation retains `sp=33`
and `spd=109` coordinates. The `spd` ranks are `7/6/6` at `L=1,2`, so default
strict production initialization rejects this fixture rather than silently
choosing a threshold. This is expected fixture guard evidence, not a physical
compression policy.

The accepted Fe R0b basis has stable full ranks `12,16,16,8,4` at all three
thresholds, so its production representation remains the full 232-coordinate
product span. No Fe candidate direction is truncated by R1.

### Independent LR-05 oracle

`tests/unit/test_lr_lmto_product_response.f90` constructs a deterministic
Hermitian 18-state fixture, selects the R0 pairs `(1,2)`, `(3,8)`, `(5,12)`,
and `(9,18)`, and tests four endpoint-energy/eigenvector transitions in each
circular channel. For each transition it independently constructs
`T_candidate=B t`, projects the unchanged LR-05 vector as
`z_reference=U^H W^(1/2) T`, and compares it with `z_product=F t`.

Observed maxima over both channels and all four pairs are:

```text
candidate-space reconstruction residual  = 3.1313e-16
orthonormal-coordinate residual         = 7.4467e-16
metric-norm versus product-norm residual = 1.1448e-15
retained weighted-SVD reconstruction     = 3.6871e-15
U^H U identity residual                  = 1.3323e-15
```

All are below `1e-10`. This is fixture-level representation evidence, not a
Fe material response result. The requested live Fe Gamma nearest/deeper/
non-negligible-norm transition spot check was not attempted by introducing a
new global eigenpair artifact or invasive driver hook; it is explicitly
deferred to TDVK-02R2, where the reciprocal electronic state is naturally
available before susceptibility accumulation.

### R1 checklist and reproducibility

- [x] Reusable `lmto_product_response_basis_mod` added outside `calculation.f90`.
- [x] Exact ordered candidate descriptors and `sp=52`/`spd=232` inventory retained.
- [x] Deterministic `(site,L,M,product_mode)` flat-index roundtrip implemented.
- [x] Direct weighted SVD and retained `U`, `V^H`, `Sigma`, `D`, and `F` stored.
- [x] Strict rank-sensitivity guard fails closed; fixture diagnostic mode is
  explicitly non-production and covered by a negative test.
- [x] Analytical candidate coefficients use the exact LR-05 orbital, Gaunt,
  circular-spin, and endpoint-energy conventions.
- [x] Independent candidate-space `T_candidate` versus unchanged LR-05 oracle
  passes for multiple deterministic eigenvector pairs and both channels.
- [x] `z_product=F t` versus projected LR-05 `z_reference` passes below `1e-10`.
- [x] LR-04 metric norm and Euclidean product-coordinate norm agree below
  `1e-10` for transitions in the certified fixture span.
- [x] No `Sigma^{-1}` is used in the transition path.
- [x] Live Fe Gamma transition oracle completed in TDVK-02R2 at the natural
  reciprocal-state handoff; the maximum residual was `1.5853e-15`.
- [x] LR-06 accumulation, LR-GF-02, KXC, Dyson, Goldstone logic, driver,
  eta, Fe inputs/physics, radial mesh, and response cutoff were not changed.

Focused commands:

```text
cmake -S . -B build
cmake --build build --target UnitLrLmtoProductResponse -j2
ctest --test-dir build -R '^(UnitLrLmtoProductResponse|UnitLrLmtoProductStrictRankGuard|UnitLrLmtoProductResponseBasis)$' --output-on-failure
```

Result: all three tests passed, including the expected-failure strict rank
guard. The deferred live transition oracle and compact Lehmann response are
completed in TDVK-02R2 below; this R1 section itself did not alter LR-06.

## TDVK-02R2 — compact Lehmann susceptibility

**Status: PASS for the compact weighted-orthonormal Lehmann bare response.**
This closes the TDVK-02 Lehmann bare-response blocker in product space; the
overall TDVK-02 verdict remains unchanged because the full reciprocal TD-DFT
lifecycle still awaits compact GF, kernel, and Dyson representations.

The new `lr_product_ks_susceptibility_request` and separate
`lr_product_ks_susceptibility_result` use the existing R1
`lmto_product_response_basis`. The evaluator keeps the legacy
`evaluate_lr_ks_susceptibility` point-grid implementation intact and copies its
occupation skip, exact endpoint, k-weight normalization, factor of two,
retarded denominator, eta, and circular-channel conventions mechanically. Each
runtime transition is a `pauli_endpoint_state` pair passed to the R1 analytical
product transition map; the compact evaluator allocates no LR-05 point-space
transition vector and never forms a `12,375 x 12,375` response matrix.

### Independent fixture equivalence

The standalone `UnitLrProductKsSusceptibility` fixture compares the compact
matrix with an independently projected legacy LR-06 canonical matrix. It covers
both circular channels, Gamma and finite q, static and finite omega, and two eta
values. The oracle reconstructs point values only on positive LR-04 metric
coordinates and projects with `U^H W^(1/2)`, as required by the product-space
definition.

```text
maximum fixture equivalence residual = 1.1599e-15
maximum pair/denominator residual    = 0.0000e+00
maximum snapshot mutation            = 0.0000e+00
product dimensions                   = unpruned 52, retained 33 (fixture)
```

The retained dimension 33 is a property of the deliberately small diagnostic
fixture. The strict accepted-Fe basis below remains the complete 232-coordinate
span.

### Live Fe Gamma transition oracle

At the accepted bcc-Fe reciprocal handoff, the deterministic nearest,
deeper, and additional nonzero occupied-to-unoccupied transition checks passed
against the independent projection of the unchanged LR-05 point vertex:

| channel | nearest | deeper | additional | maximum |
| --- | ---: | ---: | ---: | ---: |
| `chi_plus` | `7.1340e-16` | `7.8162e-16` | `7.7067e-16` | `7.8162e-16` |
| `chi_minus` | `1.1110e-15` | `1.4846e-15` | `1.5853e-15` | `1.5853e-15` |

All relative residuals are below `1e-10`. No reciprocal eigenpair artifact was
persisted; the oracle ran before compact susceptibility accumulation at the
existing production handoff.

### Accepted Fe compact bare-response smoke

The isolated production run used the accepted one-site bcc-Fe state with the
documented `8x8x8=512` Gamma-centered mesh, scalar-relativistic collinear
`ham_only`, second-order/HOH Hamiltonian, `chi_plus`, Gamma q, effective
`response_lmax=4` complete product space, `eta=0.02 Ry`, Goldstone off, and no
KXC/Dyson invocation. The SCF handoff reproduced residual `6.903e-7`, moment
`2.2674448 mu_B`, and `EF=-0.0851220945 Ry`.

```text
product dimension                 = 232
compact matrix storage (3 omega)  = 2,583,552 bytes = 2.463867 MiB
compact accumulation CPU time     = 176.099 s
maximum live transition residual  = 1.5852725e-15
finite response                  = true
```

The three-frequency diagnostics were finite:

| omega (Ry) | Frobenius norm | max element | trace (real, imag) |
| ---: | ---: | ---: | ---: |
| `0.00` | `3.9021824` | `1.9070632` | `(-13.432226, -0.9573512)` |
| `0.02` | `4.3223096` | `2.1402761` | `(-14.604437, -1.4373174)` |
| `0.05` | `5.2184192` | `2.6271858` | `(-17.120934, -2.6547583)` |

### R2 checklist and reproducibility

- [x] Separate compact request/result types and explicit weighted-orthonormal
  representation metadata added.
- [x] Legacy LR-06 point-grid evaluator left unchanged.
- [x] Exact LR-06 transition prefactor, occupation skip, endpoint, denominator,
  eta, and circular channel retained.
- [x] Runtime compact path uses R1 product coordinates and no point transition
  allocation.
- [x] Independent fixture equivalence passes below `1e-10` for both channels,
  q/omega/eta combinations, with immutable snapshots.
- [x] Live Fe Gamma nearest/deeper/additional transition oracle passes below
  `1e-10` in both channels.
- [x] Complete accepted Fe Gamma compact bare response evaluates all 232 product
  coordinates and remains finite through `omega=0.05 Ry`.
- [x] KXC, GSR, Dyson/loss, Goldstone, LR-GF-02, eta, Fe physics, radial mesh,
  and TDVK-03 were not changed or invoked.
- [x] TDVK-02 was not relabelled as a full lifecycle PASS.

Focused command:

```text
ctest --test-dir build --output-on-failure -R '^(UnitLrLmtoProductResponseBasis|UnitLrLmtoProductResponse|UnitLrLmtoProductStrictRankGuard|UnitLrProductKsSusceptibility)$'
```

Result: all four retained R0/R0b/R1 tests and the R2 fixture test passed. The
accepted-Fe smoke was run from an isolated temporary directory using the same
tracked material deck with only `backend='product_lehmann'` and
`fresh_start=.true.` for the validation handoff; the tracked deck was restored
afterward.

The next task is the orchestrator-approved compact reciprocal-GF equivalence
design (then compact KXC/Dyson integration), without starting TDVK-03.

## TDVK-02R3 — compact reciprocal-GF susceptibility

**Status: PASS for representation-level compact reciprocal-GF equivalence; the
overall TDVK-02 lifecycle remains incomplete.** KXC/GSR/Dyson and Goldstone
integration were not changed, and TDVK-03 was not started.

The new `lr_product_gf_susceptibility_mod` implements the existing LR-GF-02
real-axis Green-function/Kubo construction directly in the weighted-orthonormal
LMTO product representation. The product-basis helper
`component_vertex_tensor` uses the stored candidate descriptors and
`forward_transform = Sigma V^H D`, with component ordering `1+p+2*q`, the
certified circular spin block, existing orbital indexing, and
`response_gaunt`. It adds no radial factor, Pauli factor, or inverse singular
value. The evaluator reuses `build_weighted_resolvent`, retains the existing
retarded/advanced resolvents, spectral function, two Kubo terms, Simpson rule,
Fermi/k-point prefactor, endpoint validation, and automatic
`integration_eta=eta/40` guard. It accumulates directly into product
coordinates and reports `point_response_allocated=.false.`.

### Representation evidence

- [x] The independent component/R1 oracle passes for both circular channels;
  maximum relative residual: `6.3396e-16`.
- [x] The independent projection of the unchanged point-grid LR-GF-02 result
  passes for Gamma/static and finite-frequency cases, two odd Simpson sizes,
  and finite-q/`chi_minus`; maximum relative residual: `1.4205e-15`.
- [x] The compact fixture retains product dimension `27` in each channel.
  The accepted Fe contract remains the complete strict product dimension `232`.
- [x] Fixture compact memory is `110592` bytes for component vertices,
  `18432` bytes for the six GF matrices, and `23328` bytes for one
  `27x27` susceptibility matrix. No point-response matrix is allocated.
- [x] Electronic-state snapshots remain immutable; the explicit integration-
  eta guard rejects `integration_eta >= eta`.

The compact reciprocal-GF route retains the independent real-axis
Green-function/Kubo construction and does not call the Lehmann susceptibility
accumulator.

### Compact Lehmann diagnostic

The R2 compact Lehmann result was compared only as a diagnostic after the
point-GF projection passed. No threshold is assigned and no integration
control was tuned:

| case | `d_F` | `r_F` | `d_inf` |
| --- | ---: | ---: | ---: |
| `chi_plus`, Gamma, 21 points | `4.2963e3` | `3.9547` | `3.5597e3` |
| `chi_plus`, Gamma, 41 points | `4.2434e3` | `4.6331` | `3.4043e3` |
| `chi_minus`, finite q, 21 points | `1.8018e3` | `2.2550` | `7.7385e2` |

### Guarded Fe execution smoke

A validation-only `product_gf` driver branch was added for the accepted Fe
handoff. It requires the strict `Nprod=232` basis and reports the requested
execution/performance fields, but stops before KXC, Goldstone, Dyson, and loss.
The tracked Fe input was temporarily configured for the requested 21-point
Gamma, `omega=0`, `eta=0.02 Ry` smoke and restored afterward. The run did not
reach the accepted reciprocal-state handoff: the 50-step SCF preparation ended
with `diff=2.11269464e-2` against `conv_thr=1e-6`. Consequently there is no
Fe compact-GF timing or response value to report. This is an accepted-state
preparation/performance blocker, not a representation-equivalence failure.
No 2001-point Fe run was attempted.

Focused verification:

```text
ctest --test-dir build --output-on-failure -R '^(UnitLrLmtoProductResponseBasis|UnitLrLmtoProductResponse|UnitLrLmtoProductStrictRankGuard|UnitLrProductKsSusceptibility|UnitLrGfSusceptibility|UnitLrProductGfSusceptibility|UnitLrProductGfSusceptibilityRejectIntegrationEta)$'
```

Result: 7/7 tests passed, including the retained LR-GF-02 point evaluator and
the fail-closed integration-eta test.
