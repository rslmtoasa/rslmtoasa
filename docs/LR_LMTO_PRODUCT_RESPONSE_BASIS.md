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
