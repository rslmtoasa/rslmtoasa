# GBT WP7 / Gate G7 report — common rotating-frame density and SCF

Date: 2026-08-06. Branch `fable_v2_gbt_v2`, on top of `b0e8e01` (final G6: PASS).

## Outcome

**Gate G7: PASS.** One rotating-frame density contract now exists, both solvers
fill it, radial up/down channels are projected out of it *after* accumulation
against an explicitly stated axis, the two SCF policies are distinct, the
physicality assertions run every iteration, and the q=0 cone is invariant.

The RS/k-space SCF comparison agrees to `1.3e-4 Ry` in total energy,
`3.1e-3 mu_B` in moment, `2.8e-4 Ry` in Fermi level and exactly in charge —
in every case **smaller than each route's own convergence spread** in its own
resolution parameter, which is the bar G7 sets ("within independently
converged numerical tolerances"). §4 shows the sweeps that establish that.

One limitation is recorded and not papered over: on the single-sublattice
collinear ferromagnet used for the SCF comparison the two policies must, and
do, give identical results, so their distinctness is demonstrated at the
contract level (unit oracles, §3) rather than end-to-end. §6.

---

## 1. What was wrong before WP7

Blueprint §2.5, verified against the code at `b0e8e01`:

* **Real space** (`bands%calculate_moments`) built `dspd(l+soff, ie, na)` by
  projecting the on-site Green function onto `potential%mom` **inside the
  energy loop**, then Simpson-integrated the already-projected channel. The
  axis was whatever `mom` happened to hold at that moment in the SCF cycle and
  appeared nowhere in any interface.
* **k space** did something different: `project_dos_orbitals_*` projected each
  eigenstate onto the same implicit `potential%mom` (through
  `get_site_spin_axis`) while building a *broadened DOS*,
  `calculate_band_moments` integrated that DOS on the output energy grid to get
  `m0/m1/m2`, and `self%compute_kspace_spin_moments_spinor` rebuilt the
  Cartesian moment separately afterwards from the eigenvectors.

Three objects, one implicit axis, two integration schemes. The two routes could
not be compared, and nothing prevented an inconsistent pair of spin channels
from reaching the radial SCF.

## 2. The contract

`source/spin_density.f90` — `spin_density_mod`, class `spin_density`.

Per site `a` and angular channel `l` it retains the **complete** rotating-frame
2x2 spin matrix for the three energy-moment orders the radial SCF consumes:

```
rho(s, s', l, site, order),   order = 1,2,3  ->  int rho dE, int E rho dE, int E^2 rho dE
```

equivalently `(n, m_x, m_y, m_z)` through `sd_cartesian_from_matrix` /
`sd_matrix_from_cartesian`.

Fixed conventions (module header, and asserted by the unit oracles):

| quantity | definition |
| --- | --- |
| matrix element | `rho_{s s'} = psi_s conjg(psi_s')` |
| from a retarded G | `rho = (i/2pi)(G - G^dagger)` |
| Cartesian | `n = Tr rho`, `m_i = Tr(sigma_i rho)`; `rho_ud = (m_x - i m_y)/2` |
| lab frame | `rho_lab = D(phi) rho D^dagger(phi)`, `D = diag(e^{-i phi/2}, e^{+i phi/2})`, `phi = q.R_na` |

These are the *same* signs the pre-WP7 code used on both sides
(`-Im(G_ud + G_du)/pi` for `m_x`, `2 Im(conjg(u) d)` for `m_y`), so the contract
codifies the existing convention rather than silently changing it.

**There is no default projection axis.** `get_axis` is fatal on a site whose
axis was never `set_axis`'d, and no routine in the module reads `potential%mom`.

Key procedures: `accumulate_block` (the only writer), `set_axis`/`get_axis`,
`cartesian`/`site_cartesian`, `project_radial`, `radial_band_moments`
(`q0`, band centre, `q2` — exactly `ql(1,l,isp)`, `gravity_center` before the
`vmad` shift, and `ql(3,l,isp)`), `electron_count`, `check_physicality`,
`max_deviation`, `lab_frame_moment`, `resolve_site_axis`.

## 3. Producers, projection and policies

### 3.1 Two producers, one object

* **Real space** — `bands%accumulate_spin_density_rs` (`source/bands.f90`).
  Accumulates the four Cartesian energy channels from `green%g0` with **no axis
  in them at all**, Simpson-integrates each at orders 0/1/2, and assembles the
  2x2 block. MPI-reduced so every rank holds the same global object.
* **k space** — `reciprocal%accumulate_spin_density_kspace`
  (`source/reciprocal_spin_density.f90`). Accumulates
  `sum_k w_k sum_n f(eps) eps^p psi psi^dagger` directly from eigenvectors,
  k weights and Fermi–Dirac occupations. No energy grid, no broadening.

Both fill the same `spin_density` object (`bands%rf_density`,
`reciprocal%rf_density`).

### 3.2 Projection happens after accumulation, against a stated axis

`bands%calculate_moments` now:

1. builds `rf_density` (axis-free),
2. calls `resolve_density_axes`, which turns `potential%mom` into a *reference
   direction*, runs it through the policy, and `set_axis`es the result,
3. calls `assert_density_physical`,
4. fills `ql`/`gravity_center` from `radial_band_moments`.

`reciprocal%calculate_band_moments` does the same three steps and fills
`band_moments` from `fill_band_moments_from_spin_density`. `projected_dos`
survives as DOS output but no longer defines the SCF density, so its
`potential%mom` projection is no longer an SCF input.

The dead `bands%dspd` array (its only consumer was the loop this replaces) is
removed.

### 3.3 Policies

New `&control` key `density_policy`, validated in `control%check_all`,
default `'constrained_spiral'` (the pre-WP7 behaviour):

| policy | axis | mixed | reported |
| --- | --- | --- | --- |
| `constrained_spiral` | pinned to the imposed reference | charge + longitudinal magnitude `m.ref` | transverse residual `m - (m.ref)ref` and torque `ref x m` |
| `relaxed_reference` | follows the full rotating-frame Cartesian moment; written back to `potential%mom` | charge + full Cartesian moment | nothing left over (residual and torque vanish identically) |

Both routes emit the same per-site `DENSITY_POLICY` log line so the two can be
diffed line for line.

### 3.4 Unit oracles

`tests/unit/test_spin_density.f90` (ctest `UnitSpinDensity`, labels
`unit;gbt;wp7;density`). It calls **no production solver**: the model is a
hand-built three-band spinor system and both producers are re-implemented from
their definitions, so agreement is evidence about the convention, not a
restatement of one implementation.

| test | what it establishes |
| --- | --- |
| Cartesian round trip | matrix <-> `(n,m)` against an explicit sigma algebra |
| producer equivalence | eigenvector sum vs Green-function integral fill the same object; residual is `O(eta)` — quartering `eta` quarters it, so nothing convention-shaped survives the limit |
| spinor convention | k-space producer reproduces `m_x/m_y/m_z` from `conjg(u) d` exactly (1e-13) |
| explicit axis | same density + different *stated* axis -> different radial channels; charge is axis-independent |
| policies distinct | different axis, different `m_long`, nonzero vs identically-zero torque, same charge |
| physicality | good density passes; `\|m\|>n`, non-Hermitian, negative trace and a wrong electron count are each caught |
| policy survives re-accumulation | `zero_density` clears the density and every stated axis but keeps the configured policy (see §6) |
| q=0 cone invariance | a global SU(2) rotation of density **and** stated axis leaves `q0`, band centre, `q2` and both radial channels unchanged to 1e-13; `lab_frame_moment(phase=0)` is the identity; `phase/=0` rotates the transverse moment and preserves `m_z` and `\|m\|` |

## 4. RS vs k-space SCF comparison

Deck: bcc Fe, `alat=2.8612`, `rc=80`, `nsp=3`, `density_policy='constrained_spiral'`,
`nstep=30`, `conv_thr=1e-8`, identical `Fe.nml` start for every row. Only the
solver's own resolution parameter changes. Run outside the tracked example
tree.

| route | resolution | etot (Ry) | sumev (Ry) | charge | moment (mu_B) | E_F (Ry) | SCF converged |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| RS block | `lld=15` | -2541.98121718 | -1.98755316 | 8.000000 | 2.124549 | -0.068600 | no (3e-8) |
| RS block | `lld=21` | -2541.98162079 | -1.99701787 | 8.000000 | 2.103764 | -0.069170 | no (5e-8) |
| RS block | `lld=27` | -2541.98169457 | -1.99873957 | 8.000000 | 2.111318 | -0.069179 | no (4e-8) |
| k space | `8^3` | -2541.97923982 | -1.98654171 | 8.000000 | 2.003451 | -0.066940 | yes |
| k space | `12^3` | -2541.98186288 | -2.00000379 | 8.000000 | 2.087159 | -0.067550 | yes |
| k space | `16^3` | -2541.98120413 | -1.99901935 | 8.000000 | 2.124490 | -0.068040 | yes |
| k space | `20^3` | -2541.98182753 | -2.00313726 | 8.000000 | 2.114450 | -0.069460 | yes |

Cross-route difference at each route's best setting (`lld=27` vs `20^3`)
against the intra-route spread that bounds it:

| quantity | RS spread (lld 15..27) | k spread (8^3..20^3) | cross-route difference |
| --- | ---: | ---: | ---: |
| etot | 4.8e-4 Ry | 2.6e-3 Ry | **1.33e-4 Ry** |
| moment | 2.08e-2 mu_B | 1.21e-1 mu_B | **3.13e-3 mu_B** |
| E_F | 5.8e-4 Ry | 2.5e-3 Ry | **2.81e-4 Ry** |
| charge | 0 | 0 | **0** (8.000000 in every row) |

The cross-route difference is smaller than either route's own convergence
error in every quantity. Charge agrees exactly, which is the sharpest single
statement here: it is the one quantity with no route-dependent truncation, and
both producers hit 8.000000 electrons.

Two honest caveats, neither hidden:

* Neither sequence is monotone. The k-space sequence oscillates on this
  metallic Fermi surface and the RS sequence is bounded by the recursion depth
  limit already accepted at G5 (`lld <= 28` for this ~1000-atom cluster). The
  comparison is therefore an agreement-within-spread statement, not a
  demonstration of a common converged limit — the same standing limitation
  recorded in the WP5/G5 report, unchanged by WP7.
* The RS runs report `Not converged! Diff ~ 4e-8` against `conv_thr=1e-8` at
  30 steps. They are converged to ~5e-8, which is far below every difference
  tabulated above, so it does not affect the verdict.

## 5. q = 0 cone invariance (production)

Same deck, RS route, `nsp=3`, no SOC. Only the starting moment direction
differs: `mom = (0,0,1)` versus `mom = (sin45, 0, cos45)`.

| run | etot (Ry) | charge | moment (mu_B) | final `mom` |
| --- | ---: | ---: | ---: | --- |
| `cone0` | -2541.98162079 | 8.00000000 | 2.10376405 | (0, 0, 1.0000) |
| `cone45` | -2541.98163278 | 8.00000000 | 2.10352145 | (0.7271, 6.2e-8, 0.6866) |

`|delta etot| = 1.20e-5 Ry`, `|delta m| = 2.4e-4 mu_B`, `delta charge = 0` —
a factor 40 below the RS convergence spread of §4. The cone is invariant.

(The `cone45` axis drifts from 45.00 to 46.65 degrees over the SCF. That is
`self%mix_magnetic_moments` acting on a direction that carries no energy
gradient in an isotropic ferromagnet, i.e. pre-existing behaviour orthogonal
to WP7; the energy invariance above is the physical statement.)

## 6. Limitation: policy distinctness is not shown end-to-end

`relaxed_reference` was run through a full SCF on both routes:

| route | policy | etot (Ry) |
| --- | --- | ---: |
| RS block `lld=21` | `constrained_spiral` | -2541.9816207940894 |
| RS block `lld=21` | `relaxed_reference` | -2541.9816208073817 |
| k space `12^3` | `constrained_spiral` | -2541.9818628778057 |
| k space `12^3` | `relaxed_reference` | -2541.9818628778057 |

Both routes reproduce their constrained result under the relaxed policy
(`1.3e-8 Ry` apart on RS, bit-identical on k space). **This is correct, not a
bug**: on a single-sublattice collinear ferromagnet the accumulated moment is
already parallel to the imposed reference, so `m_transverse` and the torque
vanish and the two policies must coincide. The `DENSITY_POLICY
relaxed_reference` log lines confirm the relaxed branch is the one executing.

Distinctness is therefore established at the contract level only (§3.4:
different axes, different longitudinal moments, nonzero vs identically-zero
torque, on a density whose moment is deliberately off-reference). A production
case that separates them needs a genuinely constrained multi-sublattice state
where the density does not relax to the reference — that is WP9 validation-
ladder territory, not WP7, and no such deck was built here.

One wiring bug was found and fixed during this check: the real-space producer
called `restore_to_default`, which reset `policy` back to
`constrained_spiral`, so a selected `relaxed_reference` was silently ignored on
the RS route. Producers now call `zero_density` (clears the density and every
stated axis, keeps the policy), and `UnitSpinDensity` asserts that property so
the same mistake cannot come back.

## 7. Regression evidence

| suite | before | after |
| --- | --- | --- |
| `ctest -L unit` | 21/21 | **22/22** (new `UnitSpinDensity`) |
| `ctest -L regression` | 16/16 | **16/16** |
| `ctest -L scf` | not baselined this session | 23/33 |
| `ctest -L postproc` | not baselined this session | 11/12 |

The regression suite passing unchanged is the load-bearing result: the RS
projection moved from inside the energy loop to after integration, which is
algebraically identical (Simpson is linear and the axis is energy-independent)
and differs only at round-off — and the references confirm it.

### Every example-suite failure, attributed

The `scf`/`postproc` suites were **not** baselined at the start of this task
(only `unit` and `regression` were), so each failure was attributed
individually rather than assumed pre-existing.

| case(s) | classification | how established |
| --- | --- | --- |
| `bccFe_nsp2_block`, `bccFe_nsp2_block_hoh` | pre-existing gfortran-13 DOS tolerance (`totaldos.out` row 1500, `rel 1.9e-6`) | already documented in `GBT_WP0_G0_REPORT.md`; re-confirmed on a stashed pristine tree at `b0e8e01` |
| `bccFe_nsp4_block_spiral_qplus`/`_qminus`, `frozen_magnon_bccFe`, `_auto`, `_auto_scf` | pre-existing stale `strux_backend` decks (5 cases) | already documented in `tests/KNOWN_ISSUES.md` from WP6c; `spiral_qplus` additionally re-confirmed on the stashed pristine tree — it fails on `mom[3] = 1.0` vs a reference of `6.9e-9`, a converged out-of-plane axis, **not** a tolerance miss |
| `surface_fccCu001_block_hoh`, `bulk_Pt2MnGa_block_hoh`, `bulk_Pt2MnGa_block` | **not real** — machine contention during the batched four-suite run | re-run individually on the WP7 tree: **all three pass**, in 32.5 s, 29.2 s and 19.0 s against a 300 s limit (they had reported 349 s, 1064 s and 323 s under load) |
| `orbital_modern_bccFe` | pre-existing genuine timeout | times out at exactly 300.01 s when run alone on the WP7 tree **and** identically on the stashed pristine tree at `b0e8e01` |

So the real, attributable failure set is the 8 pre-existing cases the WP6
integration report already records. WP7 introduces none, and the
`spiral_qplus` axis failure is newly written up in `tests/KNOWN_ISSUES.md`
because it had been folded into the `strux_backend` group without its
distinct symptom being noted.

## 8. Completion checklist

- [x] One full rotating-frame density contract exists — `spin_density_mod`,
      complete 2x2 per site/l for all three energy-moment orders.
- [x] RS and k-space producers populate the same object —
      `bands%accumulate_spin_density_rs`,
      `reciprocal%accumulate_spin_density_kspace`.
- [x] Radial projection uses an explicit current axis — `set_axis`/`get_axis`,
      fatal when unset; no `potential%mom` read inside the contract; the SCF
      projects only through `radial_band_moments`.
- [x] Constrained and relaxed policies are distinct — at the contract level;
      see §6 for what was *not* shown.
- [x] Density physicality assertions pass — Hermiticity, eigenvalues, trace,
      electron count and `|m| <= n`, every iteration, `tol = 1e-6`, fatal on
      violation. No violation was observed in any run reported here.
- [x] q=0 cone invariance passes — 1e-13 in the unit oracle, 1.2e-5 Ry in
      production.
- [x] RS/k-space SCF comparison and G7 PASS/FAIL are reported — §4, **PASS**.

## 9. Files changed

| file | change |
| --- | --- |
| `source/spin_density.f90` | new — the contract |
| `source/reciprocal_spin_density.f90` | new — k-space producer + band-moment projection |
| `tests/unit/test_spin_density.f90` | new — WP7 oracles |
| `source/bands.f90` | RS producer, policy resolution, physicality assertion, projection rewired, dead `dspd` removed |
| `source/reciprocal_projection.f90` | `calculate_band_moments` fills from the contract |
| `source/reciprocal.f90` | `rf_density` member, bindings, submodule interfaces |
| `source/control.f90`, `source/include_codes/namelists/control.f90` | `density_policy` key + validation |
| `source/CMakeLists.txt`, `CMakeLists.txt` | build + test registration |

## 10. Unresolved risks / next task

* Three `scf` cases are close enough to the 300 s ctest limit that they time
  out whenever the suite runs under load (they need 19-33 s idle). That is a
  test-harness fragility, not a physics problem, but it makes batched suite
  runs hard to read. Raising their `TIMEOUT` would fix it; not done here.
* The `Example_bulk_bccFe_nsp4_block_spiral_qplus` axis failure is real,
  pre-existing and untriaged. It touches the finite-q moment direction, so it
  is squarely in GBT territory and should be triaged before WP9 validation
  leans on spiral moments. Estimated size: half a day. **Proposal, not done.**
* Policy separation needs a constrained multi-sublattice deck (§6).
* The G5 recursion-depth bound (`lld <= 28`) still caps how sharply the
  RS/k-space comparison can be stated. Unchanged by WP7.

Next allowed task: **WP8** (little-group symmetry and q lifecycle).
