# LR-METHOD-05R — 2026-10-03

**BLOCKED — ALSDA RESPONSE-REPRESENTATION CLOSURE NOT ESTABLISHED**

The physical ground-state kernel repair is implemented. The raw accepted Fe
Ward residual remains near 24.3% as eta decreases. This discrepancy is unresolved;
there is no material-tuned threshold, interaction fit, Ward correction or Goldstone
repair. The method/refactor tree is not certified closed. No commit was made.

## Identity and changed files

Repository: rslmtoasa/rslmtoasa; branch: `fable_v4`.
Starting SHA and final HEAD: `5a6ba54641a0ae22b65d16a953602ee0912c636f`.
All 05R changes are in the working tree. The material binary was built in Debug
with GNU Fortran 16.2.0 and Apple Accelerate. Its cached build-version string
(`lr-campaign-archive-55-gfd7b-dirty`) is stale; the Git SHA above is the source
base, not a claim of a committed repaired executable.

05R changes:

- `source/radial_ground_state.f90`: document the accepted SR number-density contract.
- `source/linear_response.f90`: optional response-density input, separate result provenance.
- `source/linear_response_kernel_dyson.f90`: accepted SR denominator and validation.
- `source/linear_response_run.f90`: common-kernel consumers, independent Ward action,
  SR/Pauli diagnostics, adapter electron count and projected bare selector.
- `docs/linear_response/CONVENTIONS.md`, `FORMULATION.md`, `PUBLIC_API.md`.
- `tests/unit/linear_response/kernel_dyson/test_alsda_kernel.f90`.
- `tests/unit/linear_response/kernel_dyson/test_compact_interaction.f90`.
- `tests/validation/lr_method_05r.py` and this report.
- `docs/DECISIONS.md`: one blocked-campaign record.

The five pre-existing modified example/regression input decks are outside this
change. `source/exchange.f90`, native RSGF bare formulas, basis transformations,
Pauli reconstruction, Jülich, Mills and rotation constructions are unchanged.
No existing reference files or Jülich/Mills/rotation/RSGF tests were changed.

## Physical kernel and response representation

`m_SR = n_up - n_down`, `bxc_SR = (vxc_up-vxc_down)/2`,
`Kxc_SR = bxc_SR/m_SR` from the same accepted radial functional snapshot.
Small finite denominators are evaluated directly. Active exact zeros fail closed;
the existing null-measure origin extension remains. There are no floors or fits.

Pauli density and its accepted eigensystem/large-component/frozen-core provenance
remain separate response observables. Radial and compact ALSDA use the same
central physical kernel; the compact transform remains `U^H K_point U` in weighted
orthonormal coordinates. Dyson sign/order and the separate vertex factor 2 remain
unchanged.

Old `bxc_SR/m_P` construction on a production path: **NO**.
Goldstone correction = **OFF**.

## A — accepted bcc Fe material evidence

The existing `input_compact_alsda.nml` gate supplies the unchanged lattice, SCF
controls, k mesh and temperature. Only response diagnostics/q/frequency sampling
are selected by the evidence runner. No SCF parameter is tuned.

| Quantity | Value |
|---|---:|
| Lattice | bcc Fe, alat = 2.8612 Å |
| k mesh / count | 4 × 4 × 4 / 64 |
| EF | -0.06224484675548704 Ry |
| Temperature | 300 K |
| Accepted canonical SCF moment (`kspace_scf_state.dat`) | 2.000001309354217 μB |
| Accepted radial integrated M_SR | 2.000007206020199 |
| Integrated M_P | 2.054549316876242 |
| Delta M = M_SR - M_P | -0.0545421108560430 |
| Relative L2(m_SR-m_P), full radial metric | 0.01695015460673785 |
| Magnetic-region relative L2 | 0.01655557555262725 |
| Magnetic-region definition | first radius enclosing 90% of integrated abs(m_SR) |
| Kxc_SR min / max | -4.020555710650229 / 0 Ry bohr³ |
| Compact product dimension / angular cutoff | 348 / 4 |
| Accepted reciprocal electrons / target / difference | 8.000000000002840 / 8 / 2.8404e-12 |

The zero kernel maximum includes the null-measure origin. Mesh, weight, EF,
eigenvalue, occupation and occupation-weighted projector continuity residuals
are all zero. Eigenvector unitarity maximum is 1.8874e-15.

The independent Ward source is the accepted SCF XC field, projected with the
same spherical-harmonic and radial metric map as the target Pauli density.
It is **not** constructed as Kxc times the Pauli density. The bare response is
applied unchanged. Norms are in the compact orthonormal response metric;
maximum residual is the largest absolute compact coefficient. The overlap is
`<m_P,response>/<m_P,m_P>`.

| eta (Ry) | Ward L2 | Ward relative L2 | Ward maximum | Normalized response overlap |
|---:|---:|---:|---:|---:|
| 0.01 | 0.175580908 | 0.251370583 | 0.108779774 | 0.871289738 + 0.0572422205i |
| 0.005 | 0.171219814 | 0.24512702 | 0.108953517 | 0.874393346 + 0.0287638308i |
| 0.0025 | 0.170101824 | 0.243526447 | 0.108997081 | 0.875179921 + 0.0144004664i |
| 0.00125 | 0.169820491 | 0.243123677 | 0.10900798 | 0.875377288 + 0.00720257673i |
| 0.000625 | 0.169750042 | 0.243022818 | 0.109010706 | 0.875426676 + 0.00360158209i |

These finite-eta samples approach a stable nonzero residual. They do not prove
the exact eta=0 identity. The Ward result is diagnostic only and never alters
Kxc, the field, bare response, denominator, poles or any interaction scalar.
No representation-closure tolerance was fitted to Fe. The unresolved large
residual prevents PASS.

| eta (Ry) | Denominator min singular value | Condition | Min abs eigenvalue | Dyson residual F | Relative | Infinity |
|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 0.0553459471 | 23.4923176 | 0.0736319709 | 6.44455208e-15 | 1.56246629e-15 | 3.30209888e-15 |
| 0.005 | 0.0336420937 | 38.6607909 | 0.044815229 | 1.18566972e-14 | 2.85798456e-15 | 4.01701432e-15 |
| 0.0025 | 0.0253546646 | 51.3017352 | 0.0337873873 | 1.47522576e-14 | 3.55044001e-15 | 6.72821663e-15 |
| 0.00125 | 0.0228078109 | 57.0315806 | 0.0303962369 | 1.45513413e-14 | 3.50070659e-15 | 5.77389178e-15 |
| 0.000625 | 0.0221247093 | 58.7927417 | 0.0294865336 | 1.3429783e-14 | 3.23056693e-15 | 5.6060971e-15 |

All five direct LAPACK solves succeeded. The maximum denominator singular value
is 1.30020–1.30077. The raw denominator remains finite along this ladder; no
zero singular value is required at finite eta. Interacting q/-q covariance
residual is 7.0793e-14, with no pole adjustment.

Separate real-space potential adapter smoke: state source
`accepted_realspace_potential_adapter`; fixed accepted EF = -0.37161 Ry;
reciprocal occupation count = 2, SCF target = 8, difference = -6 electrons.
An independent sum of emitted k weights and occupations also gives 2.
This significant state-representation limitation is recorded separately and
is not a reason to change the physical ALSDA kernel.

## B/C/D — tests and unresolved regressions

Commands used:

```sh
cmake --build build-lr-method-03 -j 4
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 ctest --test-dir build-lr-method-03 --output-on-failure -j 4
python3 tests/validation/lr_method_05r.py --binary build-lr-method-03/bin/rslmto.x --scratch build-lr-method-03/Testing/05r-clean --selectors
```

The broad 251-test command completed 134 tests before it was stopped after the
05R blocker to avoid expanding into unrelated material campaigns. All remaining
LR-named tests and the three existing LR material gates were then explicitly
run. A full 251-test PASS is **not** claimed.

- **B:** independent SR functional construction, Pauli negative control and no-Pauli
  input test pass. The old hybrid ratio differs measurably.
- **B:** direct local operator action and nonconstant repaired-kernel compact
  product quadrature/contraction pass; maximum compact oracle discrepancy is zero
  at printed precision (tolerance 2e-11).
- **B/C:** Dyson, compact mapping, Jülich, Mills, rotation and native RSGF algebra
  tests pass. Jülich/Mills production integrations pass unchanged.
- **C:** executed projected bare `d`, `spd`, `both` runs emit exactly their selected
  projections, respectively. All three pass.
- **A/D:** `Val22LrRadialGroundState` and `Dresp01BccFeMaterialGate` pass.
  `Val23LrPauliProjection` fails because its existing default 16-step real-space
  case does not reach outer SCF convergence. That deck is unchanged; no new
  state or convergence tuning was introduced.
- **D:** `LrBaselineNativeRsgf` fails its archived **interacting ALSDA** loss:
  old/reference -0.008699136480894; repaired -0.006887001292854.
  Its bare chiKS trace is unchanged to numerical precision:
  0.0306539533612394 + 0.0293102421139956i.
  A fresh build of the exact starting SHA passes the original archived reference.
  The changed loss is downstream of the repaired ALSDA interaction; native
  RSGF formulas and the archived reference remain untouched. No new RSGF
  material certification is claimed.
- **D:** the rotation example times out in the loaded broad CTest run. Its
  unchanged script passes with a 900-second timeout in separate scratch
  (three response rows). The ordinary timeout result is retained in the record.
- **D:** six Fe energy regressions fail both the repaired build and a fresh
  build of the exact starting SHA: `Lanczos`, `bccFe_block_fast_sp`,
  `bccFe_block_fast_dp`, `bccFe_chebyshev_fast_hoh`,
  `bccFe_chebyshev_legacy_hoh`, `bccFe_chebyshev_fast_ccor_2c`.
  Their references and tolerances were not changed.
- The broad run also encountered unrelated `Val12LmtoFieldsTorques` failure:
  global-rotation torque discrepancy 6.626e-3 T. It was not triaged or corrected
  in 05R; no attribution to the ALSDA change is asserted.
- `UnitLrEigenpairGfQuadrature` and `UnitLrEigenpairGfMixedEigenvectors` are
  already disabled in the configured suite and remain disabled.

For the 124 tests consisting of LR-prefixed tests plus the three LR material
gates: configured runs yield 119 passing, two failing, one timeout and two
disabled. The separate longer-timeout rotation rerun passes. The new selector
runs and the final four ALSDA/projection oracle checks also pass. There is no
category-E fixture used to certify the repair.

Machine-readable accepted-state evidence:
`build-lr-method-03/Testing/05r-clean/material_record.json`.
Raw Ward/Dyson artifact:
`build-lr-method-03/Testing/05r-clean/material/alsda_05r.dat`.
Test and baseline logs are copied under
`build-lr-method-03/Testing/05r-evidence/`.

05R stops with the physical kernel repaired and response closure unresolved.
