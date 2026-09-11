# TDVAL-01 — rebuilt collinear TD-DFT validation

**Campaign date:** 2026-09-11
**Status:** **BLOCKED for physical Fe/Ni validation**
**Claim level:** implementation-level algebraic and finite-fixture numerical
evidence is current; no material spectrum, stiffness, damping, or literature
agreement is claimed.

This report follows `docs/dev/RS_LMTO_TDDFT_post_LR03/09_TDVAL-01_COLLINEAR_VALIDATION.md`, `00_POST_LR03_MASTER_BLUEPRINT.md`, and `README.md` from the supplied post-LR-03 prompt pack.

## Executive result

The rebuilt LR-03 through TDDY-01 library contracts pass their current
implementation tests, including the independent reciprocal-GF cross-check,
direct ALSDA, sum-rule, optional BES, and canonical Dyson/loss fixtures.  The
live LR-01 radial and LR-02N Pauli ground-state gates also pass.

The required material campaign cannot be run on this branch.  The old
production TD-DFT backend and its validation runner were intentionally removed
by the clean-room purge.  The current executable explicitly rejects both an
old `&tddft` input and `post_processing='susceptibility'` in
[`source/calculation.f90`](../source/calculation.f90#L282).  The rebuilt
modules are library primitives with caller-supplied response spaces, radial
bases, endpoint eigenstates, susceptibilities, and interactions; no production
SCF/eigenpair-to-LR-03…TDDY-01 material adapter is present.

Consequently, no current bcc-Fe or fcc-Ni response was evaluated.  The
historical pre-purge Fe/Ni spectra are explicitly excluded and are not used to
rescue this result.

## Evidence categories

| category | result | evidence boundary |
| --- | --- | --- |
| algebraic consistency | **PASS** | LR-03 convention identities and LR-04…TDDY-01 finite fixtures |
| independent numerical cross-check | **PASS** | spectral/GF bare-response fixture and current closure tests |
| numerical convergence of Fe/Ni response | **BLOCKED** | no current production material response exists |
| physical material validation | **BLOCKED** | no Fe/Ni `chiKS`, kernel, Dyson, loss, or mode outputs |
| literature comparison | **NOT PERFORMED** | closest references recorded below; no numeric comparison is supportable |

No physics equation, sign, factor, kernel scale, moment, broadening, or
frequency shift was changed during this campaign.

## Starting provenance

The campaign started on branch `fable_v4` at exact HEAD
`e22d81d422929f4314ef55ebea9fa9f5cdca79ff`.  The worktree was dirty only
because the supplied post-LR-03 prompt-pack files were untracked; no tracked
source or test file was modified before the campaign checks.  The validation
environment was:

| item | value |
| --- | --- |
| compiler | GNU Fortran 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) |
| CMake | 3.28.3 |
| configured build | `build`, unit and regression tests enabled, libXC enabled, MPI disabled, OpenMP enabled |
| executable hash | `build/bin/rslmto.x`: `dc52d3c4db08198673107ca8485cc0870baded315bc7b8090362887f640f13ed` |

The final commit records this report and the completion record in the TDVAL
specification.  The adjacent supplied prompt-pack files remain user-owned
untracked inputs.

## Reference-state audit

No state below is accepted as a TDVAL material response reference.  The table
records the available inputs and why they cannot satisfy the rebuilt response
contract.

| material | available input and structure | magnetic/XC/electronic provenance | disposition |
| --- | --- | --- | --- |
| bcc Fe | [`tests/scf/cases/bulk/bccFe/input.nml`](../tests/scf/cases/bulk/bccFe/input.nml): `alat=2.86120`, `wav=1.40880`, `ct=3`, `r2=9`, `lmax=2` | The accepted LR-02N foundation closure is separately documented at commit `2793c1c47fd6d2669184ac53707d2186d6200919`: scalar-relativistic collinear `nsp=1`, 8³ k mesh, 300 K, `EF=-0.069249685116850188 Ry`, legacy Barth-Hedin, no additive terms; integrated SR moment `2.1038582410`, Pauli projection `2.1487266935` | ground-state/projection evidence only; no rebuilt LR response snapshot or material transition data |
| fcc Ni | [`results/validation/TDVAL-01_FE_NI/ground_state/fccNi/input.nml`](../results/validation/TDVAL-01_FE_NI/ground_state/fccNi/input.nml): `alat=3.520 Å`, `wav=1.410`, `ct=5`, `r2=25`, periodic fcc, `16³` mesh, `strux_lib`, `lmax=2` | Deck requests SOC-off, `ham_only`, 300 K, automatic Fermi level, and records input `fermi=0.101023 Ry`; the tracked artifact does not contain a current accepted LR-01 Pauli closure or a material eigenpair snapshot | legacy ground-state provenance only; no rebuilt LR response snapshot or material transition data |

The Fe SR→Pauli numbers are a foundation numerical input, not TD-DFT
validation.  The LR-02N evidence explicitly defers a comparable accepted Ni
closure.  The available Ni `Ni.nml` records orientation and radial ground-state
parameters but does not provide a current accepted response-state moment.

## Current production boundary

The clean-room purge records removal of the old `source/tddft_backend.f90`,
old `source/tddft_chi0*.f90`, old mode/Goldstone paths, the old validation
runner, and old response outputs in the [purge manifest](TDDFT_CLEANROOM_PURGE_MANIFEST.md).
The current source graph instead contains the rebuilt LR modules:

- LR-04 response-space metric/operator algebra;
- LR-05 Pauli transition vertices;
- LR-06 spectral/Lehmann `chiKS`;
- LR-GF-02 reciprocal-GF `chiKS` cross-check;
- KXC-01 direct radial ALSDA;
- GSR-01 sum-rule interaction;
- optional GCR-01 BES correction; and
- TDDY-01 canonical Dyson and loss matrix.

These modules are individually callable, but the main calculation path still
fails closed at [`calculation.f90:282`](../source/calculation.f90#L282),
[`calculation.f90:297`](../source/calculation.f90#L297), and
[`calculation.f90:380`](../source/calculation.f90#L380) for TD-DFT input.
This is an integration/availability blocker, not evidence of a physics failure
in any one of the rebuilt finite fixtures.

## Checks run

The prerequisite and closure checks were rebuilt and executed from the live
branch:

```text
python3 tests/validation/val24_lr03_conventions.py
val24_lr03_conventions: PASS (algebraic conventions only)

ctest --test-dir build --output-on-failure -R '^(Val22LrRadialGroundState($|Libxc$)|Val23LrPauliProjection|UnitLrBasisRadial|UnitLrRadialGroundState|UnitLrBasisAngular|UnitLrBasisAugmentation|UnitLrResponseBasis|UnitLrResponseSpace|UnitLrPauliTransitionVertex($|Reject)|UnitLrKsSusceptibility|UnitLrGfSusceptibility|UnitLrAlsdaKernel($|Reject)|UnitLrGoldstoneSumrule|UnitLrGoldstoneCorrection|UnitTddftDyson)$'
100% tests passed, 16 tests passed, 0 tests failed

ctest --test-dir build --output-on-failure -R '^UnitLr(AlsdaKernelReject|PauliTransitionVertexReject)'
100% tests passed, 7 tests passed, 0 tests failed
```

The 16-test run included `Val22LrRadialGroundState`,
`Val22LrRadialGroundStateLibxc`, `Val23LrPauliProjection`, the LR basis and
response-space units, LR-05, LR-06, LR-GF-02, KXC-01, GSR-01, GCR-01, and
TDDY-01.  The seven expected-failure tests covered ALSDA provenance and
zero-magnetization guards plus the five LR-05 capability guards.

The current GF fixture reports full-response, finite-fixture evidence including
analytic GF error `1.1277e-09`, q=0 static GF difference `9.2484e-10`, and
resolution samples with absolute differences `2.1443e-09` and `3.6257e-09`.
These are response-space numerical checks, not Fe/Ni material results.

## Validation hierarchy

### Stage A — static invariants

| invariant | result |
| --- | --- |
| LR-03 circular/covariance algebra | **PASS**, material-independent `val24` and unit fixtures |
| q=0 raw Ward/Goldstone residual | **PASS as a diagnostic contract only**, with synthetic residuals and explicit blocker behavior |
| rigid-rotation overlap | **PASS in finite fixtures only** |
| direct ALSDA versus sum-rule interaction | **PASS as independent finite-fixture route separation only** |
| optional BES route | **not selected for material validation**; standalone optional fixture passes |
| spectral versus reciprocal-GF static agreement | **PASS on the certified finite response fixture** |

No Stage-A item was evaluated on a current Fe or Ni response state.

### Stages B–E — material convergence

The radial/angular, k-mesh, eta, and q ladders are **BLOCKED**.  There is no
current production response output from which to measure static response norms,
Goldstone residuals, finite-q spectra, peak positions, broadenings, BZ folding,
or arbitrary-q covariance.  No reduced angular cutoff is called exact, no eta
limit is extrapolated, and no material q² fit is attempted.

## Interaction routes and collective response

The three material routes remain separate by contract:

1. direct radial ALSDA;
2. Lounis/LCMM sum-rule interaction; and
3. optional BES-corrected ALSDA.

None was run on Fe or Ni.  In particular, a passing finite-fixture GSR or GCR
test does not validate direct ALSDA, and no corrected or sum-rule fixture is
used to rescue the absent raw material route.  There are no current material
loss matrices, collective branches, damping widths, or stiffness values.  The
ferromagnetic `omega(q)=Dq²+O(q⁴)` requirement therefore remains unassessed.

## Literature comparison

The closest compatible methodological references are:

- [Buczek, Ernst, and Sandratskii, Phys. Rev. B 84, 174418 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.84.174418), which studies bulk bcc Fe and fcc Ni within ab initio linear-response TD-DFT/AL(S)DA and a KKR Green-function implementation, including Landau damping;
- [Lounis, Costa, Muniz, and Mills, Phys. Rev. B 83, 035109 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.83.035109), which defines a real-space KKR-GF TDDFT route and Goldstone-preserving interaction construction, but applies it to Cu-supported adatoms/dimers rather than the present bulk Fe/Ni target; and
- [Pajda et al., Phys. Rev. B 64, 174402 (2001)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.64.174402), a nonrelativistic TB-LMTO Green-function exchange/stiffness reference, retained only as context and not treated as a TD-DFT comparison.

The intended current setup differs from the closest bulk TDDFT reference in
basis/backend (RS-LMTO-ASA versus KKR-GF), response representation (direct
radial LR-04 space), and the explicitly documented scalar-relativistic-ground-
state to Pauli/no-SOC response approximation.  Since no current material
response was generated, no numerical literature comparison is made and no
eta, kernel scale, or moment is tuned toward a published curve.

## Failure triage and unblock condition

The observed failure is at the production integration boundary.  It should be
returned to the response-driver/interface work, not repaired in TDVAL-01 and
not addressed by changing the LR-03 equations or applying empirical
normalization.  TDVAL-01 can be reopened only after a current production path
provides, for both Fe and Ni:

- accepted LR-01 radial states and explicit Pauli magnetization provenance;
- matching `ham_only`, orthogonal, collinear, no-SOC reciprocal eigenpair
  snapshots at k and exact folded k+q;
- material calls to LR-05/LR-06, KXC-01, GSR-01, the explicitly selected GCR
  route if desired, and TDDY-01; and
- reproducible output for the Stage A–E ladders and route-separated reports.

## Completion checklist

- [x] branch, starting HEAD, worktree state, compiler, and binary recorded;
- [x] LR-03 through TDDY-01 implementation prerequisites re-audited;
- [x] LR-01/LR-02N radial and Pauli closure gates rerun;
- [x] independent spectral/GF finite-fixture evidence rerun;
- [x] historical pre-purge material response evidence excluded;
- [x] current production TD-DFT blocker verified from the live source and documented;
- [x] explicit algebraic, numerical, material, and literature evidence categories used;
- [x] no empirical tuning or physics-equation changes;
- [x] validation report written;
- [ ] clean current bcc-Fe response reference;
- [ ] clean current fcc-Ni response reference;
- [ ] Fe/Ni static invariants;
- [ ] radial/angular convergence;
- [ ] three-mesh k convergence;
- [ ] three-eta convergence and damping separation;
- [ ] commensurate and arbitrary-q production validation;
- [ ] direct ALSDA material result;
- [ ] sum-rule material result;
- [ ] optional corrected material result, if selected;
- [ ] converged quadratic branch and stiffness;
- [ ] quantitative literature agreement;
- [ ] TDVAL-01 PASS.

The unchecked items are deliberate and are the reason the campaign status is
BLOCKED rather than PASS.

## Commit

`tests: validate rebuilt collinear TDDFT response`
