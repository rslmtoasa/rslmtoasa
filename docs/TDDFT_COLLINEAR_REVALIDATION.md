# TDVAL-01R — production revalidation of the rebuilt collinear TD-DFT stack

**Campaign date:** 2026-09-11
**Branch:** `fable_v4`
**Starting HEAD:** `795b1e4` (`docs: close native RSGF capability gates`)
**Status:** **PENDING** for R4 material validation; the preflight has been reopened after TDRUN-02.

## TDVAL-01R rerun after TDRUN-02

The authoritative [`RSGF_CAPABILITY_CLOSURE.md`](RSGF_CAPABILITY_CLOSURE.md)
now records:

```text
R0 = PASS
R1 = PASS
R2 = PASS
R3 = PASS
R4 = PENDING
```

The native-RSGF production material adapter is **AVAILABLE** through the normal
accepted-ground-state driver lifecycle. The reciprocal production smoke,
native production smoke, native service/driver equivalence, feature-off path,
and unsupported SOC/generalized-overlap rejection checks all pass.

This rerun therefore clears the former R0–R3 preflight stop and hands the
campaign to the Fe/Ni material-validation stages. No TDVAL PASS is claimed:
the full Fe/Ni convergence, static-invariant, route-separation, and literature
evidence remains to be generated, so R4 stays **PENDING**.

The pre-TDRUN stopped assessment retained below is historical evidence of the
earlier gate state; its native-RSGF `BLOCKED` wording is superseded by the
closure ledger and this rerun handoff.

This section supersedes the historical TDVAL-01 report below. It records the
production revalidation decision required by TDVAL-01R, including the evidence
actually rerun from the live branch. The campaign stops when a required
capability has a final `BLOCKED` verdict. No response equation, kernel,
Green-function transformation, Goldstone correction, or empirical parameter
was changed.

The later [`RSGF_CAPABILITY_CLOSURE.md`](RSGF_CAPABILITY_CLOSURE.md) supersedes
the historical LR-REP-00 native-response blocker for R0–R2: representation,
endpoint augmentation, and finite/provider bare-response capability are pass
within their documented scopes. TDRUN-02 now also records R3 production-driver
registration as pass. This report remains **PENDING** for TDVAL R4 material
validation.

## TDVAL-01R — production revalidation

### 1. Preflight and stopping decision

The initial worktree had no tracked modifications and contained only the
supplied untracked post-LR03 prompt-pack files. The exact initial state was:

```text
## fable_v4...origin/fable_v4
?? docs/dev/RS_LMTO_TDDFT_post_LR03/00_POST_LR03_MASTER_BLUEPRINT.md
?? docs/dev/RS_LMTO_TDDFT_post_LR03/11_RSGF-01_REALSPACE_GF_RESPONSE_BACKEND.md
?? docs/dev/RS_LMTO_TDDFT_post_LR03/12_NC-00_NONCOLLINEAR_FEASIBILITY.md
?? docs/dev/RS_LMTO_TDDFT_post_LR03/13_SV-00_STERNHEIMER_FEASIBILITY.md
?? docs/dev/RS_LMTO_TDDFT_post_LR03/README.md
?? docs/dev/RS_LMTO_TDDFT_post_LR03_prompt_pack.zip
```

Final verdicts were read from the evidence documents, not inferred from
filenames:

| prerequisite | evidence | final verdict used by TDVAL-01R |
| --- | --- | --- |
| LR-GF-01 | [`LR-GF-01_GF_CONTRACT_EVIDENCE.md`](LR-GF-01_GF_CONTRACT_EVIDENCE.md) | **PASS** for the orthogonal reciprocal coefficient-space GF contract |
| LR-BASIS-00 | [`LR-BASIS-00_RADIAL_AUGMENTATION_EVIDENCE.md`](LR-BASIS-00_RADIAL_AUGMENTATION_EVIDENCE.md) | **PASS** for the stated baseline |
| LR-01 | [`LR_RADIAL_GROUND_STATE_AUDIT.md`](LR_RADIAL_GROUND_STATE_AUDIT.md) | **PASS** for the accepted scalar-relativistic collinear radial contract |
| LR-02R | [`LR_RESPONSE_BASIS_MAPPING.md`](LR_RESPONSE_BASIS_MAPPING.md) | **PARTIAL PASS**; exact scalar-relativistic response is blocked |
| SR→Pauli closure | [`LR_SR_PAULI_NUMERICAL_CLOSURE.md`](LR_SR_PAULI_NUMERICAL_CLOSURE.md) | **PASS** as a measured Fe closure diagnostic; no arbitrary physical threshold |
| LR-03 | [`LR_TDDFT_CONVENTIONS.md`](LR_TDDFT_CONVENTIONS.md) | **PASS** for locked no-SOC collinear conventions |
| LR-04 | [`LR_RESPONSE_SPACE_ALGEBRA.md`](LR_RESPONSE_SPACE_ALGEBRA.md) | **PASS** for canonical response-space algebra |
| LR-05 | [`LR_PAULI_TRANSITION_VERTEX.md`](LR_PAULI_TRANSITION_VERTEX.md) | **PASS** for the supported Pauli/no-SOC primitive and guards |
| LR-06 | [`LR_KS_SUSCEPTIBILITY.md`](LR_KS_SUSCEPTIBILITY.md) | **PASS** for spectral/Lehmann finite oracles |
| LR-GF-02 | [`LR_GF_SUSCEPTIBILITY_CROSSCHECK.md`](LR_GF_SUSCEPTIBILITY_CROSSCHECK.md) | **PASS** for the independent reciprocal-GF finite cross-check |
| KXC-01 | [`KXC_ALSDA_TRANSVERSE_KERNEL.md`](KXC_ALSDA_TRANSVERSE_KERNEL.md) | **PASS** for algebraic and finite-fixture evidence only |
| GSR-01 | [`GOLDSTONE_SUMRULE_INTERACTION.md`](GOLDSTONE_SUMRULE_INTERACTION.md) | **PARTIAL PASS**; service fixtures pass, deliberate full-response identity is blocked |
| GCR-01 | [`GOLDSTONE_EIGENVALUE_CORRECTION.md`](GOLDSTONE_EIGENVALUE_CORRECTION.md) | **PASS** for optional finite-matrix gates |
| TDDY-01 | [`TDDFT_DYSON_AND_LOSS.md`](TDDFT_DYSON_AND_LOSS.md) | **PASS** for finite-matrix Dyson/loss evidence |
| LR-REP-00 | [`LR_RS_GF_REPRESENTATION_AUDIT.md`](LR_RS_GF_REPRESENTATION_AUDIT.md) | **BLOCKED** for native RS physical radial/Pauli response promotion |
| RSGF-00 | [`RSGF_ENDPOINT_AUGMENTATION.md`](RSGF_ENDPOINT_AUGMENTATION.md) | **PASS** for the finite endpoint oracle only |
| RSGF-01R | [`TDDFT_RS_GF_BACKEND.md`](TDDFT_RS_GF_BACKEND.md) | **PARTIAL PASS** for finite/provider-oracle evidence; production material execution is not claimed |
| TDRUN-01 | [`TDDFT_PRODUCTION_DRIVER.md`](TDDFT_PRODUCTION_DRIVER.md) | **PASS** for reciprocal baseline and tiny lifecycle smoke; native RSGF is not exposed |

The LR-REP-00 final verdict is the stopping condition. RSGF-01R has a finite
response implementation and passing unit oracle, but its own claim boundary
excludes production-driver material execution. TDVAL-01R therefore does not
proceed to Fe/Ni response generation, convergence ladders, or literature
comparison.

### 2. Production-driver evidence

The focused prerequisite/closure run passed 26/26 tests, and the native finite
response oracle passed 1/1:

```text
100% tests passed, 26 tests passed, 0 tests failed
100% tests passed, 1 test passed, 0 tests failed
```

The tests include LR-01, LR-02 response basis, LR-04, LR-05, LR-06, LR-GF-02,
KXC-01, GSR-01, GCR-01, TDDY-01, TDRUN-01, capability guards, RSGF-00, and
`UnitLrRsGfSusceptibility`. `TddftProductionDriverSmoke` proves the actual
`SCF → accepted LR-01 snapshot → reciprocal eigenpair handoff → TDRUN →
LR-05/LR-06 → KXC-01 → TDDY-01` lifecycle on a tiny fixture, including
accepted-state non-mutation. It does not prove a converged Fe response. The
production `backend` contract currently accepts only `lehmann` and
`reciprocal_gf`.

The verification build was GNU Fortran 13.3.0, CMake 3.28.3, Ninja, Debug,
OpenMP enabled, MPI disabled, libXC enabled, unit tests enabled, and regression
tests disabled. The executable hash at verification was
`77e4f1a9c3d946c661b5e9f18582124c24ab6226851bd75544459e37448b57ea`.

### 3. Ground-state provenance

The canonical Fe foundation uses `tests/scf/cases/bulk/bccFe`, scalar-relativistic
collinear `nsp=1`, legacy Barth-Hedin XC, 300 K Fermi-Dirac occupations, an
8×8×8 full mesh, and the documented SR→Pauli diagnostic. Its reported
integrated SR and Pauli moments are 2.1038582410 and 2.1487266935 μB. This is
an explicit approximation diagnostic, not a tuning target.

The available Ni deck is
`results/validation/TDVAL-01_FE_NI/ground_state/fccNi/input.nml`: `alat=3.520
Å`, fcc, `wav=1.410`, `ct=5`, `r2=25`, 16×16×16 mesh, 300 K, and `ham_only`.
Its evidence explicitly says that no comparable accepted LR-01 Pauli closure
is present. It is not promoted to a TDVAL response reference.

No Fe/Ni TDVAL output was generated after the preflight stop. Per-run metadata
for this campaign is consequently **not available**: exact commit, input,
structure/ASA data, XC provenance, moment, Fermi state, temperature, meshes,
basis, response cutoff, q/ω/η, backend, interaction route, correction status,
projection metric, and build details.

### 4. Fe baseline

**BLOCKED / not run.** Fe LR-01 and SR→Pauli foundation evidence is available,
but no current TDVAL Fe response exists through all required bare backends.

### 5. Ni baseline

**BLOCKED / deferred.** A documented Ni deck exists, but no current accepted
LR-01/LR-02 Pauli response snapshot comparable to Fe is available.

### 6. Spectral versus reciprocal-GF

**[NUMERICAL CROSS-CHECK] PARTIAL PASS.** LR-06 and LR-GF-02 agree on the
complete finite response-space fixtures, including static q=0 and an off-mesh
folded-q case. No material comparison was run.

### 7. Spectral versus native-RSGF

**[NUMERICAL CROSS-CHECK] PARTIAL PASS for finite fixtures; BLOCKED for
material validation.** `UnitLrRsGfSusceptibility` passes full radial/angular
finite dense/provider, phase, and bubble oracles; RSGF-00 passes its endpoint
oracle. Native providers remain callback seams and are not selected by TDRUN.

### 8. Coefficient-GF / augmentation / bubble decomposition

| layer | finite evidence | Fe/Ni material status |
| --- | --- | --- |
| coefficient GF | dense complex-energy inverse/provider contract | **BLOCKED / not run** |
| endpoint augmentation | RSGF-00 dense, spectral, provider, and negative controls | **BLOCKED / not run** |
| real-space bubble/integration | full response-space native-versus-spectral oracle | **BLOCKED / not run** |

No susceptibility-only scalar is used to conceal an unmeasured coefficient or
endpoint error.

### 9. Recursion/Chebyshev convergence

**BLOCKED / not run.** The native block-recursion and Chebyshev types expose
independent controls, but no production material callback is wired into
TDRUN. The accepted terminator is not under audit; no three-depth or
three-order material ladder is claimed.

### 10. k-mesh convergence

**BLOCKED / not run.** No Fe/Ni response exists from which to track static
`chiKS`, Goldstone residuals, finite-q energies, or continuum features across
three meshes.

### 11. eta convergence

**BLOCKED / not run.** No material peak positions, widths, static responses, or
backend differences were generated for three broadenings. No eta→0
extrapolation or undocumented width subtraction is made.

### 12. Radial/angular convergence

**BLOCKED / not run.** The production radial mesh remains the reference. For
`spd`, formal complete response content is `L_max=4`; a reduced cutoff would
be an approximation and was not called exact.

### 13. Raw Goldstone and interaction routes

**[ALGEBRAIC] PASS; [MATERIAL VALIDATION] BLOCKED.** Finite tests exercise raw
diagnostics and separate direct ALSDA, sum-rule, and optional BES objects.
GSR-01 retains an explicit blocked full-response identity diagnostic. No Fe/Ni
q=0 residual, rigid-rotation overlap, denominator eigenvalue, or correction
effect is available.

| route | material status |
| --- | --- |
| direct ALSDA / KXC-01 | **DEFERRED** |
| Lounis sum-rule / GSR-01 | **DEFERRED** |
| optional BES / GCR-01 | **DEFERRED**; production selection remains gated |

### 14. Long-wavelength magnons

**BLOCKED / not run.** No production q sequence, loss-peak extraction,
quadratic fit, stiffness, or fit-window sensitivity is claimed.

### 15. Spectral function and loss matrix

**BLOCKED / not run for materials.** TDDY-01 finite tests validate loss
algebra and near-pole diagnostics, but no Fe/Ni Stoner continuum, collective
branch, or intrinsic-versus-eta damping result was generated.

### 16. q↔−q covariance and BZ gauge

**[ALGEBRAIC] PASS for finite contracts; [MATERIAL VALIDATION] BLOCKED.**
LR-03/LR-06/TDDY-01 cover covariance and folded-q fixtures, and RSGF covers
its finite phase. Required converged Fe/Ni commensurate and arbitrary/folded
production q checks were not run. No odd-in-q material component was
symmetrized away.

### 17. Literature comparison

**[LITERATURE COMPARISON] DEFERRED.** No converged Fe/Ni TDVAL result exists
to compare with compatible Halle/BES lineage or experimental references. No
parameter was tuned toward literature.

### 18. Failure attribution

The blocker belongs to the RSGF/representation integration boundary, not
TDVAL. The owning next layer is:

```text
native block/Chebyshev provider
    -> accepted-state production callback
    -> RSGF-00 Pauli endpoint adapter
    -> RSGF-01R bubble
    -> TDRUN backend selection and provenance
```

TDVAL must be reopened only after LR-REP-00/RSGF-01R have a non-blocked
production verdict and the reciprocal baseline is generated from the same
accepted state. TDVAL itself is not patched to bypass the blocker.

### 19. Claims discipline

- **[ALGEBRAIC]** conventions, identities, and service guards;
- **[NUMERICAL CROSS-CHECK]** independent finite spectral/GF comparisons;
- **[MATERIAL VALIDATION]** unavailable for Fe/Ni in this stopped run;
- **[LITERATURE COMPARISON]** deferred because material validation is absent.

No claim of “TD-DFT validated” is made.

### 20. Final per-route verdict

```text
spectral bare response:            PARTIAL PASS  (finite/algebraic only)
reciprocal-GF bare response:       PARTIAL PASS  (finite cross-check only)
native-RSGF bare response:         BLOCKED       (no production material route)

direct ALSDA interacting response: DEFERRED      (no Fe/Ni material run)
sum-rule interacting response:     DEFERRED      (no Fe/Ni material run)
optional BES-corrected response:   DEFERRED      (material route not selected)

bcc Fe:                            BLOCKED
fcc Ni:                            BLOCKED
production lifecycle:              PARTIAL PASS  (tiny reciprocal smoke only)
overall certified baseline:        BLOCKED       (TDVAL-01R material validation)
```

The native-RSGF blocker does not invalidate the independent finite reciprocal
cross-check. Conversely, the reciprocal finite cross-check does not validate
native-RSGF production material execution.

## Commit record

Important blockers remain; the evidence is committed with a message that does
not imply successful validation:

```text
docs: record TDVAL-01R preflight blocker
```

## Historical TDVAL-01 report (superseded)

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

The old production TD-DFT backend and its validation runner were intentionally
removed by the clean-room purge.  TDRUN-01 now supplies a clean production
SCF/eigenpair-to-LR-03…TDDY-01 material adapter through
[`source/tddft_production_driver.f90`](../source/tddft_production_driver.f90)
and the `post_processing='tddft'` hook in
[`source/calculation_preprocessing.f90`](../source/calculation_preprocessing.f90).
Its tiny magnetic lifecycle smoke test passes, so the material campaign is
unblocked.  The old `post_processing='susceptibility'` route remains rejected;
only the new minimal `&tddft` contract is accepted.

No current bcc-Fe or fcc-Ni response has yet been validated through the new
production path.  The historical pre-purge Fe/Ni spectra are explicitly
excluded and are not used to rescue this result.

## Evidence categories

| category | result | evidence boundary |
| --- | --- | --- |
| algebraic consistency | **PASS** | LR-03 convention identities and LR-04…TDDY-01 finite fixtures |
| independent numerical cross-check | **PASS** | spectral/GF bare-response fixture and current closure tests |
| production material adapter | **PASS** | TDRUN-01 accepted-state handoff and tiny SCF-to-response smoke path |
| numerical convergence of Fe/Ni response | **BLOCKED** | TDVAL rerun through the new production path is pending |
| physical material validation | **BLOCKED** | no Fe/Ni production `chiKS`, kernel, Dyson, loss, or mode evidence |
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
| executable hash | `build/bin/rslmto.x`: `eb4d0e506edb3f9f4e5959549428912a181b13d39b5bf23613ad5c6d0fbfc54b` |

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

These modules are now connected by the clean TDRUN-01 production adapter.  The
adapter accepts only the validated reciprocal collinear baseline: no SOC,
`ham_only`, orthogonal, second-order, `sp/spd`, and no unsupported additive
operators.  It consumes the accepted LR-01 snapshots and accepted Fermi state;
it does not reconstruct SCF state or search for a different Fermi level.  The
native RSGF backend remains gated pending RSGF-01.

## Checks run

The prerequisite and closure checks were rebuilt and executed from the live
branch:

```text
python3 tests/validation/val24_lr03_conventions.py
val24_lr03_conventions: PASS (algebraic conventions only)

ctest --test-dir build --output-on-failure -R '^(UnitLrBasisRadial|UnitLrRadialGroundState|UnitLrBasisAngular|UnitLrBasisAugmentation|UnitLrResponseBasis|UnitLrResponseSpace|UnitLrPauliTransitionVertex($|Reject)|UnitLrKsSusceptibility|UnitLrGfSusceptibility|UnitLrAlsdaKernel($|Reject)|UnitLrGoldstoneSumrule|UnitLrGoldstoneCorrection|UnitTddftDyson|UnitTddftProductionDriver|TddftProductionDriverSmoke)$'
100% tests passed, 15 tests passed, 0 tests failed

ctest --test-dir build --output-on-failure -R '^UnitLr(AlsdaKernelReject|PauliTransitionVertexReject)'
100% tests passed, 7 tests passed, 0 tests failed

ctest --test-dir build --output-on-failure -R '^TddftProductionDriverFeatureOff$'
100% tests passed, 1 test passed, 0 tests failed
```

The 15-test run included the LR basis and response-space units, LR-05, LR-06,
LR-GF-02, KXC-01, GSR-01, GCR-01, and TDDY-01.  The seven expected-failure
tests covered ALSDA provenance and zero-magnetization guards plus the five
LR-05 capability guards.  The TDRUN-01 additions cover feature-off parsing,
complete direct service-versus-driver reproducibility, accepted-state
non-mutation, expected SOC/generalized-overlap rejection, and the tiny
SCF-to-response smoke path.  The separate LR-01/LR-02N validation gates remain
the prerequisite evidence recorded above; their `Val22`/`Val23` CTest entries
are not registered in this build.  The separate feature-off executable
regression also passes on the same tiny ordinary SCF fixture.

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

The former production integration blocker is cleared by TDRUN-01.  TDVAL-01
must now be rerun through the new production path, not repaired by changing the
LR-03 equations or applying empirical normalization.  The rerun must provide,
for both Fe and Ni:

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
- [x] clean production material adapter exists and tiny lifecycle smoke passes;
- [x] explicit algebraic, numerical, material, and literature evidence categories used;
- [x] no empirical tuning or physics-equation changes;
- [x] validation report written;
- [ ] clean current bcc-Fe production response reference;
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

`td-dft: wire clean post-SCF response driver`
