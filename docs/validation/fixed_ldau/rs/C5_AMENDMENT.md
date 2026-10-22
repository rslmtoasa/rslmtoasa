# C5 — RS LDA+U: issue corrections, RS robustness triage, impurities (Stage RD)

**Session setup.** Sol is the orchestrator. Sol plans, delegates every execution task to Luna agents, reviews their output, and writes the final report. Sol does not run campaigns itself.

Luna agents are used as follows:

- **Workers** execute one well-defined task each. Independent runs may go in parallel.
- **One auditor** recomputes every GATE verdict and its arithmetic from the raw outputs before Sol records it, and checks each run's preflight (section 2.6).
- **Workers never change criteria, tolerances or scope.** If a worker believes a criterion is wrong, it reports to Sol, and Sol applies section 2.5.

This amendment is user-approved. Record it verbatim in `docs/validation/fixed_ldau/rs/C5_AMENDMENT.md`. It is self-contained: everything from earlier amendments that still applies is restated here.

---

## 1. Context

### Repository and branches

- Repository: `rslmtoasa/rslmtoasa`, local clone at `/Users/andersb/Jobb/rslmto_fable_v4b`.
- Integration branch: `fable_v4b_u`, pushed tip `908074f`.
- Create the work branch `fix/ldau-rs-impurity` from that tip.
- **Never push anything.** At the end, merge locally into `fable_v4b_u` (section 7) and report the merge SHA. The user reviews and pushes.
- `cleanup/` holds user-owned protocol files. Read them; never edit or commit them.
- `tests/gbt/` belongs to a different effort. Never add it.

### Scientific scope (fixed)

**Primary target:** fixed-U LDA+U in real space (RS recursion) for:

- bulk (`calctype='B'`)
- impurities (`calctype='I'`)
- surfaces (`calctype='S'`)

These are the code's unique capabilities.

**Convenience path:** KS bulk. **A KS bulk host followed by an RS impurity or surface run is a supported workflow.** KS impurity is not supported and must not be pursued.

**Guarded scope, do not touch:**

- SOC or noncollinear spin (`nsp>=2`)
- shells other than d, and `lmax/=2`
- HOH
- self-consistent U, intersite V, ACBN0

### Physics and terminology

**Local radial moments.** With `hubbard_radial_moments='local'`, the radial (atomic-sphere) step receives the band centre and second moment of the LDA part of the Hamiltonian:

- c′ = c − Tr(V_in P)/q0
- q2′ = q2 − 2 Tr(V_in M1ᶜ) + Tr(V_in² P) − t²/q0

This is implemented in `math.f90::correct_d_shell_moments`. Here:

- V_in is the Hubbard d block actually inserted into the Hamiltonian in this iteration;
- P is the occupation matrix;
- M1 is the energy-weighted occupation matrix.

In `'spectral'` mode the Hubbard shift leaks into the radial functions instead.

**Energy candidates.**

| Candidate | Definition |
|---|---|
| a | E_spectral − Tr(V_in P) + E_Hub |
| b | E_local + E_Hub |
| c | E_spectral + E_Hub |

The Hellmann–Feynman (HF) theorem is the arbiter. On a converged state of the correct functional, dE/dX = ∂E_Hub/∂X at fixed P, for X = U or J.

**RS SCF functional.** RS uses step-function occupations, so the functional is E, not E − TS.

**Base ASA non-stationarity.** At U=0, a fixed d-level shift λW gives dE/dλ ≠ Tr(WP) by −0.14 % in KS and −0.16 to −0.17 % in RS. It is a property of the ASA energy assembly and sets the precision floor of every HF test.

**Mode consistency rule.** Potential files record the radial-moment mode and the producing backend for U-carrying atoms. Every U-carrying atom in a calculation, whether a frozen host atom or a recursion atom, must use the same mode, and a mismatch is fatal. The backend may differ, since KS-produced and RS-produced files are interchangeable.

### Current status

| Configuration | Local moments | Energy (E_Hub routed) | `auto` resolves to |
|---|---|---|---|
| KS bulk | validated | validated (candidate b) | `local` |
| RS bulk | validated (RB0–RB4) | validated (C4: b passes U, J=0.5, J=0.9) | `local` |
| RS impurity (`'I'`) | not implemented for recursion atoms | none | `spectral` |
| RS surface (`'S'`) | not implemented | none | `spectral` |

### Relevant earlier findings

- **RS acceptance and stopping.** The RS outer residual plateaus around 1–2×10⁻⁸. RS-specific fixed-U stopping tolerances are in production: P ≤ 5×10⁻⁷, V ≤ 10⁻⁷ Ry, streak 3. KS is unchanged.
- **Embedding records (spectral mode).**
  - Fe impurity in an Fe host, U=0: matches RS bulk to about 10⁻⁵.
  - Fe, U=3.8 eV, J=1.1 eV: differs from bulk by **−6.6×10⁻³ μB** (about 735× the U=0 control).
  - In that Fe case, the Hubbard operators, Fermi level and Green's-function evaluation are identical between impurity and bulk. The recursion coefficients differ.
  - Cluster recipe: rc=80 for the bulk reference, rc=100 for the impurity.
  - **Ni** (U=3 eV, J=0.7 eV) embeds cleanly: its U>0 discrepancy is smaller than its U=0 control.
- **Pre-existing RS failures, all reproducible on clean `8b8f7ef` (not LDA+U):**
  - **hcp Co impurity setup:** stops before SCF on a bounds violation at `source/lattice_cluster.f90:996` ("dimension mismatch 3/50000").
  - **Co bulk at U=0** (RB6 records): goes non-finite after roughly 19–40 iterations.
  - **Fe, rc=100:** the fixed-EF impurity-style deck goes non-finite after 102 iterations at U>0 and after 34 at U=0, with pre-LU Q/mixed-q NaNs and an invalid PL. The rc=100 U=0 bulk never met acceptance in 600 steps.
  - **Fe rc=140 runs:** INVALID_SETUP (fresh atomic start, diverging trajectory). Not an instability.
- **Production recursion** has an FP32 fast path.

### Key documents

- `docs/validation/fixed_ldau/rs/`: the C1_RA4, C1_RB, C3_RC, C3_PART_A, C4_* results and their `*_PROVENANCE.json` sidecars. These hold the exact decks for the Fe, Ni and B2FeCo cases. Reuse them.
- `tests/KNOWN_ISSUES.md`
- `docs/validation/fixed_ldau/ENERGY_ACCOUNTING.md`, `REMEDIATION_GATES.md`
- `docs/source/theory/lda_u_workflow.rst`
- Fixtures: `tests/scf/cases/bulk/bccFe_hubbard_nsp1`, `tests/scf/cases/impurity/B2FeCo_nsp1`.

---

## 2. Standing rules

### 2.1 RS-accepted state

All of the following must hold:

- (a) outer Diff ≤ 3×10⁻⁸;
- (b) over the last 20 complete iterations, the spread (max − min) of the d occupation per spin is ≤ 5×10⁻⁷, and the spread of the sphere moment is ≤ 10⁻⁶ μB;
- (c) fixed-U P/V residuals are logged, not gated.

### 2.2 Fallback ladder for an SCF that is not accepted at its step cap

1. Continue from the last state for 200 steps with beta halved.
2. Then continue for 200 steps with `broyden` and beta = 0.01.
3. Then label the state PROVISIONAL, use it for RECORDS only, and move on.

Interrupted or retried attempts count toward the same budget.

### 2.3 RECORDS versus GATES

- RECORDS are measured and tabulated. They cannot fail, and they never block progress.
- Only items marked **GATE** can fail.

### 2.4 Crash triage

On a crash or a non-finite value:

1. Rerun once with the Debug, bounds-checked build and capture the file:line backtrace.
2. Rerun the same deck with the clean `8b8f7ef` binary.
3. Rerun with U=J=0.

If it also fails on `8b8f7ef` or at U=0, classify it as PRE-EXISTING and continue. STOP only if it is new relative to `8b8f7ef`, appears only with U>0, and lies in code changed since `8b8f7ef`.

### 2.5 STOP only if

- (a) no-U byte identity breaks on the standard fixtures (RS `bulk/bccFe`, KS `k_space_scf/bccFe`);
- (b) a GATE fails after one honest debugging attempt (state the hypothesis you tested);
- (c) section 2.4 says stop;
- (d) continuing would require touching guarded scope or changing a physics formula.

**If a criterion turns out to rest on a wrong error model:** finish every other item that doesn't depend on it, then report the measurement that exposes the problem together with a proposed corrected criterion. Never self-amend, and never loosen a tolerance after seeing results.

Everything else: write it down and keep going.

### 2.6 Preflight checklist (auditor verifies for every run)

- **Restart rule.** Start from the specified converged endpoint, never from fresh atomic potentials, unless a fresh start is explicitly requested.
- **Backend.** Confirm from the live log, including `use_kspace` inside the correct namelist.
- **Input parsing.** Every namelist is closed with `/` and the input parses, confirmed before counting steps.
- **Mixing.** The mixing parameters actually applied match the ladder step.
- **Resolved mode.** The `HUBBARD_RADIAL_MODE` log line matches the intended mode.
- **Executable.** Production runs use the frozen executable; diagnostic builds are separate. Record both SHA256s.
- **Processes.** Track solver PIDs, not only wrapper sessions. Kill and verify all owned processes at the end.

### 2.7 Reporting

- One markdown page per part: result tables, GATE verdicts, and at most 10 lines of interpretation.
- All hashes, commands and paths go in one JSON sidecar per part.
- Close with a 5-line summary: what passed, what was recorded, what is provisional, what stopped, and what comes next.

---

## 3. Part I — Correct KNOWN_ISSUES (documentation only; do first)

Replace the current Fe and Co entries in `tests/KNOWN_ISSUES.md` with these four. Rewrite the wording as needed, but keep the content.

1. **[OPEN] hcp Co impurity setup — bounds violation before SCF.**
   - Symptom: `lattice_cluster.f90:996`, dimension mismatch 3 vs 50000, in the Debug build.
   - The clean-baseline run exits SIGKILL before SCF.
   - It occurs at U=0, so it is not LDA+U.
   - Reference the deck and the C1 triage records. Part II/R1 resolves it.
2. **[OPEN] Co bulk U=0 non-finite.**
   - The RB6 records go non-finite after about 19–40 iterations, and this reproduces on `8b8f7ef`.
   - Reference the RB6 records.
3. **[OPEN, low priority] Fe U>0 embedding offset.**
   - An Fe impurity in an Fe host with U=3.8 eV, J=1.1 eV differs from RS bulk by −6.6×10⁻³ μB (about 735× the U=0 control).
   - Hubbard operators, Fermi level and Green's-function evaluation are identical; the recursion coefficients differ.
   - Probable cause: the cluster recipe mismatch (rc=80 for bulk, rc=100 for the impurity).
   - The matched-cluster test is blocked by entry 4. Ni does not show the effect.
4. **[OPEN] RS rc=100 runs non-finite or non-converging.**
   - The Fe rc=100 impurity-style deck goes non-finite after 102 iterations (U>0) and 34 (U=0).
   - The rc=100 U=0 bulk is not accepted after 600 steps.
   - Reproduces on `8b8f7ef`.
   - The rc=140 runs are INVALID_SETUP and are not part of this entry.

Commit.

---

## 4. Part II — RS robustness triage (timebox: 12 SCF runs plus offline analysis)

Impurity and surface runs use large clusters, so these failures sit on the critical path for Part III.

### R1. hcp Co bounds violation

1. Reproduce with the Debug build and capture the full backtrace.
2. Determine whether it is:
   - **an input error** in the deck (for example an inconsistent `ntype` / `nclu` / `inclu` / lattice entry); or
   - **a sizing bug** in `lattice_cluster.f90`.
3. Act on the answer:
   - **Code bug:** fix it with the smallest change (array sizing or index logic only; no physics). Add a regression test that builds this cluster.
   - **Input error:** add input validation with an actionable error message at read time, and a test.
4. **GATE:** the hcp Co impurity deck either runs through cluster setup or stops with the new actionable message.

Then run the U=0 Co impurity control once, as a RECORD.

### R2. Non-finite RS runs (Fe rc=100, Co bulk U=0)

1. **Locate the first non-finite quantity.** Instrument a diagnostic build (keep it separate from production) to report, per iteration, the first non-finite or out-of-range value among:
   - recursion coefficients A_n and B²_n, including any negative B²;
   - the terminator;
   - Q moments;
   - PL / `calculate_pl` inputs and outputs;
   - mixed q;
   - the radial solver.

   Start each run from its last finite state.
2. **Change one control at a time, from the same last finite state:**
   - `lld` 21 → 31;
   - the FP64 recursion path instead of the FP32 fast path, if a switch exists or can be added in a diagnostic build;
   - the energy window / `channels_ldos`;
   - beta halved.
3. **Decision:**

| Finding | Action |
|---|---|
| A specific numerical defect (for example negative B² from precision loss, an invalid `calculate_pl` branch, or a division by a vanishing denominator) | Fix it with the smallest change. Add a guard that turns silent NaN propagation into an actionable error or a safe fallback, plus a regression test. **GATE:** the Fe rc=100 U=0 bulk, restarted from the converged rc=80 endpoint, becomes RS-accepted, and the no-U standard fixtures stay byte-identical. |
| Only a parameter regime is unsafe | Document the safe ranges (rc, lld, window) in `lda_u_workflow.rst` and the input documentation. Add a log warning when a run enters the unsafe regime. |
| Neither | Record the localisation in KNOWN_ISSUES entry 4 and continue. |

4. **If R2 produces a fix and an RS-accepted rc=100 Fe bulk at U=3.8 eV, J=1.1 eV:**
   - compare the existing Fe U>0 impurity records against it;
   - report the ratio to the U=0 control (RECORD);
   - update KNOWN_ISSUES entry 3 accordingly. If the ratio is below 10, close entry 3 with "use the same rc for the bulk reference and the impurity cluster".

Commit after R1 and after R2.

---

## 5. Part III — Stage RD: impurities (`calctype='I'`), including KS-generated hosts

Choice of systems:

- **Ni in an Ni host** (fcc, U=3 eV, J=0.7 eV): the clean reference.
- **Fe in a bcc Fe host** (U=3.8 eV, J=1.1 eV): records only, because of the known offset.
- **The `B2FeCo_nsp1` fixture** (Fe impurity in B2 FeCo): the realistic heterogeneous case.

Reuse the decks from the C1/C2 provenance sidecars. Use the same rc for host and impurity wherever Part II established that it runs.

### RD1. Local moments for recursion atoms

- Apply `correct_d_shell_moments` per recursion atom in impurity mode, with:
  - V_in = the Hubbard d block actually inserted for that atom into the local cluster Hamiltonian (`hall`) in this iteration;
  - P and M1 from that atom's Green's function.
- Host atoms are frozen and never corrected; their mode is enforced by the provenance guard.
- `'local'` must be set explicitly for `calctype='I'` until RD5.
- **Unit test (GATE):** extend the two-site finite-H oracle so that one site is frozen and one is recursive. The corrected moments of the recursive site must match the direct eigenvector construction to 10⁻¹². A site-swap negative control must fail.

### RD2. Frozen replay on a converged local impurity state (Ni in Ni)

**GATE:**

- carrier ratio per spin within 1 ± 0.01;
- q2′ ≥ 0;
- the scalar-shift control is within its fitted allowance (as in RB3);
- the sign and spin-swap controls deviate from the true carrier by more than 10× the carrier error.

### RD3. Embedding equivalence (local mode throughout)

For each comparison, report:

- impurity-minus-bulk d occupations per spin and sphere moment;
- the matching U=0 control difference;
- the ratio between them;
- the maximum d-LDOS difference.

| ID | Comparison | Type and criterion |
|---|---|---|
| (a) | **RS-bulk host → same-material impurity.** Ni. | **GATE:** ratio to the U=0 control < 10. |
| (a′) | Same as (a) for Fe. | RECORD. |
| (b) | **KS-bulk host → RS impurity** (the supported mixed workflow). Ni: KS bulk (local, U=3 eV, J=0.7 eV, converged), then an RS impurity of the same material and U. The comparison reference is the KS bulk. | **GATE:** \|Δ(U>0) − Δ(U=0)\| ≤ \|δ_RB5(U)\| + 10 × \|U=0 RS control\|. Here Δ is the impurity-minus-KS-bulk difference, and δ_RB5(U) is the change in the RS–KS bulk difference between U=0 and that U, measured in RB5. Compute the U=0 RS control from (a). |
| (b′) | Same as (b) for Fe. | RECORD. |
| (c) | **Mode mismatch:** a KS-local host with an impurity forced to `'spectral'`. | **GATE:** the run stops with the provenance fatal error. |
| (d) | **Physics change:** local versus spectral for the Fe impurity in `B2FeCo_nsp1` (U on the impurity only) and for Ni in Ni. | RECORD: moment and d occupations. |

### RD4. Impurity energy (measure before adopting)

1. Determine from the code exactly which energy `calctype='I'` reports, and over which atoms. Write a one-paragraph description.
2. **W1 control:** a fixed external matrix λW on the impurity d shell only, at U=0, host frozen, λ ∈ {0, ±0.0025, ±0.005} Ry. Expectation if the embedded functional is stationary: dE_reported/dλ = Tr(W P_imp).
3. **Decision on the relative residual:**

| Residual | Action |
|---|---|
| ≤ 1 % | Run the impurity U-derivative test: vary only the impurity U (U0 = 3 eV for Ni, J = 0.7 eV fixed), at two step sizes, with candidates a, b and c. Allowance: A_U = 2 × Σ_s \|g_s · dq_s/dU\| + \|D(h1) − D(h2)\|, where g_s = (dE/dλ − Tr(W P))_s / (dq_s/dλ) comes from the W1 runs. **GATE**, on the average of the two step sizes: b within A_U, and c misses by more than 10 A_U. If it passes, route E_Hub(P_imp) into the impurity energy in local mode. |
| > 1 % | Document that the embedded-site energy is not variational. Adopt local moments for the density only, keep energy routing off for `calctype='I'`, and state this in the documentation. No correction formula may be invented. |

### RD5. Adoption

- **Who resolves to `local`.** Extend `auto` to resolve to `'local'` for eligible impurity runs:
  - `calctype='I'`;
  - fixed Liechtenstein U/J (`general` and/or `impurity`);
  - `nsp=1`, `lmax=2`;
  - no constraints;
  - every U-carrying host file stored in local mode.
- **Host files in spectral mode** (legacy) keep the run in spectral, with a log line explaining why.
- **Energy routing:** only if RD4 passed.
- **GATES:**
  - no-U byte identity;
  - explicit `'spectral'` impurity runs unchanged;
  - RS bulk and KS bulk behaviour unchanged;
  - one Ni impurity run under `auto` on an RS-local host, and one on a KS-local host, both resolve to `local` and reproduce the accepted RD3 endpoints.
- **Documentation:** update `lda_u_workflow.rst` and the input documentation with:
  - both host workflows (KS bulk → RS impurity, RS bulk → RS impurity);
  - the mode-consistency rule;
  - which energies are physical in which configuration.

  Also update `REMEDIATION_GATES.md` and `ENERGY_ACCOUNTING.md`.

Commit after each RD item.

---

## 6. Order of work

1. Part I.
2. Part II: R1 and R2 can run in parallel.
3. Part III: RD1 → RD2 → RD3 → RD4 → RD5. Within RD3, independent comparisons run in parallel.
4. Merge (section 7).

**Surfaces (Stage RE) are not part of this session.**

---

## 7. Merge (local only)

1. `git checkout fable_v4b_u && git pull --ff-only origin fable_v4b_u`
2. `git merge --no-ff fix/ldau-rs-impurity -m "Merge RS LDA+U: robustness triage and impurity stage (C5)"`
3. On the merged tree:
   - rebuild;
   - run the LDA+U unit tests;
   - run the no-U RS and KS byte-identity check against `8b8f7ef`;
   - run one RS bulk +U `auto` smoke run, confirming `HUBBARD_RADIAL_MODE mode=local ... local_energy_adopted=.true.`.
4. Do not push. Report the merge SHA for user review.
5. If the merge conflicts or a check fails, stop and report. Do not resolve conflicts in physics code on your own.

---

## 8. Final report

- **Part I:** the KNOWN_ISSUES diff.
- **Part II:** the R1 and R2 findings, any fixes with their regression tests, and the updated status of KNOWN_ISSUES entries 3 and 4.
- **Part III:** a table covering RD1–RD5, with the GATE verdicts, the RD3 comparison table and the RD4 energy decision.
- The merge SHA.
- The 5-line summary.

The next session after this one is Stage RE: surfaces, including the KS-host → RS-surface workflow and a 3d-adatom smoke test.

# C5.2 — RS LDA+U: resolve the no-U identity failure, integrate the Chebyshev hcp Co reference, finish Stage RD

**Session setup.** Sol is the orchestrator. Sol plans, delegates every execution task to Luna agents, reviews their output, and writes the final report. Sol does not run campaigns itself.

Luna agents are used as follows:

- **Workers** execute one well-defined task each. Independent runs may go in parallel.
- **One auditor** recomputes every GATE verdict and its arithmetic from the raw outputs before Sol records it, and checks each run's preflight (§2.6).
- **Workers never change criteria, tolerances or scope.** If a worker believes a criterion is wrong, it reports to Sol, and Sol applies §2.5.

This amendment is user-approved. It **supersedes C5.1** and extends C5. Record it verbatim in `docs/validation/fixed_ldau/rs/C5_AMENDMENT.md`, appended below the C5 text.

---

## 1. Context

### Repository and branches

- Repository: `rslmtoasa/rslmtoasa`, local clone at `/Users/andersb/Jobb/rslmto_fable_v4b`.
- Integration branch `fable_v4b_u` is at `908074f` (pushed).
- Work branch: `fix/ldau-rs-impurity`. It contains:
  - `d1b845a`: the C5 Parts I–II work, including the Co cluster-sizing fix (R1) in `source/lattice_cluster.f90` and the R2 diagnostics;
  - `96f1d72`: "Record validated Chebyshev RS hcp Co reference", committed by the user.
- **First action:** verify that `fix/ldau-rs-impurity` HEAD contains both commits and that the worktree is clean. Leave uncommitted generated outputs alone; do not commit them.
- **Never push.** At the end, merge locally into `fable_v4b_u` (§7) and report the merge SHA. The user reviews and pushes.
- `cleanup/` and `tests/gbt/` must never be committed.

### Read before planning

- `docs/validation/fixed_ldau/rs/C5_AMENDMENT.md`. Its §1 gives scope, physics and terminology. Its §5 gives the Stage RD specification (RD1–RD5), which this amendment executes with the changes in §5 below.
- `docs/validation/fixed_ldau/rs/C5_PART_II_RESULTS.md` and `C5_PART_I_assets/known_issues.diff`.
- **The hcp Co reference committed in `96f1d72`:** its validation summary and its four reference inputs. Note in particular the recursion method (`recur`), rc, lld, the energy window, `nsp`, the converged endpoint and the acceptance window.

### Scope (fixed)

**Primary target:** fixed-U LDA+U in real space (RS) for:

- bulk (`calctype='B'`)
- impurities (`calctype='I'`)
- surfaces (`calctype='S'`)

**Convenience path:** KS bulk. **A KS bulk host followed by an RS impurity or surface run is supported.** KS impurity is not.

**Guarded scope, do not touch:**

- SOC or noncollinear spin (`nsp>=2`)
- shells other than d, and `lmax/=2`
- HOH
- self-consistent U, intersite V, ACBN0

### Status

**Validated:**

- KS bulk and RS bulk local moments, with the energy routed as candidate b: E_local + E_Hub.
- In both, `auto` resolves to `local`.

**Not yet implemented:**

- Impurities and surfaces are still spectral, with no energy routing.

**What C5 left open:**

1. **RS no-U byte identity failed in nine scientific files** on `fix/ldau-rs-impurity`. That stopped C5.
2. **Part II, R2.** The Fe rc=100 block-recursion run puts more charge into a d channel than the shell can hold at iteration 90, then goes non-finite at iteration 102.
3. **Co controls** lacked a valid restart. The new Chebyshev hcp Co reference resolves this.

---

## 2. Standing rules (unchanged from C5; restated)

### 2.1 RS-accepted state

All of the following must hold:

- (a) outer Diff ≤ 3×10⁻⁸;
- (b) over the last 20 complete iterations, the spread (max − min) of the d occupation per spin is ≤ 5×10⁻⁷, and the spread of the sphere moment is ≤ 10⁻⁶ μB;
- (c) fixed-U P/V residuals are logged, not gated.

### 2.2 Fallback ladder for an SCF that is not accepted at its step cap

1. Continue from the last state for 200 steps with beta halved.
2. Then continue for 200 steps with `broyden` and beta = 0.01.
3. Then label the state PROVISIONAL, use it for RECORDS only, and move on.

Interrupted or retried attempts count toward the same budget.

### 2.3 RECORDS versus GATES

- RECORDS cannot fail and never block progress.
- Only items marked **GATE** can fail.

### 2.4 Crash triage

On a crash or a non-finite value:

1. Rerun once with the Debug, bounds-checked build and capture the file:line backtrace.
2. Rerun the same deck with the clean `8b8f7ef` binary.
3. Rerun with U=J=0.

If it also fails on `8b8f7ef` or at U=0, classify it as PRE-EXISTING and continue. STOP only if it is new relative to `8b8f7ef`, appears only with U>0, and lies in code changed since `8b8f7ef`.

### 2.5 STOP only if

- (a) the no-U criterion breaks. This is byte identity on the standard fixtures until §3 redefines it;
- (b) a GATE fails after one honest debugging attempt (state the hypothesis you tested);
- (c) §2.4 says stop;
- (d) continuing would require touching guarded scope or changing a physics formula.

**If a criterion turns out to rest on a wrong error model:** finish every other independent item, then report the exposing measurement together with a proposed corrected criterion. Never self-amend, and never loosen a tolerance after seeing results.

### 2.6 Preflight checklist (auditor verifies for every run)

- **Restart rule.** Start from the specified converged endpoint, never from fresh atomic potentials, unless a fresh start is explicitly requested.
- **Backend.** `use_kspace` sits in the correct namelist, confirmed from the live log.
- **Recursion method.** `recur` is confirmed from the live log.
- **Input parsing.** Every namelist is closed with `/` and the input parses, confirmed before counting steps.
- **Mixing.** The mixing parameters actually applied match the ladder step.
- **Resolved mode.** The `HUBBARD_RADIAL_MODE` log line matches the intended mode.
- **Executable.** Production and diagnostic builds are kept separate. Record both SHA256s.
- **Processes.** Track solver PIDs. Kill and verify all owned processes at the end.

### 2.7 Reporting

- One markdown page per part: tables, GATE verdicts, and at most 10 lines of interpretation.
- One JSON sidecar per part holds hashes, commands and paths.
- Close with a 5-line summary: what passed, what was recorded, what is provisional, what stopped, and what comes next.

---

## 3. Part A — Resolve the no-U identity failure (do first; timebox 8 runs)

**A1. Determinism control.**

- Run RS `bulk/bccFe` twice with the clean `8b8f7ef` binary, single thread (`OMP_NUM_THREADS=1` and the BLAS thread caps).
- Expectation: byte-identical.
- If they differ: STOP and report an environment problem.

**A2. Attribution.**

- Build each commit in `908074f..HEAD` that touches `source/`:
  - the parent of the R1 fix;
  - the R1 fix;
  - each R2 commit that touches production source.
- `96f1d72` should touch only the reference inputs and documentation; verify that.
- Run the RS no-U identity check at each build. Report the first commit that breaks it.

**A3. Magnitude.** For each of the nine differing files, report the maximum absolute and relative difference in each scientific quantity:

- d/s/p occupations;
- sphere moments;
- E_F;
- ETOT components;
- DOS values.

**A4. Decision (fixed in advance):**

| Finding | Action |
|---|---|
| (i) The cause is R1 and the max relative difference is ≤ 1×10⁻¹² (summation-order level) | Accept it and record "R1 summation-order change". For the rest of this session, the RS no-U criterion becomes: scientific quantities equal within 1×10⁻¹² relative, and KS no-U byte-identical. Regenerate the RS no-U reference outputs from the new tip and record both hashes. |
| (ii) The cause is R1 and the difference is larger | Restrict the fix so the new sizing applies only where the old sizing violates bounds, keeping the old path bit-for-bit otherwise. Re-run A2 at the restricted fix; byte identity must be restored. Record the unrestricted fix's quantitative effect on bcc Fe as a RECORD, and add a KNOWN_ISSUES entry asking the user to decide whether the old sizing affected results. |
| (iii) The cause is R2 production code | Remove the leak (diagnostics must not alter production paths) and restore byte identity. |
| (iv) Anything else | Report and STOP. |

Commit after A4.

---

## 4. Part B — Close Part II and integrate the Co reference (timebox 6 runs)

**B1. Occupation-capacity guard (production).**

- If any channel's occupation q0 per spin exceeds 2l+1 by more than 1×10⁻⁶, or a mixed q or PL value is non-finite, stop with an actionable error. The message names the atom, channel, iteration and value.
- Never propagate NaNs silently.
- Add a regression test with synthetic over-capacity input.
- The no-U criterion (as fixed in A4) must still hold.

**B2. Recursion-method attribution.** These are RECORDS, one run each, from converged endpoints and never from fresh potentials:

- **(a) Fe rc=100, Chebyshev.** Fe U=0 RS bulk at rc=100 with `recur='chebyshev'`, using the same Chebyshev settings as the Co reference, restarted from the converged rc=80 endpoint. Does it converge, or does the B1 guard fire?
- **(b) Co, block recursion.** Co U=0 RS bulk with `recur='block'`, using the Co reference geometry and its converged Chebyshev endpoint as the restart. Does the non-finite behaviour return?

**B3. KNOWN_ISSUES update.**

- **Entry 1 (hcp Co bounds violation):** CLOSED by R1, with a reference to the regression test.
- **Entry 2 (Co bulk non-finite):** state the recursion-method attribution from B2(b). Point to the Chebyshev reference as the validated Co recipe.
- **Entry 4 (rc=100 non-finite):** record the localisation (d occupation exceeds shell capacity at iteration 90, then non-finite), the B2(a) result, and that the B1 guard now reports the condition. If Chebyshev converges at rc=100, state that the issue is specific to block recursion at large rc.
- **Entry 3 (Fe embedding offset):** unchanged until Part D.

Commit.

---

## 5. Part C — Chebyshev compatibility of fixed-U LDA+U (GATE before Co enters Stage RD)

The local-moment correction and the Hubbard occupation matrix were validated only with block recursion. The Co reference uses Chebyshev, so the Chebyshev path has to be checked first.

**C1. Code audit.**

- Establish from the code whether, under `recur='chebyshev'`, the fixed-U path fills the full d-block Green's function g0(i,j) used for the matrix spectral function and the occupation matrix P.
- Establish the same for the energy-weighted occupation matrix M1, which the local moments need.
- State whether any stochastic trace estimate enters P or M1, and how many random vectors are used. A stochastic P is not acceptable for the occupation matrix.
- **Decision:**
  - P and M1 are not filled, or are stochastic → Chebyshev is unsupported for fixed U. Add an actionable fatal guard (`recur='chebyshev'` with fixed U), drop Co U>0 from this session, and continue with Ni and Fe only.
  - P and M1 are filled deterministically → continue with C2.

**C2. Block versus Chebyshev, bcc Fe RS bulk (rc=80):**

- at U=0;
- at U=2 eV, J=0 in local mode.

Compare the d occupations per spin, the sphere moment, E_F, and E_Hub.

**GATE:** the U=2 block–Chebyshev difference does not exceed the U=0 block–Chebyshev difference by more than a factor of 10 for any quantity. Report the ratio.

**C3. Carrier check on the Chebyshev U=2 local state** (frozen local-versus-spectral replay, as in RB3).

**GATE:**

- carrier ratio per spin within 1 ± 0.01;
- q2′ ≥ 0.

If C2 and C3 pass, Co enters Stage RD under Chebyshev. Otherwise, apply the guard from C1 and drop Co U>0.

Commit.

---

## 6. Part D — Stage RD: impurities (`calctype='I'`)

Execute C5 §5 (RD1–RD5) **with the following system definitions**.

### Matched recipes

- The host bulk reference and the impurity use the same rc, lld, energy window and recursion method.
- **Ni:** fcc, block recursion, rc=80, U=3 eV, J=0.7 eV. Reuse the converged rc=80 Ni bulk as the host. Ni is the gating system.
- **Fe:** bcc, rc=80, U=3.8 eV, J=1.1 eV. Use block recursion, or Chebyshev at rc=100 only if B2(a) converged. Fe is records only.
- **Co:** hcp, the committed Chebyshev reference recipe, U=3 eV, J=0.7 eV. Only if Part C passed.
- **The `B2FeCo_nsp1` fixture:** records only.
- If the impurity recipe cannot use the host's rc, report which parameter forces it to differ, then use the smallest rc that converges and passes the B1 guard, for **both** host and impurity.

### Changes to RD3 (all in local mode)

- **(a) RS-bulk host → same-material impurity.**
  - Ni: **GATE**, ratio to the U=0 control < 10.
  - Fe: RECORD. If the Fe ratio is below 10 with matched clusters, close KNOWN_ISSUES entry 3 with "use the same rc for the bulk reference and the impurity cluster".
  - Co: **GATE** if Part C passed, else omitted. It needs its own matched U=0 control first, built from the Co reference.
- **(b) KS-bulk host → RS impurity.**
  - Ni: **GATE** with the C5 formula: |Δ(U>0) − Δ(U=0)| ≤ |δ_RB5(U)| + 10 × |U=0 RS control|.
  - Fe: RECORD.
  - Co: RECORD (KS hcp Co with two atoms per cell). If KS hcp Co is impractical to converge, record why and skip it.
- **(c) Mode-mismatch guard:** **GATE**, unchanged.
- **(d) Local versus spectral:** RECORD, with Co added if Part C passed.

### RD2, RD4, RD5

- **RD2 and RD4** use Ni.
- **RD5 adoption:** `auto` resolves to `local` for eligible impurity runs as in C5. If Part C failed, this includes the Chebyshev guard. Energy routing is adopted only if RD4 passed.

Commit after each RD item.

---

## 7. Merge (local only)

1. `git checkout fable_v4b_u && git pull --ff-only origin fable_v4b_u`
2. `git merge --no-ff fix/ldau-rs-impurity -m "Merge RS LDA+U: robustness, Chebyshev compatibility, impurity stage (C5.2)"`
3. On the merged tree:
   - rebuild;
   - run the LDA+U unit tests;
   - run the no-U criterion as fixed in A4 (RS and KS);
   - run one RS bulk +U `auto` smoke run, confirming `mode=local ... local_energy_adopted=.true.`;
   - run one Ni impurity `auto` smoke run, confirming the RD5 behaviour.
4. **Do not push.** Report the merge SHA.
5. If the merge conflicts or a check fails, stop and report. Do not resolve conflicts in physics code on your own.

---

## 8. Final report

- **Part A:** attribution, magnitudes, and the A4 decision.
- **Part B:** the guard, the B2 records, and the KNOWN_ISSUES diff.
- **Part C:** the Chebyshev audit, C2/C3 verdicts, and the Co status.
- **Part D:** a table covering RD1–RD5 with GATE verdicts, the RD3 comparison table, and the RD4 energy decision.
- The merge SHA.
- The 5-line summary.

The next session is Stage RE: surfaces, including the KS-host → RS-surface workflow and a 3d-adatom smoke test.

Protocol amendment C5.3 (user-approved; extends C5.2; record verbatim in
docs/validation/fixed_ldau/rs/C5_AMENDMENT.md). Same session setup: Sol orchestrates, Luna
workers plus the Luna auditor.

1. A4 RESOLUTION (adopts the proposed no-reproduction branch)
   A4(v): not reproduced under the controlled one-thread conditions. Clean baseline, R1 parent
   and R1 were all byte-identical in A1/A2.
   - Classify the earlier mismatch as thread-order nondeterminism (OpenMP/BLAS reduction order;
     magnitudes about 4e-15).
   - Part A is CLOSED. R1 is accepted with no restriction.
   - The no-U criterion stays EXACT BYTE IDENTITY, but is always evaluated single-threaded:
     OMP_NUM_THREADS=1 and the BLAS/NUMEXPR thread caps, in both builds being compared.
     Multi-threaded runs remain allowed for SCF campaigns and records.
   - Add a KNOWN_ISSUES entry, informational: "RS/KS results are not bitwise reproducible across
     thread counts (reduction order, ~1e-15 relative). Use single-thread runs for identity
     checks."
   - Note the same requirement in the testing documentation wherever byte-identity checks are
     described.

2. RULE FOR OUTCOMES A DECISION TABLE DOES NOT COVER (applies to all remaining items)
   If a measured outcome matches no row of a decision table:
   - do not stop and do not change any criterion;
   - take the action that leaves production behaviour and defaults unchanged (no adoption, no
     fix);
   - record the outcome and that choice in the stage report, and continue with the next
     independent item.
   STOP remains limited to C5.2 §2.5 (a)-(d).

3. CONTINUE C5.2 AS WRITTEN
   Part B (B1-B3) -> Part C (C1-C3) -> Part D (RD1-RD5) -> §7 local merge.
   Use the single-thread rule from item 1 for every identity check, including those in RD5 and
   in the merge verification.
   Do not push. Report the merge SHA.


## User clarification — Part B Co continuation ladder (2026-10-20)

Context: the six-run Part B limit conflicted with the prescribed two-run continuation ladder for the unconverged Co record; five runs had been used, including the active Fe run. The clarification preserves the cap and directs that Co remain unaccepted without adoption.

Question: “The auditor found a protocol conflict: C5.2 caps Part B at six solver runs, but §2.2 requires two continuation runs for the unconverged Co record; five runs are already used, including the active Fe run. How should I reconcile this? Fe and independent preparation will continue while awaiting your answer.”

Reply: “Keep the six-run cap; record Co as finite but unaccepted, omit its ladder and continue with no adoption.”

## C5.4

C5.4 (user-approved; record verbatim in C5_AMENDMENT.md). OBJECTIVE: decide impurity adoption,
then merge. Token rules:
- Sol writes a plan of at most 1 page, dispatches the task cards below to ONE Luna worker
  sequentially, collects outputs, and runs ONE batched Luna audit of the GATE verdicts at the end.
- Sol does no per-run reviews and does not re-read protocol documents beyond the files named here.
- Every run passes the existing preflight script first. If no script exists, T0 writes one
  (python, at most 80 lines: namelist closure, use_kspace, recur, restart source, fresh run dir,
  thread caps).
- Identity checks are single-threaded (C5.3).

T1 RD2 re-judgement (offline, frozen RD2 Ni state, no SCF):
   (a) Ni scalar scan v in {+/-0.0025, +/-0.005, +/-0.01, +/-0.02} Ry; fit r(v) exactly as A2.2.
       Scalar control passes if its deviation is within the Ni-fitted allowance at its v.
   (b) Compute the predicted separation of each control from the frozen matrices. If it is
       < 10x the carrier error, mark the control NULL (A3.1).
   (c) Run the sign control (c + t/q0) if it was not already run.
   RD2 GATE = carrier ratio per spin within 1 +/- 0.01, q2' >= 0, scalar within the Ni fit,
   and every ADMISSIBLE control separates by > 10x the carrier error.
   Criterion correction (user-approved): the allowance is fitted per system/state; null
   controls are not failures.

T2 Occupation guard margin (production change):
   - Warning at q0 > 2l+1 + 1e-6.
   - Fatal only at q0 > 2l+1 + 0.05, or on any non-finite value.
   - Report the peak q0 reached in the Ni KS runs.
   - Unit test for both thresholds. Single-thread no-U identity must hold.
   - Resume the Ni KS U=0 and U=3 hosts from their last states (ladder steps 1-2 at most).

T3 Provenance guard scope: if the RD3 U=0 "mapping fatal" is the provenance guard firing for
   atoms with U=J=0, restrict the guard to U-carrying atoms and add a test. Then run the Ni U=0
   impurity control and evaluate RD3(a) Ni (GATE: ratio < 10) and, if T2 produced KS hosts,
   RD3(b) Ni (GATE as in C5).

T4 RD4 with the DEBUG build: W1 on the Ni impurity d shell, U=0, host frozen,
   lambda in {0, +/-0.0025, +/-0.005}. Apply the C5 RD4 decision table.
   Also record ONE RelWithDebInfo run (-g, -fbacktrace) of the Release crash to capture its
   backtrace, as a KNOWN_ISSUES entry; no fix.

T5 If RD2 and RD3(a) pass: finish RD5 (explicit-spectral identity, the auto cases), then do the
   C5.2 §7 local merge, report the SHA, and do not push.
   If either fails: commit the records and stop with a 5-line summary.

Timebox: 8 SCF runs in total (T2+T3+T4), plus offline replays.
Uncovered outcomes: C5.3 item 2 (no change, record, continue).

### C5.4 T1 result (2026-10-21)

The amended RD2 gate **PASSes** on the frozen Ni endpoint, with zero SCF solver invocations. Per-spin carrier ratios are 0.9914299182 and 0.9914299103; q2′ is 0.1378160163 and 0.1378160027 Ry². The eight-point per-spin A2.2 scalar fit gives an allowance of 0.0027005575 Ry at v=+0.01 Ry, including D0=5.0565×10⁻⁶ Ry; measured deviation is 0.0026958204 Ry, a 4.7371×10⁻⁶ Ry pass margin. The sign control predicts 1.4357881 Ry separation and measures 1.3962656 Ry against a 0.0615241 Ry threshold, so it is ADMISSIBLE and passes. The spin-swap control predicts 3.33×10⁻¹⁶ Ry, below that threshold, so it is NULL; its 0.0061524 Ry measured separation is not a failure. The sign control was already run and reused from the byte-identical frozen-state replay. Full fit data, coefficients, preflight and hashes are in `C5_4_T1_RD2.json`; reproduction is `C5_4_T1_scalar_scan.py`.

### C5.4 T2 result (2026-10-21)

The warning/fatal guard is implemented as specified: warning strictly above capacity + 1×10⁻⁶, fatal strictly above capacity + 0.05 or on non-finite values; call sites check active orbital channels. Four focused tests pass, including exact threshold boundaries and non-finite input. Standard full-single-thread no-U identity passes for both RS and KS against the frozen RD1 outputs, with exact artifact-name sets and byte identity after excluding only `testrun.log` and the registered runtime `cpu_s` field in `str.out`.

The two Ni KS continuations used one solver invocation each, at the next half-beta ladder steps, from the captured `Ni-imp_scf.nml` states. U=0 reached q0=5.1605390187 at iteration 2 and stopped at the revised fatal threshold; its prior peak was 5.6050891342. U=3 exited 0 as converged at iteration 42 with outer Diff 1.7×10⁻⁹ and no guard warning; its prior peak was 5.3792343437. The trace configuration retained only the rolling last iteration, so the unchanged §2.1 20-iteration occupation/moment window cannot be recomputed from full raw per-loop qP data. U=3 is therefore **not accepted**; no Ni KS endpoint is claimed. The exact new U=3 peak q0 is unavailable; the absence of warning records bounds it at ≤5.000001. Total T2 count is 4/8, including both identity candidates. Full commands, hashes, runs and limitations are in `C5_4_T2.json`; the append-only ledger is `C5_4_RUN_LEDGER.jsonl`.

### C5.4 T3 result (2026-10-21)

The prior Ni U=0 fatal was the local d-potential extraction reached for an untreated recursion atom with U=J=0. The narrow fix skips zero-U/J fixed-impurity atoms before the legacy bulk type lookup; general-U, self-consistent-U and test-operator paths remain enabled. The focused local-moments unit test passes. With the existing per-iteration watcher, the matched Ni U=0 RS impurity control passed the unchanged §2.1 criteria at its earliest complete window, iterations 7–26: Diff 1.60×10⁻⁸, d-occupation spreads 2.5655×10⁻⁷ and 2.2919×10⁻⁷, and sphere-moment spread 1.8940×10⁻⁷. Its peak q0 over 76 captured iterations was 4.7015216512. RD3(a) impurity-minus-bulk component ratios to the U=0 control are 0.64463 (d spin 1), 0.59821 (d spin 2), and 0.80695 (sphere moment), so the <10 gate passes. RD3(b) remains **NOT_EVALUATED** because T2 produced no accepted KS host endpoints. The control solver was stopped after the accepted window; this deliberate stop is exit 143. Cumulative count is 5/8. Full evidence and endpoint hashes are in `C5_4_T3_RD3.json`.

### C5.4 T4 result (2026-10-21)

The full Debug W1 series could not fit the remaining run budget together with the required RelWithDebInfo replay; no Debug λ point was launched, the RD4 table is **NOT_APPLIED**, and energy routing remains off. The one fresh, preflight-passing RelWithDebInfo λ=0 replay reproduced SIGSEGV during the first recursion, after the l=1, spin=2 Hubbard LDM physicality diagnostics. The captured `-g -fbacktrace` output contains three address-only frames and no source-symbolized frame; no crash cause or fix is claimed. The raw trace file is incomplete and no energy or accepted endpoint exists. One solver attempt used slot 6; two slots remain. See `C5_4_T4_RD4.json` and the OPEN entry in `tests/KNOWN_ISSUES.md`.

## C5.5

C5.5 (user-approved; record verbatim in C5_AMENDMENT.md). OBJECTIVE: adopt impurity local
moments (density only) and merge.
Token rules as in C5.4: Sol plans in at most 1/2 page, ONE Luna worker sequentially, ONE batched
audit at the end. Preflight script before every run. Identity checks single-threaded.
Budget: 3 SCF runs.

T1 RS-host auto GATE: Ni impurity on the accepted RS-local Ni host, with
   hubbard_radial_moments left at its default.
   GATE: the log shows mode=local, and the run reproduces the accepted RD3(a) Ni local endpoint
   (d occupations and moment within the U=0 control difference).

T2 KS-host smoke GATE (replaces the RD5 KS-host endpoint check):
   - Use the accepted KS bcc Fe U=2 eV, J=0 local bulk state from the KS validation evidence as
     the host.
   - Run an RS Fe impurity (same U, J; same rc as that host's RS counterpart, if one exists) with
     auto mode.
   - GATE: resolves mode=local; no provenance or mode-mismatch error; RS-accepted within ladder
     steps 1-2.
   - RECORD: impurity-minus-KS-bulk d occupations and moment.

T3 Adopt RD5 (density only):
   - auto resolves to local for calctype='I' when every U-carrying host file is stored in local
     mode; otherwise spectral, with a log line stating why.
   - Energy routing stays OFF for impurities.
   - GATES (single-thread): no-U RS/KS identity; explicit-spectral impurity unchanged; RS and KS
     bulk unchanged.
   - Documentation: lda_u_workflow.rst (both host workflows, the mode rule, "impurity energies
     exclude the Hubbard term").
   - KNOWN_ISSUES, OPEN entries:
       (a) KS-host -> RS-impurity numerical equivalence not yet validated (no accepted Ni KS
           host);
       (b) impurity energy (RD4) unmeasured; the Release crash in the W path is recorded.

T4 Local merge per C5.2 §7. Report the SHA. Do not push.

If T1 or T2 fails: commit the records and stop with a 5-line summary.
Uncovered outcomes: C5.3 item 2.

### C5.9 result (2026-10-22)

T1(a) passes: the regenerated T2a run reports outer Diff `2e-10` and R05 physical P/V streak 3. The d-up, d-down, and physical sphere moment `sum_l(q_l(up)-q_l(down))` from the actual written output agree with the immutable accepted-seed snapshot within `6.1e-9`; T1(c) passes because both written potential files record `local`/`ks`. The snapshot is bound to accepted restart SHA-256 `bb647254…ead6b` by the committed C5.5/C5.6 records. The original accepted raw source directory is empty; no deletion is inferred.

T1(b) is **UNEVALUATED / UNCOVERED** under C5.3 item 2. The regenerated output reports `E_F=-0.0829496275 Ry`, but no accepted endpoint EF is present in the committed proof or old checkpoint fields, and the old raw endpoint directory is empty. The initial deck value `fermi=-0.082950 Ry` is not evidence of the accepted endpoint EF and is not used as the comparison. Therefore T2b, production adoption, and merge were not performed.

The required bounded RD4 lambda-zero Release rerun passed preflight and captured one complete stable iteration/state pair before the external SIGTERM (exit 143); no crash occurred through iteration 1. This is a diagnostic regression only; no energy measurement was made. C5.9 used 0/2 campaign runs and 1 separately budgeted short invocation. Full evidence is in `C5_9_RESULTS.md` and `C5_9_RESULTS.json`.

### C5.5 clarification — run budget (2026-10-21)

Question: “Does the three-run cap cover the T1/T2 campaign attempts, with identity, unit,
and §7 verification tracked separately?”

Reply: “Yes—three campaign runs; verification separate”.

Count T1/T2 campaigns, including failed starts and continuations, against the three-run cap.
Track required identity, unit, and §7 verification invocations separately; they do not consume
the campaign-run cap.

## C5.6

C5.6 (user-approved; record verbatim in C5_AMENDMENT.md). OBJECTIVE: finish T2-T4 of C5.5.
Token rules as in C5.4. Budget: the 2 remaining SCF runs.

T0 Disk: delete raw run directories under /private/tmp that belong to CLOSED items
   (RA4, RB, C3, C4, C5 Parts A-B, RD1-RD3) whose compact evidence is committed.
   Keep:
     - the accepted checkpoints that later items use as hosts or references (KS Fe U=2,
       Ni RS-local host, the Ni U=0 control);
     - the user's Co files.
   Report the freed space.

T2a Regenerate the KS Fe host:
   - KS bcc Fe U=2 eV, J=0, auto mode, restarted from the accepted KS Fe state with the current
     binary.
   - Continue until KS-accepted (expect a few steps; ladder steps 1-2 at most).
   - Confirm the written potential file carries mode=local and backend=ks.

T2b KS-host smoke GATE as in C5.5 T2, using the regenerated host.

T3 As in C5.5 T3, plus: for any U-carrying host file without provenance, log a WARNING stating
   that its mode is unknown, that it is treated as spectral, and that the user should restart the
   host once with the current version to record its mode. Add a unit test for the warning.
   Add the same note to lda_u_workflow.rst.

T4 Local merge per C5.2 §7. Report the SHA. Do not push.

If T2a or T2b fails: commit the records and stop with a 5-line summary.
Uncovered outcomes: C5.3 item 2.

## C5.7

C5.7 (user-approved; record verbatim in C5_AMENDMENT.md). OBJECTIVE: fix the provenance-guard
scope, then finish T2-T4 of C5.6.
Token rules as in C5.4. Budget: 3 SCF runs.

T1 Guard scope fix (production change, design correction):
   - Enforce the mode-consistency check ONLY for atoms that are frozen in the current run: host
     atoms in calctype 'I'/'S' whose potentials are not updated by the SCF.
   - For atoms the SCF updates (all atoms in a bulk run; recursion atoms in 'I'/'S'): if the
     stored mode is missing or different, log INFO "restart seed in <mode>, run uses <mode>; seed
     will be rewritten", and proceed. The written file carries the run's mode and backend.
   - Apply the legacy-file WARNING from C5.6 T3 only to frozen host atoms.
   - Unit/fixture tests:
       (a) a bulk +U restart from a legacy file (no provenance) runs and writes provenance;
       (b) a bulk restart from a spectral-stored file with auto -> local runs (INFO logged);
       (c) an impurity with a frozen local host and a forced spectral impurity -> fatal
           (existing gate);
       (d) an impurity with a legacy frozen host -> WARNING, run is spectral.
   - GATE: single-thread no-U RS/KS identity.

T2a/T2b, T3, T4 exactly as in C5.6, using the fixed binary.
   (T3's legacy warning is now scoped as described in T1.)

If T2a or T2b fails: commit the records and stop with a 5-line summary.
Uncovered outcomes: C5.3 item 2.

## C5.8

C5.8 (user-approved; record verbatim in C5_AMENDMENT.md). OBJECTIVE: find and fix the
Release-only pre-SCF crash, then finish T2a-T4.
Token rules as in C5.4.
Budget:
  - up to 12 SHORT reproduction invocations (stop after the first SCF iteration or at the crash);
    these do not count as SCF runs;
  - 3 SCF runs.

T1 Reproduce and classify (short invocations, the T2a deck):
   (a) Release build at the current tip: confirm the segfault, with OMP_NUM_THREADS=1.
   (b) RelWithDebInfo (-O2 -g -fbacktrace -fcheck=bounds,pointer -finit-real=snan
       -finit-integer=-999999): capture the backtrace or the first check failure.
   (c) Debug build: does it pass?
   (d) Release builds at 908074f (fable_v4b_u tip) and at 8b8f7ef, same deck.
       Whatever needs adapting (for example provenance fields unknown to old binaries): remove
       only those lines, and record what was removed.

T2 Bisect: if (d) passes at 908074f, build Release at each commit in 908074f..HEAD that touches
   source/ (binary search) and report the first crashing commit.

T3 Fix:
   - Smallest change that removes the cause (uninitialised variable, out-of-bounds access,
     unallocated or unassociated array, etc.). No physics change.
   - Regression: the T2a deck passes in Release.
   - Also rerun the RD4 lambda=0 deck once in Release and record whether that crash is gone.
   - GATE: single-thread no-U RS/KS identity; the guard fixtures (13) still pass.
   If (d) also crashes at 908074f/8b8f7ef: classify as PRE-EXISTING, record the backtrace in
   KNOWN_ISSUES, and continue T2a with the Debug build instead.

T4 Continue C5.6 T2a, T2b, T3, T4 with the fixed binary (or the Debug binary per T3).
   Local merge only; report the SHA; do not push.

If T2a or T2b then fails for any other reason: commit the records and stop with a 5-line summary.
Uncovered outcomes: C5.3 item 2.

C5.9 (user-approved; record verbatim in C5_AMENDMENT.md). OBJECTIVE: accept the regenerated
KS Fe host, then finish T2b-T4.
Token rules as in C5.4. Budget: 2 SCF runs.

T1 KS host acceptance. This defines "KS-accepted" for host regeneration and replaces the RS
   window rule for this case. Evaluate from the existing T2a run; no new run.
   (a) It reported production KS convergence: outer threshold met and R05 P/V streak of 3.
   (b) d occupations per spin and the sphere moment agree with the previously accepted KS Fe
       U=2 eV, J=0 state within 1e-6, and E_F within 1e-6 Ry.
   (c) The written potential file carries mode=local and backend=ks.
   If (a)-(c) hold, the host is accepted. If (b) fails, record the differences and stop with a
   5-line summary.

T2 C5.6 T2b: KS-host smoke GATE (RS Fe impurity on this host, auto mode).
   RS acceptance (20-iteration window) applies to the impurity run as before.

T3, T4 as in C5.6. Local merge only; report the SHA; do not push.
   Also rerun the RD4 lambda=0 deck once in Release (short invocation) and record whether its
   crash is gone. No RD4 measurement.

Uncovered outcomes: C5.3 item 2.

C5.10 (user-approved; record verbatim in C5_AMENDMENT.md).

T1 KS Fe host: ACCEPTED. C5.9 T1 (a) and (c) pass; (b) passes on the d occupations and the
   sphere moment. The E_F comparison is N/A (no stored reference value); record the new E_F as
   the reference for future use.

Rule for the rest of C5: a sub-criterion whose reference data does not exist is recorded as N/A.
It does not block, provided the gate's remaining quantitative sub-criteria pass.

Continue with C5.9 T2, T3, T4 (budget 2 SCF runs). Local merge only; report the SHA; do not push.
Uncovered outcomes: C5.3 item 2.

C5.11 (user-approved; record verbatim in C5_AMENDMENT.md). OBJECTIVE: close the impurity auto
workflow. Base: fable_v4b_u at b04bcba; work branch fix/ldau-impurity-auto.
Token rules: Sol plans in at most 1/2 page; ONE Luna worker sequentially; ONE batched audit at the
end. Preflight script before every run. Identity checks single-threaded.
Budget: 2 SCF runs.

T1 Apply docs/validation/fixed_ldau/rs/C5_10_AUTO_CANDIDATE.patch: auto -> local for
   calctype='I' also when hubbard_u_impurity is set, provided every U-carrying frozen host is
   stored local. Keep its unit test (UnitAutoImpurityLocalMode).
   Add a fixture test: a hubbard_u_impurity deck on a local U-carrying host resolves to local with
   no provenance error.
   GATES: single-thread no-U RS/KS identity; the 13 guard fixtures; the existing RD3(a) Ni
   auto case still resolves to local.

T2 KS-host smoke GATE, redone. Build the deck by COPYING the existing accepted Fe RS impurity
   deck (RA4/C3 rc=80 recipe, ct = 3.0 per type) and changing ONLY:
     - the host potential file -> the accepted KS Fe host from C5.8/C5.10;
     - U = 2 eV, J = 0 on host and impurity.
   Before launch, diff the deck against its source; the diff must show only those changes.
   GATE: resolves mode=local; no provenance error; RS-accepted within ladder steps 1-2.
   RECORD: impurity-minus-KS-bulk d occupations and moment.

T3 Documentation (lda_u_workflow.rst): both host workflows, both U-input styles, the mode rule,
   "impurity energies exclude the Hubbard term", and the legacy-file restart note.
   KNOWN_ISSUES: close the impurity-auto entry.

T4 Merge locally into fable_v4b_u with --no-ff. Report the SHA. Do not push.

If T2 fails: commit the records and stop with a 5-line summary.
Uncovered outcomes: C5.3 item 2.

C5.12 (user-approved; record verbatim in C5_AMENDMENT.md). OBJECTIVE: resolve the C5.11 T2
failure. Work branch fix/ldau-impurity-auto.
Token rules as in C5.11. Budget: 2 SCF runs.

T0 (no runs) Compare the T2 deck's &energy fermi value (and its energy window) with the KS host's
   endpoint E_F, -0.082949627451487762 Ry. Report both values.

T1a If T0 shows a mismatch:
   - Correct the T2 deck: fermi = the KS host E_F. Keep the window width, shifted by the same
     amount as E_F.
   - Diff against the C5.11 deck; it may differ ONLY in fermi and the window.
   - Rerun T2 (GATE and RECORD as in C5.11).

T1b If T0 shows no mismatch:
   - Run ONE U=0 control: same deck, same KS host, U=J=0.
       * converges  -> the problem is specific to +U on a KS host: RECORD and stop.
       * oscillates -> KS-host -> RS-impurity embedding is generally unstable for this recipe:
         KNOWN_ISSUES entry and stop.
   - Do not merge in either case.

T2 If T1a passes: complete C5.11 T3 (documentation) and T4 (local merge, report the SHA, do not
   push).

Deck rule for all future impurity decks: copy from a working deck, then set fermi and the energy
window to the HOST endpoint values, and diff before launch. Add this rule to the preflight script
as a check: fermi within 1e-6 Ry of the host E_F recorded in the host's run output.
Uncovered outcomes: C5.3 item 2.
