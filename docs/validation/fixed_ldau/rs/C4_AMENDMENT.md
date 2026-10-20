# C4 — RS LDA+U: finish the bulk energy stage, close the Fe embedding question, merge

**Session setup.** Sol is the orchestrator. Sol plans, delegates every execution task to Luna agents, reviews their output, and writes the final report. Sol does not run campaigns itself. One Luna agent acts as an **independent auditor**. Before any GATE verdict is recorded, the auditor recomputes the verdict and its arithmetic from the raw outputs.

This amendment is user-approved. Record it verbatim in `docs/validation/fixed_ldau/rs/C4_AMENDMENT.md`. It is self-contained: the earlier amendments C1–C3 are summarised in section 2 where they still apply.

---

## 1. Context

### Repository and branches

- Repository: `rslmtoasa/rslmtoasa`, local clone at `/Users/andersb/Jobb/rslmto_fable_v4b`.
- Integration branch: `fable_v4b_u`. Work branch: `fix/ldau-rs-embedding`.
- Never push the work branch. At the end, merge it into `fable_v4b_u` (section 6).
- `cleanup/` holds user-owned protocol files. Read them, never edit them, never commit them.
- `tests/gbt/` belongs to a different effort. Do not add it to this branch.
- First action: check `git status` on `fix/ldau-rs-embedding`. If anything is uncommitted, commit it in logical pieces before starting:
  1. source and unit tests;
  2. fixtures and validation scripts;
  3. `docs/validation/fixed_ldau/rs/`.

### Scientific scope (fixed)

**Primary target:** fixed-U LDA+U in real space (RS recursion) for:

- bulk (`calctype='B'`)
- impurities (`calctype='I'`)
- surfaces (`calctype='S'`)

**Convenience path:** KS bulk. A KS bulk host followed by an RS impurity or surface run is a supported workflow. KS impurity is not.

**Guarded scope, do not touch:**

- SOC or noncollinear spin (`nsp>=2`)
- shells other than d, and `lmax/=2`
- HOH
- self-consistent U, intersite V, ACBN0

### Physics and terminology

**Local radial moments.** With `hubbard_radial_moments='local'`, the radial step receives the band centre and second moment of the LDA part of the Hamiltonian:

- c′ = c − Tr(V_in P)/q0
- q2′ = q2 − 2 Tr(V_in M1ᶜ) + Tr(V_in² P) − t²/q0

This is implemented in `math.f90::correct_d_shell_moments`. In `'spectral'` mode the Hubbard shift leaks into the radial functions instead.

**Energy candidates.**

| Candidate | Definition |
|---|---|
| a | E_spectral − Tr(V_in P) + E_Hub |
| b | E_local + E_Hub |
| c | E_spectral + E_Hub |

The Hellmann–Feynman (HF) theorem is the arbiter. On a converged state of the correct functional, dE/dX = ∂E_Hub/∂X at fixed occupation matrix P, for X = U or J.

**KS bulk status.** Done. Candidate b is adopted, and `auto` resolves to `local` for clean-bulk KS runs.

**Base ASA non-stationarity.** At U=0, a fixed external d-level shift λW gives dE/dλ ≠ Tr(WP) by:

- −0.14 % in KS;
- −0.16 to −0.17 % in RS.

It is a property of the ASA energy assembly, not of either backend, and it sets the precision floor of every HF test.

### Results so far on this branch

**RA1–RA3: done.** Impurity diagnostic indexing, moment-mode provenance in potential files with a mismatch guard, `nsp=1` fixtures with the `nsp=2` decks kept as guard tests, and VAL-04 updated.

**RA4 (records, spectral mode).**

- Fe, U=0: the impurity reproduces bulk to 2×10⁻⁶ (d up), 7.5×10⁻⁶ (d down) and 9.0×10⁻⁶ μB (moment).
- Fe, U=3.8 eV, J=1.1 eV: the impurity differs from bulk by **−6.6×10⁻³ μB**, about 735× the U=0 control, and the same from both seeds. V. Frilén's 2025 report shows the same sign and size.
- Ni, U=3 eV, J=0.7 eV: impurity-vs-bulk discrepancies are *smaller* than its U=0 control. Ni embeds cleanly.
- hcp Co: the impurity setup fails with non-finite values at U=0, and this reproduces on clean `8b8f7ef`. It is a pre-existing issue, and Co has been replaced by Ni.

**RB0–RB4 gates: passed.**

- RS-specific R05 stopping tolerances: P ≤ 5×10⁻⁷, V ≤ 10⁻⁷ Ry, streak 3. KS is unchanged.
- No-U identity holds.
- Three local Fe states are accepted: U=2 eV with J=0, U=4 eV, and U=2 eV with J=0.9 eV.
- Carrier ratios: 0.993 and 0.9998.
- Restart passes; serial and MPI runs agree exactly.

**RB5 and RB6 (records).** Done. One deviation is recorded: a U=0 continuation used unchanged beta.

**C3 Part B (RS energy).**

- **RC0.** RS SCF uses step-function occupations integrated to E_F, so the functional is E, not E − TS.
- **RC1.**
  - h=0.005 Ry, accepted endpoints: residual −0.171 %.
  - h=0.0025 Ry, one endpoint provisional: −0.159 %.
- **RC2, U=2 eV, J=0.** Candidate b is within the C3 allowance (3.205×10⁻⁴ Ry/eV) at both step sizes, and c misses by more than 10× that allowance.
- **RC2, J=0.5 eV.**
  - C3 allowance: 3.409×10⁻⁴ Ry/eV.
  - Candidate b error: 3.413×10⁻⁴ at h=0.05 and 2.988×10⁻⁴ at h=0.10.
  - Candidate c misses by 0.117 Ry/eV.
  - The run stopped on a 4.45×10⁻⁷ excess. That allowance formula had the wrong error model; it is replaced in section 3.
- **RC2, J=0.9 eV.** One accepted local endpoint exists; the rest is incomplete.
- **RC3.** Not run.

**C3 Part A (Fe embedding localisation).**

- **A1.** The Fermi level is correct (U>0 −0.163441 Ry, U=0 −0.069252 Ry).
- **A2.** All 15 cluster on-site blocks, including the Hubbard blocks, and the first- and second-shell hoppings are identical to bulk.
- **A3, partial.**
  - The d-subblock A_n (n=0–19) and B²_dd (n=0–20) show no level above 10× the U=0 control.
  - The terminator differences are smaller at U>0 than at U=0.
  - Full blocks and exact B_dd were not compared.
  - The production recursion has an FP32 fast path.
- The C2 frozen comparison showed that the first-iteration Green's function already differs between impurity and bulk.
- **Cluster mismatch.** The bulk reference uses `rc=80`; the impurity cluster uses `rc=100`.
- **The `rc=140` runs are INVALID_SETUP.** Both started from fresh atomic potentials, not the converged endpoint. The bulk run diverged (last complete Diff 0.59) and went non-finite at step 184.

### Key documents

- `docs/validation/fixed_ldau/rs/C1_RA4_RESULTS.md`, `C1_RB_RESULTS.md`, `C3_RC_RESULTS.md`, `C3_PART_A_RESULTS.md`, plus their `*_PROVENANCE.json` sidecars.
- `cleanup/RS_LMTO_ASA_LDAU_RS_Embedding_Protocol.md`, the original stage plan.
- `docs/validation/fixed_ldau/ENERGY_ACCOUNTING.md`, `LOCAL_MOMENTS.md`, `REMEDIATION_GATES.md`.

---

## 2. Standing rules (carried over from C1–C3; these do not change)

### 2.1 RS-accepted state

All of the following must hold:

- (a) outer Diff ≤ 3×10⁻⁸;
- (b) over the last 20 complete iterations, the spread (max − min) of the d occupation per spin is ≤ 5×10⁻⁷, and the spread of the sphere moment is ≤ 10⁻⁶ μB;
- (c) fixed-U P/V residuals are logged, not gated.

### 2.2 Fallback ladder for an SCF that is not accepted at its step cap

1. Continue from the last state for 200 steps with beta halved.
2. Then continue for 200 steps with `broyden` and beta = 0.01.
3. Then label the state PROVISIONAL, use it for RECORDS only, and move on.

Steps from interrupted or retried attempts count toward the same budget.

### 2.3 RECORDS versus GATES

- RECORDS are measured and tabulated. They cannot fail, and they never block progress.
- Only items marked **GATE** can fail.

### 2.4 Crash triage (C2)

On a crash or a non-finite value:

1. Rerun once with the Debug, bounds-checked build and capture the file:line backtrace.
2. Rerun the same deck with the clean `8b8f7ef` binary.
3. Rerun with U=J=0.

If it also fails on `8b8f7ef` or at U=0, record it as a PRE-EXISTING ISSUE and continue. STOP only if it is new relative to `8b8f7ef`, appears only with U>0, and lies in code this branch changed.

### 2.5 STOP only if

- (a) no-U byte identity breaks;
- (b) a GATE fails after one honest debugging attempt (state the hypothesis you tested);
- (c) section 2.4 says stop;
- (d) continuing would require touching guarded scope or changing a physics formula.

Everything else: write it down and keep going.

### 2.6 Preflight checklist for every run (lessons from C1–C3)

- **Restart rule.** Start from the specified converged endpoint, never from fresh atomic potentials, unless a fresh start is explicitly requested.
- **Backend.** Confirm it from the live log (`use_kspace` inside the correct namelist).
- **Input parsing.** Confirm the input parses (every namelist closed with `/`) before counting any steps.
- **Mixing.** Confirm the mixing parameters actually applied match the ladder step.
- **Executable.** Production runs use the frozen executable, and diagnostic builds stay separate. Record both SHA256s.
- **Processes.** Track solver PIDs, not only wrapper sessions. Kill and verify all owned processes at the end.

### 2.7 Reporting

- One markdown page per item: result tables, GATE verdicts, and at most 10 lines of interpretation.
- All hashes, commands and paths go in one JSON sidecar per stage. Do not repeat them in prose.
- Close with a 5-line summary: what passed, what was recorded, what is provisional, what stopped, and what comes next.

---

## 3. Part B — RS bulk energy (execute first)

### B1. RC1 verdict (decided)

RC1 = **PASS** on the accepted h=0.005 pair (−0.171 %, gate 1 %). Keep h=0.0025 as a supporting record. No reruns.

### B2. RC2 allowance and gate

This replaces the C3 formula, which scaled the error with ∂E_Hub/∂X instead of with the occupation response.

For X ∈ {U, J}:

```
A_X = 2 × Σ_s |g_s · dq_s/dX|  +  |D_X(h1) − D_X(h2)|
```

where:

- **g_s** (Ry per electron) is the RS gradient mismatch per spin from the RC1 W1 runs:
  `g_s = (dE/dλ − Tr(W P))_s / (dq_s/dλ)`.
  If the spin split is not available, use the total residual with dq_tot/dλ.
- **dq_s/dX** is taken from the same converged series used for the derivative.
- **D_X(h)** is the central-difference derivative at step h.

The decision this protects is identifying the correct energy candidate. The b-versus-c separation is about 0.1 Ry/eV.

**GATE**, evaluated on the **average of the two step-size estimates**:

- candidate b is within A_X;
- candidate c misses by more than 10 A_X.

Candidate a is a RECORD.

**Execution:**

1. Re-evaluate U (U0=2 eV, J=0) and J=0.5 eV from the existing runs. No reruns.
2. Complete J=0.9 eV (U=2 eV) under the same rule, resuming from its accepted endpoint and following section 2.2. Use the same step sizes as J=0.5.

The auditor recomputes g_s, dq_s/dX, A_X and every candidate error from the raw outputs before the verdict is recorded.

### B3. RC3 adoption (only if B2 passes for U, J=0.5 and J=0.9)

Extend the eligibility logic (`is_eligible_local_hubbard_ks` or an RS twin) so that `auto` resolves to `local` and E_Hub is routed into `physical_total_energy` for clean RS bulk:

- `calctype='B'`, `nbulk=0`, `nrec=ntype`
- `nsp=1`, `lmax=2`
- fixed Liechtenstein U/J, no constraints

Impurity and surface runs stay spectral, with no energy routing.

**GATE:**

- no-U output byte-identical to `8b8f7ef` (RS and KS);
- explicit `'spectral'` U>0 unchanged;
- KS behaviour unchanged;
- one RS bulk U=2 eV run under `auto` reproduces the accepted local endpoint and reports candidate b as its physical total energy.

Update `ENERGY_ACCOUNTING.md`, `REMEDIATION_GATES.md` (R04 scope: KS and RS bulk) and `docs/source/theory/lda_u_workflow.rst`. Commit.

---

## 4. Part A — Fe embedding discrepancy (after Part B)

**Timebox:** 3 SCF runs plus offline replays. No production change unless A4.3 finds an evaluation inconsistency.

### A4.1 Matched-cluster bulk

- Run RS bulk bcc Fe at U=3.8 eV, J=1.1 eV with `rc=100`, the same cluster recipe as the impurity.
- **Restart from the converged `rc=80` +U bulk endpoint.**
- Accept per section 2.1.
- Compare the existing accepted Fe U>0 impurity records against this bulk, and report the ratio to the U=0 control.
- If a matched U=0 reference is needed, run the `rc=100` U=0 bulk with the same restart rule.

### A4.2 Swap test (offline, no SCF)

Use the frozen iteration-1 captures, in production precision.

1. Evaluate the bulk run's recursion coefficients plus terminator with the **impurity-path** Green's-function routine.
2. Evaluate the impurity coefficients plus terminator with the **bulk-path** routine.
3. Compare the d-LDOS and P_out of each against the native results.

Interpretation:

- **The swap reproduces the other path** → the coefficients differ. Compare the **full** blocks A_n and B_n (all orbitals, all levels; an FP64 replay is allowed as a diagnostic) and report the first differing level and block.
- **The swap does not reproduce** → the evaluation differs. Change one input at a time (terminator, energy mesh or window, precision path) and report which one accounts for the difference.

### A4.3 Decision

| Outcome | Action |
|---|---|
| A4.1 brings the ratio below 10 | Document "use the same rc for the bulk reference and the impurity cluster in +U comparisons" in `lda_u_workflow.rst`. CLOSE. |
| A4.2 finds an evaluation inconsistency between the paths | Fix it with the smallest change. Add a regression test (same coefficients through both paths give the same Green's function). Rerun the two Fe U>0 impurity records. CLOSE. |
| Neither | Add a KNOWN_ISSUES entry with the A4.1 and A4.2 results. CLOSE. |

### A4.4 Documentation

- **`KNOWN_ISSUES.md`, entry "RS impurity numerical instability".** Covers only the hcp Co non-finite failure reproduced on `8b8f7ef`, with the deck and the triage results.
- **`rc=140` runs.** Reclassify them as INVALID_SETUP in the Part A record. They are not an instability.
- **Base non-stationarity.** Note in the docs: 0.14 % in KS and 0.16–0.17 % in RS, a property of the ASA energy assembly.

Commit.

---

## 5. Order of work

1. Clean up the working tree (section 1).
2. Part B: B1 → B2 → B3.
3. Part A: A4.1 and A4.2 in parallel → A4.3 → A4.4.
4. Merge (section 6).

Commit after each completed item, with a descriptive message, on `fix/ldau-rs-embedding`.

---

## 6. Merge

1. `git checkout fable_v4b_u && git pull --ff-only origin fable_v4b_u`
2. `git merge --no-ff fix/ldau-rs-embedding -m "Merge RS LDA+U: bulk energy (RC), Fe embedding (C4 Part A)"`
3. On the merged tree:
   - rebuild;
   - run the LDA+U unit tests;
   - run the no-U RS and KS bcc Fe byte-identity check against `8b8f7ef`.
4. Push `fable_v4b_u` only if all of those pass.
5. If the merge conflicts or a check fails, stop and report. Do not resolve conflicts in physics code on your own.

Never push `fix/ldau-rs-embedding`.

---

## 7. End of session

Stop after the merge, or at a failed GATE.

Final report:

- the Part B table: RC1, U, J=0.5 and J=0.9, each with A_X, the b and c errors and the verdict, and RC3;
- the Part A result and its decision;
- the KNOWN_ISSUES entries added;
- the merge commit SHA;
- the 5-line summary.

The next stage after this session is RD/RE: local moments and energy for impurities and surfaces, including the KS-host → RS-impurity workflow. Do not start it.