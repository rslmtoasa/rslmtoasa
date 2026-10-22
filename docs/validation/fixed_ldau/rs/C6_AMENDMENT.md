C6 (user-approved; record verbatim in docs/validation/fixed_ldau/rs/C6_AMENDMENT.md).
OBJECTIVE: LDA+U local moments for surfaces (calctype='S'), density only.
Base: fable_v4b_u at fe4ae42 (or later). Work branch: fix/ldau-surface.
Token rules: Sol plans in at most 1/2 page; ONE Luna worker sequentially; ONE batched audit at the
end. Preflight script (including the host-E_F check) before every run. Identity checks
single-threaded.
Budget: 4 SCF runs.
Standing rules: C5.2 §2 and C5.3 item 2. Smoke gates use the C5.13 smoke criteria (final Diff
<= 1e-7; 20-step spreads <= 1e-5); the strict RS window applies only to references.

T0 (no runs) Identify the existing surface fixture(s) under tests/scf/cases/surface/. Report:
   - how calctype='S' sets up host (frozen) and recursion layers;
   - whether E_F is fixed to the host value as in impurity runs;
   - the layer count;
   - which host potential files it reads.
   Pick ONE fixture for this stage (prefer Fe(001), or the simplest bcc/fcc surface present).
   Decks are built by COPYING that fixture and changing only: the host file, U/J, and
   fermi/window set to the host endpoint. Diff before launch.

T1 Code: extend auto -> local to calctype='S' under the same conditions as 'I' (every U-carrying
   frozen host stored local). Confirm is_frozen_hubbard_host classifies the surface host layers
   correctly. Energy routing stays OFF for 'S'.
   Add a unit/fixture test: a surface deck on a local U-carrying host resolves to local.
   GATES: single-thread no-U RS/KS identity; the 13 guard fixtures; the impurity auto fixture.

T2 RS-host surface, using the accepted RS-local bulk host of the same material, with U on all d
   sites. Two runs:
   (a) U=0;
   (b) U>0 in auto mode.
   GATE (smoke, on b): mode=local; no provenance error; converges under the smoke criteria.
   RECORD: layer-resolved d occupations and moments for (a) and (b); the deepest-layer minus
   bulk difference at U>0 relative to the same at U=0.

T3 KS-host surface smoke: the same deck with the accepted KS-local host of the same material
   (regenerate it once with the current binary if it lacks provenance; that uses 1 run).
   GATE: smoke criteria as in T2. RECORD: layer moments.

T4 Documentation: the surface workflow (both hosts), "surface energies exclude the Hubbard term",
   and the deck rule. KNOWN_ISSUES: new OPEN entries for anything recorded but unresolved.
   Local --no-ff merge into fable_v4b_u. Checks on the merged tree: identity, guard fixtures,
   the impurity and surface auto fixtures. Report the SHA. Do not push.

If T1 or a smoke GATE fails: commit the records and stop with a 5-line summary.
Uncovered outcomes: C5.3 item 2.

C6.1 (user-approved; record verbatim in C6_AMENDMENT.md). OBJECTIVE: complete Stage RE smoke
validation. Branch fix/ldau-surface (fa68c56c).
Token rules as in C6. Budget: 4 SCF runs. Smoke criteria from C5.13. Standing rules: C5.2 §2,
C5.3 item 2.

T1 (no SCF) Rerun the gates for the C6 T1 code change: single-thread no-U RS/KS identity, the 13
   guard fixtures, the impurity auto fixture, and the surface auto unit test. Any failure: stop.

T2 Ni(001), RS host. Build the deck by COPYING tests/scf/cases/surface/fccCu001 and changing only:
   - alat and wav, taken from the accepted Ni bulk deck;
   - labels and potential files: host = the accepted RS Ni bulk; surface/subsurface seeds =
     copies of that Ni bulk potential; ES stays ES;
   - U/J on the Ni d channels only;
   - fermi and window = the Ni host endpoint.
   First verify that layer positions are in alat units. If they are absolute, STOP and report;
   do not edit geometry. The worker lists every changed line with its reason; the auditor checks
   the diff.
   Runs:
   (a) U=0 with the accepted Ni U=0 bulk as host;
   (b) U=3/J=0.7 auto with the accepted RS-local Ni U=3/J=0.7 host.
   GATE (smoke, b): mode=local; no provenance error; smoke convergence.
   RECORD: layer-resolved d occupations and moments for (a) and (b), and surface minus bulk at
   U>0 versus U=0.

T3 Cu(001), KS host (workflow smoke):
   (a) KS fcc Cu bulk, U=3/J=0, auto (writes local provenance); KS production convergence.
   (b) Unchanged fccCu001 deck with that host, the same U, and fermi/window = the KS host
       endpoint.
   GATE (smoke, b): as in T2. RECORD: layer d occupations.

T4 Documentation and KNOWN_ISSUES (close the C6 OPEN entry, or update it). Local --no-ff merge
   into fable_v4b_u, with the merged-tree checks as in C6 T4. Report the SHA. Do not push.

If a GATE fails: commit the records and stop with a 5-line summary.
Uncovered outcomes: C5.3 item 2.

C6.2 (user-approved; record verbatim in C6_AMENDMENT.md). Branch fix/ldau-surface (945ff47).
Token rules and budget as in C6.1 (4 SCF runs).

T1 MANDATORY, before anything else: rerun, on the current branch tip, the exact C5.13 gate
   scripts:
   - single-thread RS identity (27 cases) and KS identity (24 cases);
   - the 13 guard fixtures;
   - the impurity auto fixture;
   - the surface auto unit test.
   Report counts. Prior evidence does NOT substitute for this. Any failure: stop.

T2/T3 as in C6.1, with one additional allowed deck change in every copied deck: nsp = 1.
   - Ni: the surface/subsurface seeds are copies of the accepted nsp=1 Ni bulk potential, as
     already planned.
   - Cu: if the fixture's ES.nml / Cu-S.nml seeds fail to read under nsp=1, replace them with
     copies of the converged Cu bulk host potential (ES stays an empty sphere) and record the
     substitution.
   - No other deck changes beyond the C6.1 list.

T4 As in C6.1.

Uncovered outcomes: C5.3 item 2. Do not stop at an uncovered outcome if independent items remain.

C6.4 (user-approved; record verbatim in C6_AMENDMENT.md). Branch fix/ldau-surface (1ebee1e).
Token rules as in C6. Python: .venv. Budget: 4 new SCF runs.

T1 Code fix: in bands.f90 local-moment correction, use the recursion-site -> cluster-site -> type
   map (get_rd1_recursion_site_type) and the inserted-V_in path for calctype 'S' exactly as for
   'I'. Generalise rd1_local_impurity_eligible() to embedded runs ('I' and 'S') under the same
   conditions; keep the bulk path unchanged. Atoms whose type carries no U/J (e.g. ES) are skipped
   before extraction, as in the impurity path.
   Unit test (extend the frozen-plus-recursive-site oracle): a surface-like layout with
     - one frozen host type with U,
     - one empty-sphere recursion type without U,
     - one recursion type with U;
   the mapping returns the correct types, the U site's corrected moments match the direct
   construction to 1e-12, and the ES site is skipped.
   GATES (same scripts as C6.2 T1): RS 27/27, KS 24/24, guards 13/13, impurity and surface auto
   fixtures, plus the new test.

T2/T3 as in C6.2, with additional allowed deck changes:
   - nstep = 200;
   - &mix parameters copied from the accepted bulk host deck of the same material.
   Rerun T2(a) and T2(b), then T3(a) and T3(b). Smoke criteria from C5.13.

T4 As in C6.1. Uncovered outcomes: C5.3 item 2.

C6.5 (user-approved; record verbatim in C6_AMENDMENT.md). OBJECTIVE: enforce preflight
mechanically. No SCF runs.

T1 Add tools/run_rslmto_checked.py (repository .venv):
   - it runs the existing preflight on the target run directory;
   - it writes preflight.json;
   - it launches rslmto.x ONLY if preflight passes, otherwise it exits non-zero without
     launching.
   All future campaign runs must be launched through this wrapper. The auditor rejects any run
   without a matching preflight.json.

T2 Unit test: a deck with an unclosed namelist and a deck with a fermi/host-E_F mismatch are both
   refused; a valid deck launches (dry-run mode).

T3 Mark the surface stage PAUSED in C6_AMENDMENT.md, pending a user-provided converged no-U
   Ni(001) or Fe(001) reference (nsp=1). Commit. No merge.

Surface stage status: PAUSED pending a user-provided converged no-U Ni(001) or Fe(001) reference
(nsp=1).
