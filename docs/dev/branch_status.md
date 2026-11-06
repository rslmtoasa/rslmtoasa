# Branch status of rslmtoasa/rslmtoasa and main / fable development inventory

Written 2026-11-06 from a read-only fetch of all remote branches. Branch tips change; re-run before acting on it.

Snapshot: 2026-11-06, all 54 remote branches fetched (read-only). Nothing was pushed, deleted or merged.
Heads: main 35ead8b6, fable_v4c 9a87c095, cleanup/v4c-slim-rc a807ebb8. Merge base of main and fable_v4c: 9651a900 (2025-10-12).
Method: ancestry (git merge-base --is-ancestor), git cherry (patch ids), GitHub REST pull-request list (18 PRs), and, because several
PRs were squash merges and the fable line was partly hand-merged, content checks (new files present, identifiers present, added-line presence).
Complete commit-level digest by theme: docs/dev/branch_developments.md.

## 1. Headline findings

1. **Every fable_v* branch from fable_v1_test to fable_v4b_u / fable_v4b_gbt is an ancestor of fable_v4c.** So are belem25, gershgorin, gershgorin_v2,
   lda_u, libXC and print_hamiltonian. fable_v4c = v4b + v4b_u (LDA+U) + v4b_gbt (GBT); it has no commits of its own. cleanup/v4c-slim-rc contains fable_v4c.
2. **The belem25_strux* / strux / ldim / hall / sphinx line is in the fable line by squash or hand-merge, not by ancestry.** PR #6 squashed
   belem25_strux_ldau_kspace_nc into belem25 (tree diff between the two is empty). strux and ldim were merged by hand into belem25 (commit 7d6e682d).
   All 24 new source files of that line exist in fable_v4c; theta_ss (the 'hall' cone angle) exists in fable_v4c.
3. **PR #7 (fable_v2-main-merge) brought only 5 files of main into the fable line** (bands.f90 +214, calculation.f90, control.f90, exchange.f90,
   namelists/control.f90). Of main's 34 non-merge commits since the merge base, only three have most of their source changes in fable_v4c
   (quadrupole implementation eb8966e2 88%, exchange output format 5f426578 90%, bands print fix 71c97ab4 50%). The Kubo rewrite of 2026-06-11
   (0a7d6843, +1546 lines hamiltonian.f90, +503 recursion.f90) is not in fable_v4c (5% of its added lines), and neither is anything after it
   (table in section 4). pol_alpha and setup_kubo_operators exist only on main; js_alpha only in fable_v4c.
4. **Work found on no integration branch (candidates to lose if branches are deleted):**
   - main side: feat/force-nonmagnetic (1 commit, per-atom forced nonmagnetic SCF), fix/control-lattice-memory-cleanup (8 commits, control/lattice/charge/hamiltonian
     allocation fixes, tests/memory), codex/mn2f2-* (PR #10 and #11 open: magnetic GF symmetry diagnostics, source/gf_audit_io.f90),
     audit/conductivity-legacy-20261103 (9 CI benchmark commits).
   - fable side: fable_v1_kernels_hamgpu (15 commits, Hamiltonian build on GPU; fable_v4c uploads a CPU-built Hamiltonian), fable_v1_kernels (1),
     fable_v1_ctest (2, CUDA kernel tests), fable_v1_flat (2, flattening hoh into H), convoluted (5 commits on top of belem25_strux_ldau_kspace: convolution/sparse backend plugin),
     reciprocal (1 commit: S(r)->S(k) k-space strategy; its routines are absent from belem25 and fable_v4c, possibly superseded by the strux_lib based k-space (inferred)).
5. **Both sides developed transport independently**: main = Kubo operator rewrite, quadrupole/octupole operators, full HOH actions, physical units (Ramon Cardias, Jun-Nov 2026);
   fable = KPM / Kubo-Bastin transport with CPU/GPU kernels (conductivity.f90 grew to 780 lines; recursion split into recursion_core/chebyshev/haydock/transport;
   Aug 2026). 15 source/build files were changed on both sides (list in section 5); fable_v4c keeps the test UnitKpmTransport and the validation test Val09KuboBastinTransport.

## 2. Pull requests (REST)

PRs 1-9 and 12-18 are merged; 10 and 11 are open. Into main: #1 init_scf, #2 actions, #3-5 update_tests, #8 tb-orbital-effective-field, #9 octupole-accumulation,
#12-#18 (transport, units, labels). Into belem25: #6. Into fable_v2: #7. #10 (codex/mn2f2-gf-symmetry-audit -> main) and #11 (codex/mn2f2-sbar-reciprocity-test -> the audit branch) are open.
Only PRs #1, #2, #3-5, #6, #7, #9 were squash merges (branch tips are not ancestors of their targets).

## 3. All branches by status

### A. Merged into main

| branch | tip date | author | merged into / how | note |
|---|---|---|---|---|
| fix/response-label-filenames | 2026-11-05 | Ramon | main (merge commit, PR #18) |  |
| fix/conductivity-units | 2026-11-04 | Ramon | main (merge commit, PR #17) |  |
| fix/kubo-periodic-velocity | 2026-11-04 | Ramon | main (merge commit, PR #16) |  |
| examples/fccpt-current-conductivity-inputs | 2026-11-03 | Ramon | main (merge commit, PR #15) |  |
| fix/masked-hoh-original-routines | 2026-11-03 | Ramon | main (merge commit, PR #14) |  |
| refactor/transport-existing-style | 2026-11-03 | Ramon | main (merge commit, PR #13) |  |
| fix/transport-multipoles-hoh | 2026-11-03 | Ramon | main (merge commit, PR #12) |  |
| feature/octupole-accumulation | 2026-08-31 | Ramon | main (PR #9, squash: patch present in main) |  |
| fix/tb-orbital-effective-field | 2026-08-24 | Ramon | main (merge commit, PR #8) |  |
| update_tests | 2026-03-12 | Anders | main (PRs #3, #4, #5, squash: 29/29 files present in main) |  |
| actions | 2025-02-11 | Anders | main (PR #2, squash) | CI workflow update |
| init_scf | 2025-02-03 | Anders | main (PR #1, squash) | cold-start flag is in main and in fable_v4c |

### B. Merged into fable_v4c

| branch | tip date | author | merged into / how | note |
|---|---|---|---|---|
| fable_v4c | 2026-11-05 | Anders | slim branch cleanup/v4c-slim-rc (ancestor) | the integration branch: v4b + v4b_u + v4b_gbt |
| fable_v4b_u | 2026-11-04 | Anders | fable_v4c (ancestor) |  |
| fable_v4b_gbt | 2026-11-04 | Anders | fable_v4c (ancestor) |  |
| fable_v4b | 2026-10-16 | Anders | fable_v4c (ancestor) |  |
| fable_v4 | 2026-10-03 | Anders | fable_v4c (ancestor) |  |
| fable_v3 | 2026-08-26 | Anders | fable_v4c (ancestor) |  |
| fable_v2 | 2026-07-28 | Anders | fable_v4c (ancestor) |  |
| fable_v1 | 2026-06-29 | Anders | fable_v4c (ancestor) |  |
| fable_v1_hoh | 2026-06-15 | Anders | fable_v4c (ancestor) |  |
| fable_v1_batched | 2026-06-13 | Anders | fable_v4c (ancestor) |  |
| fable_v1_test | 2026-06-12 | Anders | fable_v4c (ancestor) |  |
| belem25 | 2026-06-10 | Anders | fable_v4c (ancestor) |  |
| gershgorin_v2 | 2025-10-25 | Anders | fable_v4c (ancestor) |  |
| lda_u | 2025-10-15 | Anders | fable_v4c (ancestor) |  |
| libXC | 2025-10-11 | Anders | fable_v4c (ancestor) |  |
| print_hamiltonian | 2025-10-09 | Anders | fable_v4c (ancestor) |  |
| gershgorin | 2025-02-04 | Anders | fable_v4c (ancestor) |  |

### C. Content in the fable line by hand-merge or squash (tip is not an ancestor)

| branch | tip date | author | merged into / how | note |
|---|---|---|---|---|
| fable_v2-main-merge | 2026-07-05 | Anders | fable_v2 (PR #7, squash a8ec82e9): brought main up to 2026-06-15 into the fable line |  |
| belem25_strux_ldau_kspace_nc | 2026-06-10 | Anders | belem25 (PR #6 for the _nc tip, squash); hand-merge of strux and ldim, commit 7d6e682d | its 24 new source files all exist in fable_v4c |
| belem25_strux_ldau_kspace | 2026-05-31 | Anders | belem25 (PR #6 for the _nc tip, squash); hand-merge of strux and ldim, commit 7d6e682d | its 24 new source files all exist in fable_v4c |
| belem25_strux_ldau | 2026-05-26 | Anders | belem25 (PR #6 for the _nc tip, squash); hand-merge of strux and ldim, commit 7d6e682d | its 24 new source files all exist in fable_v4c |
| belem25_strux | 2026-05-25 | Anders | belem25 (PR #6 for the _nc tip, squash); hand-merge of strux and ldim, commit 7d6e682d | its 24 new source files all exist in fable_v4c |
| strux | 2026-05-15 | Anders | belem25, hand-merged (7d6e682d); 16 new source files exist in fable_v4c |  |
| ldim | 2026-04-08 | Anders | belem25, hand-merged (7d6e682d); flexible basis present in fable_v4c |  |
| sphinx | 2026-03-03 | Anders | belem25 line (doc commit 552c11a4 is on the _nc tip, which PR #6 squashed into belem25) |  |
| hall | 2025-12-12 | Anders | belem25 line ('Improved HALL treatment. Added simple spin-spiral support', 6ca7ff45); theta_ss exists in fable_v4c | not in main |

### D. Not merged anywhere (decision needed)

| branch | tip date | author | merged into / how | note |
|---|---|---|---|---|
| reciprocal | 2025-11-05 | Anders | not merged | 1 commit: k-space built as S(r)->S(k) instead of H(r)->H(k) (reciprocal.f90); its routines fourier_transform_structure_constants and mount_kspace_hamiltonian_from_sk are absent from belem25 and fable_v4c; may be superseded by the strux_lib based k-space (inferred, not checked) |
| audit/conductivity-legacy-20261103 | 2026-11-03 | Ramon | not merged | 9 CI benchmark commits (Pt conductivity, old vs new Kubo), 5 files |
| codex/mn2f2-gf-symmetry-audit | 2026-10-21 | Ramon | not merged; PR #10 open into main | magnetic GF symmetry diagnostics (source/gf_audit_io.f90 is new) |
| codex/mn2f2-sbar-reciprocity-test | 2026-10-21 | Ramon | not merged; PR #11 open into the audit branch | 7 more commits on top of the audit branch |
| feat/force-nonmagnetic | 2026-09-30 | Ramon | not merged | 1 commit: per-atom forced nonmagnetic SCF (force_nonmagnetic); not in main or fable_v4c |
| fix/control-lattice-memory-cleanup | 2026-09-15 | Ramon | not merged | 8 commits: control/lattice/charge/hamiltonian allocation fixes + memory audit doc; not in main or fable_v4c |
| fable_v1_kernels_hamgpu | 2026-07-02 | Anders | not merged | 15 commits: Hamiltonian build on GPU phases 1-3; fable_v4c uploads a CPU-built Hamiltonian instead |
| fable_v1_kernels | 2026-06-30 | Anders | not merged | 1 commit: k-space GPU kernel exploration; contained in fable_v1_kernels_hamgpu |
| fable_v1_flat | 2026-06-13 | Anders | not merged | 2 commits: experiment 'flattening hoh into H', 'removing U from hoh' |
| fable_v1_ctest | 2026-06-12 | Anders | not merged | 2 commits: CUDA kernel tests (2026-06-12) |
| convoluted | 2026-06-09 | Anders | not merged | belem25_strux_ldau_kspace plus 5 commits on a convolution/sparse backend plugin (2026-06-09) |

### E. Stale or WIP, not merged

| branch | tip date | author | merged into / how | note |
|---|---|---|---|---|
| j_representation | 2025-02-04 | Anders | not merged | WIP 2025-02: TB parameters for JIJ |
| backtick | 2024-07-03 | Ramon | not merged | 2024, 7 commits; 28 of 41 changed files already match main |
| prettify | 2023-11-27 | Anders | not merged | 2023, 11 commits; stale |

Notes: the Kubo-related fix/* and examples/* branches are ancestors of main and of each other; main itself is 43 commits (9 merges) ahead of fable_v4c by ancestry.
Branches in class A and B are fully contained in their targets and are candidates for deletion by the repository owners; class C branches are redundant with their squash targets.
Class D and E should be kept until a decision is taken.

## 4. main commits since the merge base: presence of their source changes in fable_v4c
(added source lines of length >= 12 that appear verbatim in the same file on origin/fable_v4c; n/a = no source lines changed)

    sha      date        src +lines  present  subject
    f53cdece 2025-10-24           0      n/a  Graphene example added
    b4e5aada 2025-10-24           0      n/a  graphene clear dir
    7c9f1b82 2025-11-15           0      n/a  Remove unused parameters from C1_out.nml
    384c6b2a 2025-11-15           0      n/a  Clean up zero values in C2_out.nml
    0a6262a8 2026-03-07           0      n/a  Update tests (#3)
    79f41aea 2026-03-08           7       0%  Fixed pre-allocation before namelist reads for increased compiler stability
    847802ed 2026-03-08           3       0%  Improving the preallocation from previous commit
    f9eadd7c 2026-03-09          37      16%  small updates on dimensions and prints
    09a2d207 2026-03-11           1       0%  Update tests (#4)
    c183dadc 2026-03-12           0      n/a  Update tests (#5)
    eb8966e2 2026-06-11          96      88%  quadrupole implementation
    00835589 2026-06-11           5       0%  small modification on surface clust calculation
    5f426578 2026-06-11          10      90%  updated format of exchange output
    2e1aabef 2026-06-11           0      n/a  self changes
    0a7d6843 2026-06-11         464       5%  update on kubo linear response formalism. Need do update conductivity examples s
    abb8ceba 2026-06-11           0      n/a  teste
    69863d94 2026-06-11           0      n/a  teste done
    71c97ab4 2026-06-15           2      50%  fix print bug on bands.f90
    f951873b 2026-07-06         111       5%  added quadrupole current and accumulation to the kubo operators
    495709e7 2026-08-21          48       0%  lsham build fixed to accept negative xi
    22b9de89 2026-08-24          55       2%  Fix TB orbital-resolved effective field
    49c0031b 2026-08-31          73       4%  Add orbital octupole accumulation operators (#9)
    e604f402 2026-10-01           0      n/a  docs: add comprehensive RS-LMTO-ASA README
    c1ffdfd3 2026-11-03           7       0%  Fix quadrupole normalization and traceless diagonal components
    1522cea4 2026-11-03          85       2%  Add pure orbital octupoles and conventional current response
    cd0a2432 2026-11-03         214       0%  Implement full HOH actions for all Kubo response operators
    c9165744 2026-11-03         145       5%  Keep corrected transport products in existing recursion and BLAS style
    e692578f 2026-11-03         319       3%  Restore explicit legacy-style transport loops with complete HOH corrections
    41f3a8ac 2026-11-03           0      n/a  Update fccPt conductivity examples to current response inputs
    a741176f 2026-11-04          50       2%  Fix periodic bond displacements in matrix-free Kubo velocities
    a2469d77 2026-11-04         216       1%  Fix conductivity normalization and add explicit physical response units
    dbf81309 2026-11-04          36      19%  Preserve atomic unit-cell traces and legacy numerical conventions in units outpu
    dc54ef0c 2026-11-05           1       0%  Use atomic labels for resolved response filenames
    9072004c 2026-11-05           0      n/a  Fix Fortran MPI CMake integration and document oneAPI LP64 builds
    

Reading: 0-20% means the change is not in fable_v4c (it can still be superseded by a different implementation; that was not checked). The four 2026-03 preallocation/dimension
commits and the 2026-06-11 surface-cluster tweak (00835589) are not in fable_v4c and touch the code fable rewrote for ldim/strux (lattice.f90, charge.f90, hamiltonian.f90).

## 5. Conflict zones for a main merge

Files changed on both sides since the merge base: 96 (main changed 132 files in total, fable_v4c 3469). In source/ and the build:
CMakeLists.txt, cmake/SetFortranFlags.cmake, source/bands.f90, source/calculation.f90, source/charge.f90, source/conductivity.f90, source/control.f90, source/exchange.f90, source/hamiltonian.f90, source/include_codes/namelists/control.f90, source/include_codes/namelists/hamiltonian.f90, source/include_codes/namelists/self.f90, source/lattice.f90, source/recursion.f90, source/self.f90.
A trial merge of main into the slim branch conflicted in 43 files, 52 hunks in 12 source files (self.f90 among them). Highest overlap: hamiltonian.f90 (fable split it into hamiltonian_build/ccor/hubbard/paoflow_io),
recursion.f90 (fable split it into recursion_core/chebyshev/haydock/transport), conductivity.f90, bands.f90, lattice.f90, and the namelists (hamiltonian: pol_alpha/pol_beta on main vs js_alpha/jl_alpha, q_ss, theta_ss, ccor_2c, Hubbard keys on fable).
Files main changed that fable did not touch (these should transfer without textual conflict): docs/octupole_response.md, example/transport/*, example/conductivity/fccPt inputs, source/include_codes/namelists/conductivity.f90,
tests/benchmarks/periodic_velocity and response_units, tests/regression/test_{octupole,quadrupole,hoh_transport,response_units}*.py and their drivers.
The slim branch deleted tests/benchmarks/ and most tests/scf/references, so main's edits to those paths will appear as modify/delete conflicts.

## 6. Caveats

- "Not merged" means not reachable and not patch-equivalent; content can still exist by hand-merge. Hand-merges were checked by new-file and identifier presence for the belem25 line, hall, reciprocal, init_scf and force-nonmagnetic, not line by line.
- Section 4 compares added lines only; it cannot show that a change was replaced by another implementation.
- The REST pull-request list shows 18 PRs; forks and branches without PRs are not visible. Tags (8, including archive/* of 2026-11-05 and lr-campaign-archive*) were not analysed.
- fable_v4c was not slimmed by the owners; the slim branch is cleanup/v4c-slim-rc.
