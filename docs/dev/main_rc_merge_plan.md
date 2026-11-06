# Plan: merge `main` into the slim fable line (`main_rc`)

Agreed with the owner on 2026-11-06. Inputs: `docs/dev/branch_status.md` (which branches are merged where), `docs/dev/branch_developments.md` (what each side developed).

## Decisions

1. **Branch.** New branch `main_rc` from `cleanup/v4c-slim-rc` (`88a860eb`; check the head first, nobody should push there meanwhile). Push only to new `*_rc` branches, never to `main` or `fable_v4c`.
2. **Unmerged branches are settled.** `feat/force-nonmagnetic` is already on `main` (PR #19, `main` = `2fccbffc` when written). No other unmerged branch is merged.
3. **Transport:** the framework comes from `main` (it has priority); the Chebyshev/KPM kernels come from the slim branch where they apply. Where both sides implement the same quantity, `main`'s definition wins. Everywhere else the slim layout wins, with `main`'s semantic fixes ported into it.
4. **`pol_alpha`/`pol_beta` replacing `js_alpha`/`jl_alpha` is deliberate** (more general than S and L). Migrate every deck and test on the slim side; do not add compatibility keys.
5. **Keep all four fable-only transport features if they cause no redundant functionality:** Kubo-Bastin GPU reconstruction, route-agnostic estimators (the triad tests), k-space conductivity moments, GPU stochastic transport batching. Where one duplicates something `main` now does, stop and ask.
6. **Tolerances.**
   - `main`'s loops against an equivalent framework or kernel: start at 1e-8 (relative, per moment and per conductivity element), expect better.
   - A single-precision (SP) kernel replacing a double-precision one: 1e-5.
   - Imaginary part of the total conductivity (`cond_total`) for periodic systems: about 1e-14. Test case: fcc Pt with periodic boundaries, `example/conductivity/fccPt`.

## Method

- Port `main`'s final state function by function into the slim file layout. Do not replay commits: `main`'s history refactors and then restores (`c9165744`, `e692578f`). Do not run a plain `git merge`: a trial merge gave 43 conflicting files and 52 hunks in 12 source files, because the slim branch split `recursion.f90` and `hamiltonian.f90`.
- Where the code lives. `main`: operator application in `source/recursion.f90` (`setup_kubo_bonds`, `setup_kubo_operators`, `apply_kubo_operator`, `kubo_ham_vec_matmul`, `kubo_operator_vec_matmul`, `compute_moments_stochastic`). Slim: `recursion_transport.f90` (`velo_vec_matmul`, `velo_hoh_vec_matmul`, `ham_vec_matmul`, `compute_moments_stochastic`), `chebyshev_fast.f90` (`cheb_moments_*_fast*`), `recursion_chebyshev.f90`, `recursion_core.f90`, `kpm_profile.f90`, `source/cuda/`, `conductivity.f90` (Kubo-Bastin contraction).
- **Open technical question (spike first).** `cheb_moments_stochastic_fast` runs the whole double Chebyshev loop itself and takes velocity-operator arrays (`v_a`, `v_b`, `vo_a`, `vo_b`); it has no hook to apply an operator in the middle. `main` applies operators through generic slots (quadrupole, octupole, conventional current, full HOH actions). Check whether `main`'s operators reduce to bond-resolved arrays (`setup_kubo_bonds` suggests so). If not, "kernels from the slim branch" is limited to the `H|psi>` application and the recurrence.

## Steps and gates

0. **Baselines.** Build `main` and record its transport outputs: tests in `tests/regression/test_{hoh_transport,quadrupole_operators,octupole_operators,response_units}.py`, benchmarks `tests/benchmarks/periodic_velocity` and `response_units`, the fccPt conductivity examples. Record the slim branch's known failures (`tests/KNOWN_ISSUES.md`) so they are not taken for merge damage. Use the CI configuration, `-DENABLE_MARCH_NATIVE=OFF`, for anything compared against references.
1. **Small independent `main` changes:** March 2026 preallocation commits (`79f41aea`, `847802ed`, `f9eadd7c`), surface-cluster tweak (`00835589`), negative-xi `lsham` (`495709e7`), TB orbital-resolved field (`22b9de89`), `force_nonmagnetic` (#19), README, MPI CMake (`9072004c`). `self.f90` stays off-limits in `CLAUDE.md`, but the owner allows resolving its merge hunks for this task (2026-11-06): apply only what the port needs, make no other edit in the file, and list every `self.f90` hunk and how it was resolved in the report.
2. **Transport framework from `main`** on the slim branch's `legacy` Chebyshev backend. Gate: `main`'s transport tests and benchmarks pass on `main_rc` within the tolerances above; the slim quick and standard tiers still pass.
3. **Kernel substitution behind a backend switch.** The default stays `main`'s loops until a kernel matches them within tolerance and a timing comparison on one problem shows it is faster. Report both numbers.
4. **Register `main`'s transport tests** with `rslmto_add_test` in a tier by measured runtime (`quick` stays <= 2 min at `-j4`; each quick test <= 10 s). Keep the benchmark data small.
5. **Record ancestry last.** Only after the gates pass: `git merge -s ours origin/main` on `main_rc`, so later `main` merges conflict only on new changes.

## Rules for references

`main`'s new conductivity normalisation and physical units move the numbers behind the slim branch's triad sigma test, `Val09KuboBastinTransport` and `UnitKpmTransport`. Rebase those references only with the transport developer's recorded approval and a one-line reason in the commit and in the reference file. Never edit a reference or tolerance to make a check pass (`CLAUDE.md`).

## Environment notes (cloud container)

- No toolchain by default: `apt-get install gfortran liblapack-dev libopenblas-dev libxc-dev libfftw3-dev libsymspg-dev strace`; `pip install f90nml numpy pyyaml`.
- Configure tests with `-DRUN_UNIT_TESTS=ON -DRUN_EXAMPLE_TESTS=ON -DRUN_REG_TESTS=ON -DENABLE_LIBXC=ON`. A `-march=native` build on FMA hardware misses some bitwise and 1e-10 checks (`tests/KNOWN_ISSUES.md`).
- Pushing needs the repository attached with push access; push only the `*_rc` branch. A stop hook asks to push `cleanup/v4c-slim` (a branch name that does not exist on the remote); it can be ignored.
- `self.f90` and `symbolic_atom.f90`: the exception above covers `self.f90` merge hunks only; nothing else in either file.
- Commit messages: one line, no co-author lines (`CLAUDE.md`).
