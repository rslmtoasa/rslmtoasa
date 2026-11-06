# Developments digest: `main` vs the fable family (fable_v4c)

Sources: per-layer commit logs (non-merge commits), `git show --stat` on a few commits, `docs/DECISIONS.md`, `ROADMAP.md`, `DEVELOPER_MAP.md` of `origin/fable_v4c`.
Counts per theme are tallies of commit subjects and are approximate (`~`) where themes overlap. Anything not stated in a log, doc or diff stat is marked "(inferred)".
Layer commit counts: L00 106, L01 3, L02 15, L03 11, L04 12, L05 6, L06 108, L07 192, L08 214, L09 143, L11 179, L12 116 (sum 1105 in fable_v4c, excluding merges), main 34. L13 (fable_v4c own commits) is empty.

Key finding for the merge decision (measured on `origin/fable_v4c`): PR 7 ("Fable v2 main merge", a8ec82e9) changed 5 files only (bands.f90 +214, calculation.f90, control.f90, exchange.f90, namelists/control.f90). The Kubo linear-response rewrite on main dated 2026-06-11 (0a7d6843: hamiltonian.f90, recursion.f90, namelists) is NOT visible in fable_v4c: `pol_alpha`, `setup_kubo_operators` are absent from `source/` there, and the `&hamiltonian` namelist still lists `js_alpha, jl_alpha`. So "main up to 2026-06-15 is in fable_v2" holds for the quadrupole/bands/exchange-output/lattice small fixes, not for that Kubo rewrite.

---

## 1. Layer by layer

### L00 pre-belem25 stack (gershgorin_v2, lda_u, libXC, print_hamiltonian, gershgorin)
2024-06-26 .. 2025-10-25, 106 commits (authors in the log: ViktorFr and Emil for LDA+U in 2024, Anders Bergman for 2025).
- LDA+U (DFT+U) first implementation: hubbard_u module, U/J inputs, Slater integrals (d, then f), double counting, density matrix, +V intersite correction, self-consistent U, U in Ry; ~54 commits, 2024-06-26 .. 2025-01-09. SHAs a0b5e57f, 786fe0bc, 81920120, c195047d.
- LDA+U on impurities (Co host / Fe imp, atom.nml read/write): ~8 commits, 2024-12 .. 2025-01. b447d85b, 70b70e70, 28f5c124.
- K-space / reciprocal code starts: basic k-space, spglib interface, symmetry.f90, band-structure FT, tetrahedron DOS, band/DOS namelist tasks, examples; ~20 commits, Sep-Oct 2025. 87726190, 09f3c0c9, c88107f0.
- K-space SCF (single- then multi-site, Fermi-level fixes, Blochl correction/integration, LDA+U in k-space, linear mixing of density matrices): ~25 commits, Oct 2025. 52321898, b3e3bef4, d8a0d76d, 30054302.
- libXC integration (LDA/GGA, X and C functionals, gradients, CMake): ~12 commits, 2025-09-21 .. 2025-10-11. 72cb73ee, 2bd6a9f5, 6ee108ec.
- Hamiltonian bounds (Gershgorin) and Chebyshev bounds/energy scaling: 5 commits. 4f2f629c, e28d94f5, 6da80bff.
- Performance: optimised reciprocal routines, OpenMP-enabled H_all mounting, build_locham optimisation: 5 commits. 04f14b1b, 049b5de5, df87f18e.
- Spin-spiral vector introduced (`q_ss`, first appearance of the spiral input): 1 commit. 32ace7f9.
- Hamiltonian printing routine, FeCo / NiO / Co-Fe examples: ~5 commits. e3d06f60, 5a129cf6.

### L01 belem25 (3 commits)
2025-11-25 .. 2026-06-10.
- One large squash "Belem25 strux ldau kspace nc (#6)": 169 files, +27347/-7461 (stat). Adds source/include_codes/strux_lib (strux_api, bessel, geometry, lmto47_mode, tb, ...), docs/ Sphinx tree (keywords, theory, user_guide), tests/ harness (run_test.py, tests/scf cases + references), `.github/workflows/tests.yml`, and large rewrites of hamiltonian.f90, bands.f90, reciprocal.f90, lattice.f90, recursion.f90, self.f90, green.f90; basis.f90 (+129) added. 1355a50c. Contents beyond the stat are (inferred) from the PR title: strux structure constants, LDA+U, k-space, non-collinear.
- Housekeeping: untracked files removed (82dc905d), quick fix (aa0840df).

### L02 fable_v1_test (15 commits)
2026-06-10 .. 2026-06-12. Sparse initiative and accelerator backends.
- Restart of the sparse matrix work, sparse.f90 / sparse_matrix_mod, hamiltonian_cpu_mod, interface drivers, FFT backend hook, cpp backend tests: ~6 commits. c3d5a167, b95ea55a, 66773b89.
- First CUDA backend (AtomicAdd, low precision kernels), CUDA DOS routine wired in: ~7 commits. 2fe16378, df3eeccf, bed1378a.
- Fortran Chebyshev kernel in half precision, commit message cites 30% gain: 1 commit. 6ee9ca02.

### L03 fable_v1_batched (11 commits)
2026-06-13.
- Batched CPU kernels, MKL batched Chebyshev (chebyshev_mkl_batch.cpp), ifx + MKL compile fixes, memory-allocation improvements for batched kernels: ~8 commits. 410a676e, f5476b32, 3b36fdb1.
- CUDA kernel tuning (non-Chebyshev), CMake build modes, original Chebyshev kernel added to the new backends: 3 commits. 46f76099, 0f74f779, bc4c98cd.

### L04 fable_v1_hoh (12 commits)
2026-06-14 .. 2026-06-15.
- `hoh` support in the new CPU/GPU kernels (fast CPU block recursion hoh, GPU block kernels, fp32 CUDA kernel): ~7 commits. a015b053, 1ba2fd4d, ff0e06df.
- Removal of cpp_plugin, c-library stub for non-CUDA builds, Chebyshev routines moved into a driver: 3 commits. e5ac5fc9, 9284e898, cfb9ae28.
- GPU port of DOS/Green-function routines (bgreen on CUDA); hoh support in reciprocal mode: 3 commits. ae7981bd, be9060ee, 95972953.

### L05 fable_v1 (6 commits)
2026-06-15 .. 2026-06-29.
- CCOR (ccor_2c) new attempt, workflow clean-up and docs; CCOR reported working for occupied bands: 3 commits. d762e8eb, a5f99a0d, b6299e23.
- Tetrahedron k-space algorithm improvements, k-space energy windowing reverted: 3 commits. 1da58a8d, df483442.

### L06 fable_v2 (108 commits)
2026-07-03 .. 2026-07-28.
- Refactor T1-T13 (backend regression coverage, dead recursion code removed, Chebyshev cache moved into handle, shared CPU kernel loop, green-function driver dedup, CUDA plugin context, large modules split into submodules, MKL options collapsed, CMake hygiene): ~15 commits. 542e95a7, cbb53457, 864f8242.
- Test/CI/docs plan P1-P10 (Lanczos in SCF suite, PAOFLOW import, orbital moments, Hubbard SCF, k-space CCOR/HOH DOS, serial-vs-MPI schema, GPU/MKL CI strategy, quick/full split, nightly, Windows job, DEVELOPER_MAP.md and root CLAUDE.md, reference-regeneration workflow): ~27 commits. a85ba580, ae4261b3, eced2068.
- Squash merge of main into fable (PR 7): 1 commit. a8ec82e9.
- GBT / spin-spiral work in real space and first k-space GBT transform (first and HOH/CCOR), `q_ss` unit fixes, frozen_magnon post-processing: ~11 commits. cbbb3208, aabb842d, 130b20dc, 13f3d171.
- B2/B3/B5 k-space Green functions (reciprocal_green skeleton, Lehmann and Dyson backends, eta ladder, torque families, route dispatch), Bloch spectral functions A(k,E), k-space conductivity moment generator, route-agnostic J_ij/sigma tests: ~20 commits. 5e78e9cc, def7d8b6, b1c9b4ff, 52374509, 75aed0b1.
- Gilbert damping (do_damping wiring, SOC-in-G alpha), simpson_f out-of-bounds fix (cause of J_ij NaN), Chebyshev on-site G_ii normalisation fix: 5 commits. 22b42c85, 8b429287, dff13698.
- B6/B7 surface and interface electrostatics: l=1 dipole electrostatics via Madelung, madl2d overflow guard, per-site region registry, `interfacepot` two-sided electrostatics, vacuum lead parameter generator, alignment solver, `calctype='L'` interface mode, Cu oracle: ~25 commits. 8c28ef62, 449b5363, 82334a1e, 97f1e0e4, 8883c608.
- Executable-stack fix in Chebyshev kernels; docs for B1-B7 keywords, module overview, code-features inventory: ~6 commits. 4a345938, 00abdf58.

### L07 fable_v3 (192 commits)
2026-08-08 .. 2026-08-26.
- GBT work packages WP1-WP10 (operator gates, solver unification, HOH/CCOR covariance, rotating-frame SCF, q little-group BZ reduction, legacy deletion): ~15 commits. 56f8ca7b, fa8d7fc4, 98df57f7, a6f935cd.
- TDDFT / linear-response (TDDFT-00..11: response conventions, arbitrary-k eigenpair service, transverse Kohn-Sham susceptibility, ALSDA kernel, Goldstone handling, enhanced susceptibility and magnon analysis, GF-based backend, pair-potential operators, Ward diagnostics): ~41 commits, almost all 2026-08-10. 007f72fe, 25c4ae9a, 806c8190, 38d5bac4.
- Reciprocal-code refactor RF-01..06 and GC-01..06 (parallel state, k-point ownership, batched H(k) workspaces, execution backends, tile scratch reuse, GPU gate review) plus Fermi-exponent overflow guard, zheev safeguards: ~15 commits. 41fc840f, c57ae2f6, 5339fffe.
- Accelerators ACC and KPM (reciprocal CUDA backend skeleton and Hermitian eigensolver, batched CUDA eigensolver, CUDA Lehmann Green functions, MPI rank-to-device, persistent KPM CUDA paths, GPU BLAS Kubo-Bastin reconstruction, stochastic-trace batching, CPU/GPU benchmark harness and campaigns, SCF performance campaign): ~42 commits. f015b45e, 412c703a, 65ac455b, 5a1ace1a, 348e8a2a.
- Phase II validation records (reciprocal SCF, libXC SCF, Lehmann/Dyson, advanced exchange, damping and inertia, spin/orbital Kubo-Bastin, constraining fields, spin dynamics, interface/vacuum electrostatics, GBT vs supercells, bcc Fe and fcc Ni TDDFT): ~20 commits. 656eb66a, ebf35196, 6d3e0064, 124163e8.
- Test framework: manifest-driven tests with labels/tiers, strux_lib as default test structure backend, legacy-vs-strux structure-constant equivalence tests, lean Si sp fixtures, orbital-basis dimension tests, LDA+U algebra and libXC equivalence tests, lean spin/orbital KPM transport tests: ~13 commits. ebf2b59c, 3a07a3c6, 5c0643e7, 1d2a061a, 0efbd3ae.
- CI: runner-native reference artifacts, macOS/Windows (UCRT64)/MPI/spglib fixes, unit tests enabled in CI: ~19 commits. b2e2ab20, 615501c2, db02480a.
- Submodule splits (calculation preprocessing/reciprocal, green lifecycle, Chebyshev and recursive Green, self modern reciprocal helpers, legacy 2D Madelung isolated, interface electrostatics out of charge, damping/inertia out of exchange): 9 commits. 809b9a75, 585fb75c, ff787b7e.
- Structure constants and basis: strux custom-supercell fix and faster local solves, orbital basis generalisation in k-space (STAB-04), Lanczos spin-polarised fix, reciprocal DOS normalisation fix: ~6 commits. 66f3d1bd, 9cf59cdb, ec5cf632, ad092861, b63160b2.

### L08 fable_v4 (214 commits)
2026-08-30 .. 2026-10-03. Linear-response campaign, later archived.
- XC reconciliation (legacy vs libXC selectors, fixed-density tests, radial GGA and libXC derivatives repaired, zero-density guard, production libXC contract, atomic SCF controls restored): ~15 commits. 91be9db2, 80390cc9 (DECISIONS: shifts bcc-Fe total energies by 6e-6 to 1.5e-5 Ry), 65dc136d.
- TDDFT framework consolidation (chi0 backend unification, k-space and real-space Green-function chi0, contour integration, loss-matrix magnon analysis, Fe/Ni, ferri/AFM chirality, longitudinal charge-spin response, Ward diagnostics, Goldstone correction): ~65 commits, Sep 1-15. d1ddd3ec, 5e79297c, 8d2cb284, 6ab01357.
- DRESP series (projected reciprocal bare susceptibility, KL/LKAG bridge, finite-q exchange `exchange_q`, Mills RPA, Juelich response, ALSDA spatial comparison, Ward/field-mapping and rigid-rotation checks, product basis, scalar-relativistic vertex audits): ~53 commits, Sep 16-23. d139bde2, a7a8ebea, fecc4f2e, 70eb73b4.
- Refactor LR-REF-01..05 (DRESP/Ward/tangent modules deleted, `linear_response_mod` with basis/bare/kernel_dyson/rotation submodules, `&tddft` replaced by `&linear_response`, test-oracle library, merged multi-case unit tests): ~30 commits, 2026-09-27/28. 4691b25c, 63cf9cd3, 5a0b47fc, 9ae10a9e.
- LR-CLOSE-01..18 and LR-METHOD (rotation workflow, baseline harness, Mills 1U replaces fitted Mills route, public API): ~21 commits. e6180ad3, 6a87b9b2, 3fa43a83.
- Docs reorganisation DOC-00: docs/conventions, docs/validation, ROADMAP.md, DECISIONS.md created; prompts and campaign reports deleted; docs policy added to CLAUDE.md: 3 commits. 12c89fae, 592ebc4d, 9d6ffa3d.
- Recovered GBT/reciprocal and SCF changes from pending state: 2 commits. 44b69a00, 483e656c.
- Last commit: ALSDA kernel contract repair, Ward closure still diagnosed (not established per DECISIONS). e95a17cf (tag lr-campaign-archive-2026-10).

### L09 fable_v4b (143 commits)
2026-10-03 .. 2026-10-16.
- Removal of the archived linear-response campaign from fable_v4b (source, tests, decks, `post_processing` values stopped at startup), docs update, new "lr-archive" test label and isolation check: ~10 commits. 63202725, 955f6ced, 54dd21b3, dd3d8ca1.
- Test-suite triage: `quick` label and per-commit quick tier, disabled/failing regression tests recorded in KNOWN_ISSUES, legacy Lanczos/Block/Chebyshev harness retired and replaced by run_matrix cases, Triad J_ij -O0 vs -O2 diagnosis, references regenerated, CI Python selection: ~35 commits. 32772694, 06407dad, ea98f1b9, b925efae.
- Exchange energy integral fix: integrate to grid point enpt+1 instead of fermi; `&energy temperature` (default 0 K) with a simpson_f Sommerfeld unit test: 3 commits. 2cb6f554, 7be98f7a.
- New `spin_response` (B11 restart): C0 kernels, prepare/interaction/run drivers, bcc Fe deck and LSWT oracle, several sites with site-dependent Bloch phase, 24^3 reproduction baselines, OpenMP over omega, fcc Ni deck and LSWT oracle, channels s/p/d/f, experimental ALSDA kernel, multi-site adiabatic dispersion, AFM/ferrimagnet oracles, bcc Cr deck, force-theorem frequencies, SOC kernel (NC-2 negative result recorded): ~65 commits. 1cc43b20, 6fe83577, add66d2a, 02582da1, c24fc42e, 36e1cb3c.
- On-site L.S added to the first-order k-space Hamiltonian for nsp = 2 and 4: 1 commit (with a re-baselined reference). f005b770, 72473205.
- Fixed-U LDA+U remediation (spherical tensor, complex occupation pipeline, energy functional, frozen-H projector oracle, local radial moments, A1/A2/M0/M2/W1 diagnostics): ~28 commits. f95ade1e, 32352799, f0ed469d.
- Reciprocal Fermi-Dirac occupation exposed as a public routine; path variables lengthened to 1024 characters: 2 commits. 65680213, cb1037d5.

### L11 fable_v4b_u (179 commits, LDA+U line)
2026-10-17 .. 2026-11-04. About 90 of 179 commit subjects are "Record/docs" evidence entries.
- Fixed-U energy routed in eligible KS and RS bulk runs; local Hubbard moments; DM restart files (dm sidecars); LDA+U diagnostics extracted: ~12 commits. 8b8f7ef7, 506f9c47, 29bdfa1d.
- RS LDA+U impurity and cluster work (campaign stages C4/C5/C6: impurity residual-mapping out-of-bounds fix, hcp Co cluster cutoff size fix, impurity auto local mode, surface auto local moments, preflight enforcement, zero-U mapping guard): ~60 commits, mostly records. 9fd24bf2, 2686e4e2, a2972286, 7333c275.
- Spinor (nsp=2) d-shell Liechtenstein/FLL functional kernel, spinor occupations, spinor V insertion, spinor local radial moments, routed spinor Hubbard energy; SOC+U bulk by default at nsp=2 (`hubbard_soc_experimental` becomes a no-op): ~30 commits, 2026-10-26/27. 4545a325, d311668b, e7093e2c, aea6366d.
- SOC+U impurities (calctype I) at nsp=2 with host nsp provenance; spin-frame-aware local-moment kernel, `hubbard_axis` provenance, tilted axis needs Chebyshev in RS: ~20 commits. 0e3bfb51, 610f4b72, d0fa6c7e.
- SOC+U surface (calctype S) at nsp=2, Ni(001) decks, frozen-host handling (host must be converged, no dm rewrite), nsp=4 guard test, TESTING_POLICY for U+SOC: ~20 commits. c87851eb, f28186a6, e8a82c77, 66923e3b.
- Fixes found along the way: ldm_beta applied to spinor P (c0cd3825), spinor moments trace uses recursion-site V for impurities (4dfcd374), POISS0 export for Linux link (e809ffb6), `constraints_mom_ref` array read (c1028530), MPI argument/buffer fixes, MPI_BCAST in dm reader (a449c0e3).
- Build/CUDA: native CUDA architecture detection, Chebyshev spectral positivity enforced, mixed MPI buffers in GNU Fortran: 4 commits. dbb5b948, a4b76428, e4c7a474.
- fccNi001 example and regression rows for nsp=1/nsp=2, SOC-to-zero and snapshot comparisons: ~10 commits. d77e8dbd, 43917142, 24837efc.

### L12 fable_v4b_gbt (116 commits, GBT / generalized Bloch theorem line)
2026-10-17 .. 2026-11-04. Restart of GBT in k-space after deleting the old tests.
- Deletion of old GBT source-grep tests, numpy re-implementation tests, validation scripts and fixtures (kept two supercell lattices): 5 commits at start, plus 4 later removals (k-space branches of the RS frozen-magnon routines, little-group mesh code, `q_symmetry_policy` key, spglib little-group wrapper). 9f1ebdf3, ef017163, 93777d87.
- K-space GBT kernels at first and second order (one-band, assembly, Commensurate/Gauge/Antiparallel tests), grand_potential, force-theorem test, `frozen_magnon_kspace` driver and symmetry test, `mode='scf'`: ~20 commits. ad1865a6, 95f9cfe4, c7c1fe1e, d6c1315d.
- Developer-owned GBT oracles O1-O3 (rigid rotation of ASA potential for bcc Fe, band moments O3b): ~5 commits. 31f8a86e, d2f0fb82.
- `kspace_ham_order='orthogonal'` (GBT and also outside GBT: periodic_nc and collinear k-space SCF at q = 0), orthogonal constrained-spiral tests: ~17 commits. 111faed2, 28d30a58, e2518349.
- Constraining field in k-space GBT (local-frame field, Newton update with measured response, density projected on local axis, reference orientation kept in written potential) and constrained-spiral SCF identity tests (bcc Fe, Cr, doubled cell): ~21 commits. 957222ef, d5f4ba49, 39571863.
- CCOR (`ccor_2c`) in k-space GBT (kernel gbt_kspace_hcc, gauge and commensurate identities, NiO AFM fixture; `ccor_vmt_mode pair_surface` refused): ~21 commits. c1b41b75, 33a6810e, 4bdb8d23.
- hoh with ccor_2c in fast Chebyshev and block Lanczos CPU kernels and CUDA plugin, plus exchange/conductivity/orbital cases with legacy references: ~8 commits. d90ef7a8, 010b4c68, 0ab69895, 1b8cf332.
- Hardening: invalid `tetra_symmetry_mode`, `reciprocal_backend`, `green_backend`, `reciprocal_mode` now fatal; `gpu_precision` selects arithmetic of every GPU kernel; cmplx-without-kind precision fix in build_obarm/build_enim; constraints namelist read error fixed (fixed-size buffers, iostat): ~10 commits. 1467633c, f2bc31ce, ab489f52, 5a906684.
- Surface/strux: sbar indexed by type in write_strux_block; finite-cluster strux branch used for surface builds with named crystal_sym: 2 commits. 716411a6, d2b86b88.
- Example `example/k_space_spiral` (bcc Fe, fcc Ni; MFT paths, constrained and unconstrained cone SCF, scf_sweep); Linux RsFreeze references: ~8 commits. b22b9791, cecc48a6, 87787e9c.
- KNOWN_ISSUES entries (S9, self.f90 cancellation, screening tables, SCF stalls, round-off failures): ~23 commits, record only. cb7f563f, 1dacdb8d.

### M main since 9651a900 (34 commits)
2025-10-24 .. 2026-11-05. Authors: Ramon Cardias (most), Anders Bergman (tests, preallocation).
- Graphene example, nml clean-ups: 4 commits (2025-10-24, 2025-11-15). f53cdece.
- Test framework on main ("Update tests" #3/#4/#5): 79 + 31 + 29 files; tests/ grows from 17 files at base to 113 on main (ls-tree count): 3 commits, 2026-03-07..12. 0a6262a8, 09a2d207, c183dadc.
- Pre-allocation before namelist reads (lattice.f90) and dimension/print fixes (charge, exchange, hamiltonian, recursion, self): 3 commits, March 2026. 79f41aea, 847802ed, f9eadd7c.
- Orbital quadrupole operator in bands.f90 (+216), surface cluster calculation tweak, exchange output format, bands print fix: 5 commits, June 2026 (in PR 7 per its stat). eb8966e2, 00835589, 5f426578, 71c97ab4.
- Kubo linear-response rewrite (setup_kubo_operators, `pol_alpha/pol_beta` replace `js_alpha/jl_alpha`, hamiltonian.f90 +1546, recursion.f90 +503), 2026-06-11: 1 commit. 0a7d6843.
- Quadrupole current and accumulation in Kubo operators (2026-07-06), orbital octupole accumulation (#9, 2026-08-31), pure orbital octupoles + conventional current response and quadrupole normalisation/traceless fixes (2026-11-03): 5 commits. f951873b, 49c0031b, c1ffdfd3, 1522cea4.
- Transport rewrite: full HOH actions for all Kubo operators (recursion.f90 +792, new transport_operators.f90), explicit legacy-style loops with HOH corrections, periodic bond displacements in matrix-free Kubo velocities, conductivity normalisation and explicit physical response units, resolved-response filenames by atom label: ~9 commits, Nov 2026. cd0a2432, e692578f, a741176f, a2469d77, dc54ef0c.
- Tight-binding/lsham fixes: `lsham` accepts negative xi, TB orbital-resolved effective field: 2 commits. 495709e7, 22b9de89.
- Build and docs: comprehensive README (e604f402), fccPt conductivity examples (41f3a8ac), Fortran MPI CMake integration and oneAPI LP64 notes (9072004c): 3 commits.
- Benchmarks/regressions added by main: tests/benchmarks/periodic_velocity, tests/benchmarks/response_units, tests/regression/test_{quadrupole,octupole,hoh_transport}_*.py.

---

## 2. Capabilities in fable_v4c (with first appearance)

- LDA+U (U, J, +V intersite, f-orbitals, impurity U, self-consistent U): L00, 2024-06 (a0b5e57f).
- libXC functionals (LDA/GGA, X and C): L00, 2025-09 (72cb73ee); XC semantics reconciled: L08 (91be9db2).
- Gershgorin and Chebyshev spectral bounds: L00 (4f2f629c, e28d94f5).
- K-space (reciprocal) bands, DOS (tetrahedron), symmetry/spglib, k-space SCF: L00, Sep-Oct 2025 (87726190, 52321898).
- Spin-spiral vector `q_ss` input: L00 (32ace7f9).
- strux_lib structure constants (library + default for tests): L01 squash (1355a50c); default test backend L07 (5c0643e7).
- Flexible orbital basis / ldim (basis.f90, orbital-basis dimension tests, k-space orbital basis generalisation): basis.f90 in L01 (1355a50c); tests and k-space generalisation L07 (1d2a061a, ec5cf632). Details of ldim are (inferred) from subjects and stats only.
- Non-collinear k-space (nc) in squash: L01 (inferred from PR title).
- Sparse matrix module and backend drivers: L02 (b95ea55a).
- CUDA recursion kernels, DOS and bgreen on GPU: L02-L04 (bed1378a, ae7981bd); reciprocal CUDA backend, batched eigensolver, KPM on CUDA: L07 (f015b45e, 412c703a, 65ac455b).
- MKL batched / sparse Chebyshev kernels: L03 (cc4bfc20, f5476b32).
- Chebyshev, Block Lanczos, Lanczos recursion backends with fast CPU kernels: Chebyshev fast L02-L03; recursion core split L06 (848ad2f2).
- HOH in new CPU/GPU kernels and reciprocal mode: L04 (1ba2fd4d, 95972953); hoh+ccor_2c in fast kernels and CUDA: L12 (d90ef7a8, 0ab69895).
- CCOR (ccor_2c): L05 (d762e8eb); in k-space GBT L12 (33a6810e).
- GBT / spin spirals: first k-space transform L06 (aabb842d); WP1-10 L07 (fa8d7fc4); restart in k-space with orthogonal ordering, constrained SCF, oracles L12 (ad1865a6, 111faed2, 957222ef).
- Frozen-magnon post-processing: L06 (cbbb3208); k-space driver L12 (5fcb53c7).
- k-space Green functions (Lehmann/Dyson), Bloch spectral functions, route-agnostic moments/J_ij/sigma: L06 (def7d8b6, b1c9b4ff, 52374509, 75aed0b1).
- Gilbert damping and magnetic inertia: L06 (22b42c85); damping/inertia split out of exchange L07 (c9dcb1d4).
- Surface/interface/vacuum electrostatics, `calctype='L'`: L06 (8c28ef62, 97f1e0e4, 8883c608).
- Constraining fields: pre-existing in squash/abspinlib (constrain.f90 stat in 1355a50c); k-space GBT field L12 (957222ef); constraints_mom_ref read fix L12 (5a906684).
- Spin response `spin_response` (Juelich-d, channels s,p,d,f, multi-site, AFM, force-theorem): L09 (1cc43b20, 6fe83577, 02582da1).
- Fixed-U LDA+U with local moments and routed energy (KS and RS): L09 (f95ade1e) and L11 (8b8f7ef7, 506f9c47).
- SOC+U at nsp=2 for bulk, impurity, surface: L11 (aea6366d, 0e3bfb51, c87851eb).
- Transport / Kubo-Bastin with spin and orbital operators, KPM transport: L06 `fix: restore conductivity postproc outputs` (48c18c29), L07 KPM CUDA (5a1ace1a), recursion_transport.f90 present in v4c file list.
- Orbital quadrupole operator in bands.f90: via PR 7, L06 (a8ec82e9).
- Test/CI infrastructure (manifest-driven tests, quick tier, reference workflow, CI matrix, GPU manual matrix): squash L01, L06 P1-P10, L07, L09 quick tier (32772694).
- Docs: Sphinx tree L01; DEVELOPER_MAP/CLAUDE.md L06 (eced2068); ROADMAP/DECISIONS L08 (12c89fae); TESTING_POLICY L11 (e8a82c77).
- Removed in fable_v4c lineage: linear-response/TDDFT campaign code (L09, 63202725; archive tag lr-campaign-archive-2026-10); old GBT tests (L12, 9f1ebdf3); legacy Lanczos/Block/Chebyshev harness (L09, 06407dad); `exchange_q`, `linear_response`, `pauli_projection` post_processing values (stop at startup; DEVELOPER_MAP).

## 3. Capabilities only on main (not in fable_v4c)

Main commits dated after 2026-06-15 (not in fable_v4c):
- Quadrupole current and accumulation in Kubo operators (f951873b, 2026-07-06).
- `lsham` build accepting negative xi (495709e7, 2026-08-21); TB orbital-resolved effective field fix (22b9de89, 2026-08-24).
- Orbital octupole accumulation operators (49c0031b, 2026-08-31); pure orbital octupoles and conventional current response (1522cea4, 2026-11-03).
- Quadrupole normalisation and traceless diagonal components fix (c1ffdfd3, 2026-11-03).
- Full HOH actions for all Kubo response operators; transport_operators.f90 (cd0a2432, c9165744, e692578f, 2026-11-03).
- Periodic bond displacements in matrix-free Kubo velocities, with periodic_velocity benchmark (a741176f, 2026-11-04).
- Conductivity normalisation and explicit physical response units, atomic-label filenames (a2469d77, dbf81309, dc54ef0c, 2026-11-04/05).
- fccPt conductivity examples updated (41f3a8ac); comprehensive README (e604f402); Fortran MPI CMake integration and oneAPI LP64 build notes (9072004c, 2026-11-05).

Main commits dated up to 2026-06-15 that are not visible in fable_v4c by the checks above:
- Kubo linear-response rewrite 0a7d6843 (2026-06-11; `pol_alpha`/`setup_kubo_operators` absent from fable_v4c source). Reason is (inferred): PR 7 stat does not include hamiltonian.f90/recursion.f90.
- Pre-allocation edits in lattice.f90 (79f41aea, 847802ed) and the f9eadd7c dimension/print edits: whether they are in fable was not checked; belem25 squash rewrote lattice.f90 (+2071 lines changed), so overlap is likely (inferred).
- Graphene example: fable_v4c has 4 paths containing "graphene" in its tree (ls-tree grep), so it appears to be present; not diffed.

## 4. Possible overlaps between main and fable

Files both sides changed since base 9651a900 (from stats; textual conflict is (inferred) where noted):
- `source/hamiltonian.f90`: main rewrote Kubo operator setup, quadrupole/octupole operators, lsham, HOH Kubo actions (0a7d6843, f951873b, 495709e7, 49c0031b, c1ffdfd3, cd0a2432, c9165744). Fable rewrote it (1355a50c, +2597 lines changed) and split it into `hamiltonian_build.f90`, `hamiltonian_ccor.f90`, `hamiltonian_hubbard.f90`, `hamiltonian_paoflow_io.f90` (L06/L07, file list in fable_v4c). Main's edits will not apply textually (inferred).
- `source/recursion.f90`: main added Kubo/HOH transport (+792 lines in cd0a2432, +743 in e692578f); fable split recursion into `recursion_core`, `recursion_chebyshev`, `recursion_haydock`, `recursion_transport` and added fast/batched/MKL/CUDA kernels (L02-L06). Highest overlap with fable's transport path.
- Transport / Kubo / conductivity: main = Kubo operators incl. quadrupole/octupole, physical units, periodic velocity (`source/conductivity.f90` +310 net vs base); fable = KPM transport and Kubo-Bastin on CUDA/CPU (L07: 5a1ace1a, 348e8a2a, 434ab35e "Update transport calculation and conductivity paths"), `restore conductivity postproc outputs` (48c18c29), `&reciprocal` moments (75aed0b1). Both changed `conductivity.f90`, namelists `conductivity.f90`.
- Namelists: `namelists/hamiltonian.f90` (main: pol_alpha/pol_beta replace js_alpha/jl_alpha; fable: q_ss, theta_ss, magnetic_representation, ccor_2c, Hubbard keys), `namelists/control.f90` (both, incl. PR 7), `namelists/self.f90` (main), `namelists/conductivity.f90` (main).
- `source/bands.f90`: main quadrupole (+216, 2026-06-11, in PR 7), print fix (71c97ab4) and TB effective field (22b9de89, +99 later); fable rewrote bands.f90 in the squash (+2130 lines changed). The later TB orbital field fix is not in fable.
- Quadrupole/octupole operators: main only; the 2026-06 quadrupole piece reached fable_v2 via PR 7, octupole and quadrupole-current did not.
- Tests framework: main "Update tests" #3-#5 (tests/ 17 -> 113 files, run via its own scripts) vs. fable tests/ harness (run_test.py, tests/scf/cases.json, manifests, quick/full tiers; 821 files in v4c). Both touch `.github/workflows/tests.yml` (79f41aea on main; squash and L06/L07 on fable). Main's added tests/regression/test_*.py and tests/benchmarks/* sit in a tree that fable reorganised (inferred).
- Pre-allocation / ldim / dimensions: main 79f41aea, 847802ed, f9eadd7c (lattice.f90, charge.f90, exchange.f90, hamiltonian.f90, recursion.f90, self.f90) vs. fable basis.f90 and flexible basis tests (L01, L07), "Corrected atomic namelist reading" (e1e8c7a6), path length fix (cb1037d5). Same concern: dimension handling and namelist reads.
- `source/lattice.f90`: main pre-allocation and surface cluster tweak (00835589), periodic bond displacements (a741176f); fable squash rewrote it (+2071 lines changed), added strux_lib calls, finite-cluster strux for surfaces (d2b86b88), sbar indexing (716411a6).
- `source/exchange.f90`: main exchange output format (5f426578, in PR 7); fable exchange energy integration to enpt+1 (2cb6f554), damping/inertia split (c9dcb1d4), temperature (7be98f7a).
- `source/self.f90` (off-limits per CLAUDE.md for new edits): main touched it in 2e1aabef, f9eadd7c, 22b9de89; fable touched it in the squash and L11/L12 notes. Any merge needs the owner's attention.
- `source/charge.f90`: main f9eadd7c (6 lines); fable split interface electrostatics from charge (b559b2c3).
- `CMakeLists.txt` / MPI: main 9072004c (Fortran MPI CMake integration, oneAPI LP64) vs fable CMake changes (L02 build modes, L06 T13 hygiene, L11 CUDA architecture detection dbb5b948, MPI fixes 6b0bf345, e4c7a474, MPI CI 99d28cd5).
- Docs: main README rewrite (e604f402) vs fable Sphinx docs/, DEVELOPER_MAP, ROADMAP and CLAUDE.md; main added docs/octupole_response.md.
- Examples: `example/conductivity/fccPt` (main 41f3a8ac), fable reorganised examples and deleted some decks (L09 dd3d8ca1); not diffed.
- Not an overlap but related: fable deleted linear-response/TDDFT (L09) and old GBT tests (L12); main has no equivalent work, so nothing from main depends on them (inferred).
