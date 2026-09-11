# TD-DFT clean-room purge manifest

Status: LR-00 inventory recorded before purge

Starting commit: `012bdf2be1ca9d1a79d53a04d88f53d5b589b687`

Archive tag: `pre_tddft_cleanroom_20260909`

The active response implementation at the starting commit is archived by Git
and is not retained as a production fallback.  This manifest classifies the
files and routines that were inspected before the LR-00 removal.

## DELETE_PHYSICS

These files encode the old site/scalar, pair-Xi, response-vertex, Dyson,
Goldstone, mode, longitudinal, non-collinear, or TD-DFT-specific response
contracts.  They are removed from the active build and source tree.

### Response implementation

- `source/tddft_backend.f90` — backend polymorphism and backend selection.
- `source/tddft_chi0.f90` — eigenpair Kohn–Sham response and static response.
- `source/tddft_chi0_green.f90` — Green-function response adapter.
- `source/tddft_chi0_realspace.f90` — native real-space response backend.
- `source/tddft_circular.f90` — circular-channel conventions and routing.
- `source/tddft_config.f90` — `&tddft` parser and validation.
- `source/tddft_conventions.f90` — old response factors and denominators.
- `source/tddft_dyson.f90` — legacy susceptibility layer was removed; the
  active file is now the clean TDDY-01 canonical Dyson/loss service.
- `source/tddft_four_component.f90` — old four-component response path.
- `source/tddft_goldstone.f90` — old Goldstone diagnostics and corrections.
- `source/tddft_longitudinal.f90` — old charge/longitudinal response.
- `source/tddft_modes.f90` — old loss/mode analysis.
- `source/tddft_occupation.f90` — TD-DFT occupation state.
- `source/tddft_performance.f90` — TD-DFT work decomposition/profile objects.
- `source/tddft_transition_engine.f90` — old transition and pair-operator engine.
- `source/tddft_ward.f90` — old Ward/sum-rule/projection diagnostics.
- `source/tddft_xi.f90` — direct old Xi construction.
- `source/response_basis.f90` — old site response basis.
- `source/response_components.f90` — old response-channel labels.
- `source/response_vertices.f90` — old site/projected Pauli vertices.
- `source/xc_response_kernel.f90` — old site-projected ALSDA kernel provider,
  circular kernel, longitudinal derivative, and Goldstone kernel data.

### Pair-potential and response tangent implementation

- `source/lmto_pair_potential.f90` — pair-potential Xi, endpoint phases, and
  response spinor unfolding.
- TD-DFT-only routines removed from `source/lmto_magnetic_tangent.f90`:
  `lmto_bond_tangent`, `lmto_make_endpoint_record`, and
  `lmto_ordinary_tangent_supported`.
- TD-DFT-only routines removed from `source/hamiltonian_build.f90`:
  `ham0m_nc_tangent`, `ham0m_nc_endpoint_tangents`,
  `lmto_tangent_capability`, and `lmto_has_active_soc`.
- TD-DFT-only routine removed from `source/reciprocal_fourier.f90`:
  `build_lmto_pair_potential_at_kpoint`.
- TD-DFT-only reciprocal binding removed from `source/reciprocal.f90`:
  `build_lmto_pair_potential_at_kpoint` and its pair-transition metadata.
- TD-DFT-only `lmto_pair_operator_tile_source`, pair-Xi builders, pair
  diagnostics, metadata writers, and susceptibility helpers removed from
  `source/calculation.f90`.

### TD-DFT-specific input, tests, and validation artifacts

- `source/include_codes/namelists/tddft.f90` — the legacy parser was removed;
  TDRUN-01 now uses this path for the minimal clean production namelist.
- `source/self_xc_response.f90` — removed with the TD-DFT response-provider
  synchronization hooks.
- `tests/unit/test_tddft_*.f90` and `tests/unit/test_tddft_*.py` — expected
  behavior of the removed response layer.
- `tests/unit/test_response_conventions.f90` and
  `tests/unit/test_response_vertices.f90` — old response convention/vertex
  behavior.
- `tests/unit/test_lmto_pair_potential.f90` and
  `tests/unit/test_lmto_magnetic_tangents.f90` — pair-Xi/tangent behavior.
- `tests/unit/test_magnetic_scf_feedback.f90` — response-provider projection
  behavior removed with the provider; the ordinary `tests/magnetic_scf/`
  campaign is retained separately.
- `tests/regression/tddft_validation/` — old TD-DFT validation driver,
  fixtures, expected response outputs, and campaign data.
- `example/susceptibility/` — old active TD-DFT input/output examples.
- Old TD-DFT validation/user-facing documents under `docs/` and `docs/dev/`
  are removed where they claim that the purged response implementation is
  available or validated.  The literature-locked cleanroom prompt package is
  retained as the governing redevelopment record.

## RETAIN_GENERIC

These components do not encode TD-DFT response physics and remain available
to ordinary SCF, band, DOS, transport, exchange, and Green-function paths.

- `source/dyson_kernel.f90` — generic dense resolvent inversion primitive used
  by reciprocal Green-function and band-structure consumers; it is not the
  removed susceptibility Dyson layer.
- `source/lehmann_kernel.f90` — generic resolvent/Lehmann matrix algebra.
- `source/green*.f90` — ordinary real-space Green functions and recursion
  infrastructure.
- `source/reciprocal*.f90` — reciprocal mesh, Hamiltonian Fourier, eigenpair,
  DOS, band, and Green infrastructure after pair-Xi hooks are removed.
- `source/hamiltonian*.f90` — ordinary LMTO Hamiltonian construction,
  including the value map used by SCF; TD-DFT tangent entry points are removed.
- `source/lmto_magnetic_tangent.f90` — retained only for the ordinary magnetic
  bond-value algebra and spinor conversion used to assemble the ground-state
  Hamiltonian and GBT path; no response tangent is retained.
- `source/math.f90`, `source/array.f90`, `source/sparse.f90`, `source/safe_alloc.f90`,
  `source/mpi.f90`, `source/kpoint_workset.f90`, and Fourier helpers —
  physics-neutral numerical infrastructure.
- `source/spin_dynamics.f90`, `source/exchange*.f90`, `source/conductivity.f90`,
  and `source/frozen_magnon.f90` — independent non-TD-DFT consumers.

## RETAIN_GROUND_STATE

- `source/self.f90` — retained SCF radial Poisson/XC evaluation and potential
  construction.  The TD-DFT-only XC response-provider fields, projection
  accumulation, longitudinal derivative requests, and response diagnostics
  are removed.
- `source/xc.f90` and `source/xc_radial.f90` — retained ground-state XC
  functional evaluation and radial XC support.
- `source/potential.f90`, `source/charge*.f90`, `source/symbolic_atom.f90`,
  and SCF preprocessing modules — retained ground-state state and radial data.

## REVIEW_DEPENDENCY

These files contain ordinary paths plus TD-DFT dispatch or dependency edges;
the listed TD-DFT portions are removed while unrelated paths are retained.

- `source/calculation.f90` — remove the old susceptibility dispatch and
  response helper types/routines; retain only the TDRUN-01 minimal production
  input hook.
- `source/calculation_reciprocal.f90` — remove native real-space response-pair
  setup; retain the shared ordinary post-processing stack.
- `source/hamiltonian.f90` — remove tangent interfaces and pair-response type
  dependencies; retain ordinary `ham0m_nc`.
- `source/reciprocal.f90` and `source/reciprocal_fourier.f90` — remove
  pair-potential construction; retain normal Fourier/eigenpair operations.
- `source/CMakeLists.txt` and top-level `CMakeLists.txt` — remove deleted
  sources and TD-DFT tests from the build/test graph.
- `docs/DEVELOPER_MAP.md` and other cross-references — remove or revise stale
  claims that the purged TD-DFT path is active.
- `tests/KNOWN_ISSUES.md`, test indexes, and validation indexes — remove stale
  active TD-DFT validation claims while retaining unrelated SCF/Green/GBT
  evidence.

## User-facing clean-room contract

Any `post_processing='susceptibility'` request, or an unsupported/legacy
`&tddft` form, fails explicitly.  The clean minimal `&tddft` form is accepted
only with `post_processing='tddft'`.  The old route fails with:

`post_processing='susceptibility' is the removed legacy TD-DFT route; use post_processing='tddft' with the minimal &tddft input.`

No old response equation, sign, factor, basis, radial approximation, kernel,
Goldstone correction, or mode-analysis implementation is replaced or inferred
by LR-00.
