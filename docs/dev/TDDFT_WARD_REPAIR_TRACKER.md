# TDDFT Ward-repair tracker

- [x] WR-00 — Freeze pair-potential conventions and executable identities
- [x] WR-01 — Add generic weighted-vertex and direct-Xi machinery
- [x] WR-01b — Preserve endpoint-resolved LMTO magnetic tangents
- [x] WR-01c — Validate finite-q endpoint phases and pair-potential gauge
- [x] WR-02 — Implement the LMTO pair-potential operator with rotation oracle
- [x] WR-03 — Integrate direct Xi in production shadow mode
- [x] WR-04 — Implement the real static limit and ground-state provenance
- [x] WR-05 — Repair Goldstone option semantics and controlled correction
- [x] WR-06 — Repair Cartesian/circular factors and unsupported-route gates
- [ ] WR-07 — Replace heuristic mode classification
- [ ] WR-08 — Pass Fe/Ni material gates and switch the validated default

## WR-00 evidence

- `UnitTddftWardConventions`: analytic uniform and unequal-orbital
  pair-potential oracles, including the scalar-kernel negative control.
- Existing `UnitResponseConventions`, `UnitResponseVertices`, `UnitTddftChiKS`,
  and `UnitTddftGoldstone` remain required regression evidence.

## WR-01 evidence

- `UnitTddftDirectXi`: uniform scalar reduction, unequal-orbital weighted
  oracle, complex nonzero-`q` scalar/batched equivalence, metadata propagation,
  and wrong-right-vertex-order negative control.
- Regression set: `UnitResponseConventions`, `UnitTddftWardConventions`,
  `UnitResponseVertices`, `UnitTddftChiKS`, `UnitTddftGoldstone`,
  `UnitTddftDysonModes`, `UnitTddftConfig`, and `UnitKspaceOccupations`.

## WR-01b evidence

- `UnitLmtoMagneticTangents`: old/reference value preservation, complex spinor
  map, three-step central finite differences for either endpoint and a common
  rotation, dot/cross negative controls, reverse-bond covariance, same-type
  distinct-site identity, and capability rejection.
- Focused regression: `UnitGbtStructure`, `UnitResponseConventions`,
  `UnitTddftWardConventions`, and `UnitTddftDirectXi`.

## WR-02 handoff

- `UnitLmtoPairPotential`: q=0 LMTO circular pair-potential construction,
  full-normal-builder finite rotation, unequal-orbital scalar negative
  control, signed-moment normalization, adjoint relation, normal Bloch phase,
  rigid spinor covariance, and rejected zero moment.
- The q=0 rotation oracle rotates the two response moments, rebuilds every
  directed block through production `ham0m_nc`, and compares its Cartesian
  central difference with the production reciprocal Q service.  A consistent
  rigid rotation also satisfies `H(theta)=U_y(theta) H(0) U_y(theta)^dagger`.
- The service remains restricted to `ham_only`; no assembled-`hxc` or scalar
  site-kernel substitute is used.

## WR-01c closure evidence

- `UnitLmtoPairPotential`: endpoint phase isolation at `q=1/2` and `q=1/3`,
  q=0 reduction, and unfolded/folded two-sublattice site-gauge metadata.
- `UnitLmtoPairPotential`: independent explicit 2- and 3-cell cosine/sine
  texture oracle. It rebuilds the ordinary LMTO bond Hamiltonian from rotated
  moments, projects the central difference from `k=0` to `k+q`, and verifies
  the analytic circular Q operator over a three-theta O(theta^2) ladder.
- The reciprocal service now assembles `Qplus` from a separate reverse
  transition accumulator; `UnitLmtoPairPotential` checks its adjoint identity
  from distinct forward/reverse inputs.
- `UnitLmtoPairPotential` now drives the production reciprocal service on a
  completed ordinary-LMTO fixture.  It verifies q=0 reduction, q=1/2 and
  q=1/3 analytic endpoint-phase operators, reverse Qplus, and named
  zero-moment/overlap/HOH/GBT/CCOR/SOC rejection paths.
- The same service fixture now includes two response sites with identical
  potential parameters but distinct site identities and absolute positions;
  its onsite blocks prove that perturbing one site does not leak onto the
  other.
- The closure oracle now directly compares the actual reciprocal service with
  an independently rebuilt two-sublattice N-cell normal Hamiltonian at q=1/2
  and q=1/3, for each of the two response-site channels.  It uses normalized
  cosine/sine textures, absolute-position projection, three theta values, and
  production `ham0m_nc` only; it cannot consume endpoint tangents or the
  analytic Q assembler.  The service `Qplus`, assembled through its separate
  reverse accumulator, agrees with the supercell `Qminus` adjoint.
- gfortran-13 evidence (2026-08-10): final two-sublattice absolute errors are
  2.973e-10 (q=1/2, site 1), 1.486e-10 (q=1/2, site 2), and 2.229e-10 for
  each q=1/3 site; every theta ladder shows the expected O(theta^2) reduction.
  `UnitLmtoPairPotential` reports a maximum error of 8.856e-10.
- Implementation evidence has been accepted by the designated TDDFT/LMTO
  physics maintainer; WR-01c and WR-02 are closed.

## WR-03 evidence

- `xi_backend='legacy_site_scalar'` remains the default comparison route and
  retains its existing `*_dyson.dat` filename.  `pair_potential` writes
  `*_pair_dyson.dat`; `compare` writes separate `*_legacy_dyson.dat` and
  `*_pair_dyson.dat`, plus correspondingly named mode files.
- The pair route builds `Q_a^-(k,q)` for every response site and k point from
  the validated ordinary-LMTO `ham_only` rotation service, accumulates direct
  Xi from the same eigenpairs/occupations/k weights/band window as chiKS, and
  solves `(I-Xi) chi=chiKS` without constructing a response-space kernel.
  Generalized-overlap, Green, longitudinal, and full-component pair routes
  are explicitly rejected pending their own metric/derivative work.
- `*_goldstone.dat` retains the legacy raw diagnostic and, for pair/compare,
  adds separately named legacy and pair raw residuals with pair-potential
  representation, provenance, and signed-moment source.  No correction is
  applied by WR-03.
- gfortran-13 build and focused evidence (2026-08-10):
  `cmake --build build-gfortran13 --target rslmto UnitTddftDirectXi UnitTddftDysonModes UnitTddftConfig -j2`;
  `ctest --test-dir build-gfortran13 --output-on-failure -L tddft` (16/16
  passed); `python3 tests/unit/test_tddft_dispatch.py` and
  `python3 tests/regression/tddft_validation/test_validation.py` (both pass).
  `UnitTddftDirectXi` includes the k-resolved-Q reduction and nonuniform-k
  retention checks; `UnitTddftDysonModes` checks direct-Xi versus kernel Dyson
  equivalence in the uniform oracle.  No MPI-enabled build was available in
  this workspace, so serial-versus-multirank equivalence remains unchecked.

## WR-04 evidence

- `build_static_chi_ks_from_eigenpairs` and the k-resolved
  `build_static_direct_xi_from_k_dependent_eigenpairs` use the real q=0,
  omega=0 Lehmann divided difference.  Exact and near degeneracies use the
  analytic derivative of the same Fermi function; the static implementation
  has no dynamic-eta input.
- The production driver inherits reciprocal ground-state Fermi level,
  temperature, and electron count unless an `&tddft` value is explicit.  It
  records both provenance values and override flags, recomputes the response
  count from the response eigenpairs, and fails when `|dN| > 1e-8 max(1,N)`.
  Its transverse response moment is the signed occupied `P_site sigma_z`
  population.
- At Gamma, `*_goldstone.dat` additionally records observed raw loss-grid
  maxima for bare and available legacy/pair enhanced spectra.  These dynamic
  measurements are explicitly separate from the real static Ward operator
  and are not fitted shifts or corrections.
- gfortran-13 evidence (2026-08-10):
  `cmake --build build-gfortran13 --target UnitTddftChiKS UnitTddftDirectXi UnitTddftGoldstone UnitTddftConfig -j4`;
  `ctest --test-dir build-gfortran13 --output-on-failure -R 'UnitTddft(ChiKS|DirectXi|Goldstone|Config)'`
  (4/4 passed).  The focused tests cover exact/near-degenerate divided
  differences, zero-temperature and finite-temperature limits, static
  eta-independence, the two-orbital pair-potential Ward identity, reversed
  one-site and signed ferro/antiferromagnetic two-site moments, and explicit
  provenance overrides. `python3 tests/unit/test_tddft_dispatch.py` confirms
  the material electron-count mismatch is a fatal production path.
- No material Gamma pole, raw material Ward convergence, MPI comparison, or
  eta-ladder material evidence is claimed here; those global gates remain
  open for WR-08.

## WR-05 evidence

- `goldstone_mode='off'` constructs no Goldstone diagnostic or correction;
  `diagnose` writes raw diagnostics and uses raw direct pair-potential Xi;
  `correct` requires `xi_backend='pair_potential'` or `compare`.  Deprecated
  `sum_rule` inputs migrate to `correct` with an explicit warning and output
  provenance flag.
- `build_goldstone_column_correction` solves the real-static constrained
  column problem `Re[Xi_static diag(M)] s=M` with a rank-revealing SVD and
  unit minimum-change weights. It rejects material static imaginary content,
  rank deficiency, small moments, nonfinite scales, and any scale change over
  25 percent. No finite-dynamic-eta inverse is used.
- Raw `*_pair_dyson.dat`, corrected `*_pair_corrected_dyson.dat`, static raw
  and corrected residuals, SVD rank/condition, every scale, decision, Gamma
  loss maxima, and the corrected spectral-weight nonnegativity check are
  emitted independently. A rejected correction writes raw diagnostics then
  terminates clearly; it never falls through to diagnose.
- gfortran-13 evidence (2026-08-10):
  `cmake --build build-gfortran13 --target UnitTddftGoldstone UnitTddftConfig -j4`;
  `ctest --test-dir build-gfortran13 --output-on-failure -L tddft` (16/16
  passed); `python3 tests/unit/test_tddft_dispatch.py` and
  `python3 tests/regression/tddft_validation/test_validation.py` (both pass).
  `UnitTddftGoldstone` covers accepted one-/multi-site scales, dynamic-Xi
  change, rank-deficient/ill-conditioned, small-moment, complex-static, and
  negative-spectrum controls; `UnitTddftConfig` covers explicit and migrated
  option semantics.

## WR-06 evidence

- `k_perp_circular` and `circular_transverse_kernel` are reserved for the
  unhalved circular vertices.  `cartesian_transverse_kernel` is the separately
  named `2*k_perp_circular` derivative used by the full Cartesian x/y block;
  the old direct reuse of the circular scalar is removed.
- `UnitTddftWardConventions` proves bare and dressed Cartesian/circular
  equivalence in the collinear no-SOC two-level oracle and includes the
  wrong-Cartesian-half negative control.  `UnitTddftFourComponent` verifies
  the Cartesian factor and the complete selected-XC derivative capability
  gate.  Unvalidated or incomplete derivative slots are rejected before a
  production `channel='full'` response starts.
- Production metadata labels `chi0_backend='green'` as the
  eigenpair-resolvent reference and states that a native RS Green provider is
  unavailable.  The longitudinal production route now fails explicitly until
  its static calibration uses the WR-04 real static limit; no LLB readiness is
  claimed.  GPU work remains deferred.
- gfortran-13 evidence (2026-08-10):
  `cmake --build build-gfortran13 --target UnitTddftWardConventions
  UnitTddftFourComponent UnitResponseConventions UnitTddftGoldstone
  UnitTddftConfig -j4`; `ctest --test-dir build-gfortran13
  --output-on-failure -L tddft` (16/16 passed); `python3 tests/unit/test_tddft_dispatch.py` and
  `python3 tests/regression/tddft_validation/test_validation.py` (both pass).

## Global gates

- [ ] Raw no-SOC `Xi(0,0) M = M` converges before correction
- [ ] Raw no-SOC Gamma pole is below 1 meV for converged Fe and Ni runs
- [ ] Static result is stable under the defined k-mesh, temperature, band, basis, and static-limit sweeps
- [ ] Small-q stiffness is stable and compared independently with corrected GBT/frozen-magnon and `Jij`
- [ ] Bare Stoner, enhanced, Xi-eigenvalue, and mode-projected outputs agree on the interpretation
- [ ] Generalized-overlap response is either metric-correct and tested or explicitly unsupported
- [ ] No GPU work has been started before all preceding CPU physics gates pass
