# TDDFT Ward-repair tracker

- [x] WR-00 — Freeze pair-potential conventions and executable identities
- [x] WR-01 — Add generic weighted-vertex and direct-Xi machinery
- [x] WR-01b — Preserve endpoint-resolved LMTO magnetic tangents
- [ ] WR-02 — Implement the LMTO pair-potential operator with rotation oracle
- [ ] WR-03 — Integrate direct Xi in production shadow mode
- [ ] WR-04 — Implement the real static limit and ground-state provenance
- [ ] WR-05 — Repair Goldstone option semantics and controlled correction
- [ ] WR-06 — Repair Cartesian/circular factors and unsupported-route gates
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

- [ ] Consume endpoint Cartesian tangents with endpoint-position Bloch phases,
  establish signed `M_a` provenance, form `Q_a^-`, and validate against a
  commensurate-supercell finite-rotation oracle. No assembled-`hxc` or scalar
  site-kernel reconstruction is allowed.

## Global gates

- [ ] Raw no-SOC `Xi(0,0) M = M` converges before correction
- [ ] Raw no-SOC Gamma pole is below 1 meV for converged Fe and Ni runs
- [ ] Static result is stable under the defined k-mesh, temperature, band, basis, and static-limit sweeps
- [ ] Small-q stiffness is stable and compared independently with corrected GBT/frozen-magnon and `Jij`
- [ ] Bare Stoner, enhanced, Xi-eigenvalue, and mode-projected outputs agree on the interpretation
- [ ] Generalized-overlap response is either metric-correct and tested or explicitly unsupported
- [ ] No GPU work has been started before all preceding CPU physics gates pass
