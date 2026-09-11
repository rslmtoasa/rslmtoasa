# TDRUN-01 — clean post-SCF TD-DFT production driver

Status: implemented for the validated reciprocal collinear baseline.

This document describes the production seam added after the LR-00 purge.  It
does not add response equations or a second TD-DFT implementation.  The
driver assembles requests for the rebuilt response services and owns the
accepted-state lifecycle, capability gate, and result provenance.

## Preflight and input contract

LR-00 removed the old `tddft_backend`, `tddft_chi0*`, response-object,
transition-engine, mode, Goldstone, and susceptibility-dispatch paths.  The
surviving generic calculation namelists were retained, but an old parser
remnant was not used as a reason to resurrect the old object graph.

TDRUN-01 adds the minimal parser at
[`source/include_codes/namelists/tddft.f90`](../source/include_codes/namelists/tddft.f90)
and the production config/driver at
[`source/tddft_production_driver.f90`](../source/tddft_production_driver.f90).
The active calculation contract is:

```text
&calculation
  pre_processing = 'bravais'
  post_processing = 'tddft'
/
&tddft
  enabled = .true.
  channel = 'chi_plus'                 ! chi_plus or chi_minus
  n_q = 1
  q_list = 0.0, 0.0, 0.0               ! reduced coordinates, repeated by q
  n_omega = 3
  omega_min = 0.00
  omega_max = 0.20
  eta = 0.01
  response_lmax = -1                   ! default: complete 2*lmax product
  interaction_route = 'direct_alsda'
  goldstone_correction = .false.
  backend = 'lehmann'
  write_full_matrix = .true.
  output_file = 'tddft_response.dat'
/
```

`use_omega_grid=.true.` with `omega_grid(:)` is also supported.  The optional
`backend='reciprocal_gf'` selects the validated reciprocal-GF susceptibility
service; it does not select the native real-space RSGF route.  The initial
interaction route is `direct_alsda`.  `goldstone_sumrule` remains an explicit
service route, while the `goldstone_correction` switch is rejected until its
separate TDVAL evidence is complete; it is never silently applied.

The old `post_processing='susceptibility'` spelling is rejected with a
migration error.  `&tddft` is feature-off when absent.  Ordinary calculations
therefore do not enter the driver.

## Accepted ground-state handoff

The call is made only for the accepted bulk SCF object graph.  The driver
consumes:

- the accepted Hamiltonian/potential used by SCF;
- the LR-01 radial snapshots in
  `lattice%symbolic_atoms(nbulk+1:nbulk+nrec)%radial_ground_state`;
- XC functional, backend, TXC, density convention, and accepted radial mesh
  provenance from those snapshots;
- the accepted `energy%fermi` and temperature/smearing;
- a dedicated reciprocal eigenpair service bound to the accepted Hamiltonian;
  and
- the accepted lattice, basis, and structural metadata.

It requires `valid`, `accepted`, and `pauli_basis_valid` LR-01 snapshots.  It
does not reconstruct a second ground state, invoke SCF a second time, or find a
different Fermi level.  Reciprocal eigenpairs are generated from the accepted
potential in the dedicated reciprocal service cache, with the accepted Fermi
level fixed before occupations are snapshotted.  The accepted radial state is
passed to response services through immutable references; the driver copies
only the Pauli radial arrays into the response-basis view required by the
service API.

The response-space origin is an exact zero-measure coordinate.  The rebuilt
vertex service already supplies its regular s-s origin limit and zeroes
non-spherical origin products; no epsilon or origin regularizer is introduced
by the driver.

## Lifecycle and orchestration boundary

The bravais path now has this ordering:

```text
input parse
    -> capability gate
    -> SCF / ATOMSC
    -> accepted LR-01 snapshot
    -> reciprocal eigenpairs and exact k+q endpoint snapshots
    -> TDRUN-01 response request batch
    -> LR-05 vertex -> LR-06/LR-GF-02 chiKS
    -> KXC-01 or explicit GSR-01 interaction
    -> TDDY-01 Dyson and loss
    -> provenance/result output
    -> ordinary report/save cleanup
```

The hook is immediately after `self_obj%run()` in
[`source/calculation_preprocessing.f90`](../source/calculation_preprocessing.f90)
and before the SCF-owned objects are released.  Driver input objects are
`intent(in)` except for the dedicated reciprocal service cache.  No SCF object
is modified by response evaluation.

## Capability gate

The gate runs before reciprocal mesh construction and response work.  The
validated baseline is:

| feature | accepted baseline |
| --- | --- |
| spin | collinear (`nsp=1`) |
| SOC | absent |
| reciprocal representation | `ham_only` |
| overlap | orthogonal; no generalized overlap |
| Hamiltonian | second order |
| operators | no unsupported additive operators |
| basis | `sp` or `spd` (`lmax=1` or `2`) |
| geometry | bulk/bravais |

Each rejection names the feature, for example `spin-orbit coupling` or
`generalized overlap / non-orthogonal representation`.  This gate is backed
by the live control, Hamiltonian, lattice, and reciprocal objects rather than
by a caller-provided claim.

## Initial response route

The driver invokes existing services only:

1. `evaluate_pauli_transition_vertex` through the selected susceptibility
   service;
2. `evaluate_lr_ks_susceptibility` for spectral/Lehmann `chiKS`, or the
   explicitly selected `evaluate_lr_gf_susceptibility` reciprocal-GF route;
3. `evaluate_lr_alsda_kernel` for the direct interaction, or the explicit
   `evaluate_lr_goldstone_sumrule` route when selected; and
4. `evaluate_tddft_dyson` for enhanced susceptibility and loss output.

No susceptibility, augmentation, Kxc, Goldstone, or mode-extraction formula
is present in the calculation layer or in the driver.  The native RSGF
backend is not exposed by this task; it remains gated until RSGF-01 passes.

## Output provenance

The default output retains the complete canonical matrices.  For a large
response space, `write_full_matrix=.false.` writes the metric trace per q and
frequency while the complete matrices remain in the driver result object;
this is the mode used by the small integration smoke fixture.

Every output contains the following metadata before any result row:

- Git/version field (or an explicit unavailable marker);
- structure/calctype, lattice constant, atom-type/site count, and symbol;
- XC identity, backend, and TXC;
- accepted magnetic moment;
- accepted Fermi level and temperature/smearing;
- response capability;
- q list and omega list;
- eta;
- response angular cutoff;
- radial mesh identity (`n`, `a`, `b`, `rmax`, and mesh sum);
- interaction route and service provenance;
- Goldstone-correction status; and
- backend.

The q block records the actual reduced q coordinates, and each result row
records the actual omega.  A result is not emitted without these fields.

## Integration evidence

The TDRUN-01 tests are registered in the top-level CTest graph:

- `UnitTddftProductionDriver`: absent-input feature-off behavior, direct
  service-versus-driver reproducibility, and accepted-radial-state
  non-mutation;
- `UnitTddftProductionDriverRejectSoc`: expected SOC failure;
- `UnitTddftProductionDriverRejectGeneralizedOverlap`: expected generalized
  overlap failure; and
- `TddftProductionDriverSmoke`: a tiny magnetic bcc-Fe-shaped `sp` fixture
  through SCF, accepted LR-01 snapshot, reciprocal eigenpairs, one Gamma q,
  and one frequency.  It is an integration/lifecycle smoke test, not a
  material-accuracy test;
- `TddftProductionDriverFeatureOff`: the same ordinary SCF fixture with no
  `&tddft` group, confirming the feature-off path does not enter the driver.

Useful commands from the configured build are:

```text
ctest --test-dir build --output-on-failure -R 'UnitTddftProductionDriver'
ctest --test-dir build --output-on-failure -R 'TddftProductionDriverSmoke'
```

The direct reproducibility test compares the complete `chiKS`, enhanced, and
loss matrices against direct calls to the same validated services on the same
prepared state.  The material-validation ladder is intentionally not part of
this task.

## TDVAL handoff

The production material adapter now exists and has a passing tiny lifecycle
smoke path.  This unblocks material validation work, but it does not validate
Fe or Ni.  TDVAL-01 must be rerun for bcc Fe and fcc Ni through
`post_processing='tddft'`, with its stated q/mesh/eta convergence, static
invariant, route separation, and literature-comparison evidence.  No Fe/Ni
spectrum, stiffness, damping, or literature agreement is claimed by TDRUN-01.

## Commit

`td-dft: wire clean post-SCF response driver`
