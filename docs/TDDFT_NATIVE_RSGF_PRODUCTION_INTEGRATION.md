# TDRUN-02 — Native RSGF production-driver integration

Audit date: 2026-09-11
Branch: `fable_v4`
Starting HEAD: `795b1e4` (`docs: close native RSGF capability gates`; working tree dirty from TDRUN-02)

## Verdict

**PASS for TDRUN-02 at the registered production-integration scope.**

The normal post-SCF TD-DFT driver accepts the explicit `native_rsgf`
selector, hands the accepted SCF/LR-01 state to the registered native provider,
returns the public LR-04 response object, and sends that object through the
existing KXC/Dyson path. R4 remains pending: no Fe/Ni material-validation claim
is made by this task.

## Authoritative preflight

`docs/RSGF_CAPABILITY_CLOSURE.md` is the gate for this run. Its entry state was:

```text
R0 = PASS   native coefficient-GF representation
R1 = PASS   Pauli radial endpoint augmentation
R2 = PASS   native finite/provider bare-response backend
R3 = PENDING production-driver registration
R4 = PENDING Fe/Ni material validation
```

The historical `BLOCKED` statements in LR-REP-00 and the earlier TDVAL
preflight were not treated as campaign vetoes because the closure ledger
explicitly records their R0–R2 owners as closed/pass at the certified scope.

## Integration scope and call graph

TDRUN-02 adds orchestration and registration only. It does not add a response
equation, radial augmentation formula, representation conversion, recursion or
Chebyshev algorithm, XC/kernel/Dyson implementation, or spectrum/mode reducer.

```text
pre_processing_bravais
  -> accepted self-consistent state
  -> run_tddft_production
       -> validate production capability
       -> accepted reciprocal eigenpair/EF/T handoff
       -> tddft_native_rsgf_provider%initialize
       -> evaluate_tddft_production_sweep
            -> evaluate_lr_rs_gf_susceptibility
                 -> provider%get_pair
                 -> existing LR-04 endpoint augmentation
                 -> existing native RSGF bubble and q phase
                 -> LR-04 canonical chiKS object
            -> existing KXC service
            -> existing Dyson service
       -> provenance-bearing output
```

The native route is selected independently from the existing routes:

| selector | route |
| --- | --- |
| `lehmann` or `spectral` | LR-06 spectral service |
| `reciprocal_gf` | reciprocal-GF service |
| `native_rsgf` | registered coefficient-GF provider and native RSGF service |

No selector aliases the native route to reciprocal GF, and feature-off input
continues to leave TD-DFT disabled.

## Production capability gate

The registered baseline fails closed unless all of the following hold:

- collinear, no SOC, orthogonal `ham_only` state;
- second-order/HOH Hamiltonian and supported `sp`/`spd` radial basis;
- no additive/local-axis/unsupported Hubbard or generalized-overlap operators;
- complete `lattice%ijpair` pair workset;
- serial replicated baseline workset (`numprocs = 1`);
- `control%recur = block` or `chebyshev`, matching the provider selector;
- accepted LR-01 radial/Pauli snapshots and a valid common EF/temperature state.

The production input is explicit:

```text
backend = 'native_rsgf'
native_rsgf_provider = 'block' | 'chebyshev' | 'auto'
```

`auto` follows the existing recursion selector. The block provider reuses the
accepted on-site block coefficients for an all-on-site pair workset and uses the
existing pair-recursion service for off-site worksets. It supplies directed
coefficient-space `G_ij`, `G_ji`, and `hgamma` blocks to the public RSGF
provider contract; it does not reduce the response to site-only data.

## Accepted-state handoff and response identity

The driver consumes the converged structural data, accepted potential,
LR-01 radial/Pauli snapshot, XC provenance, Fermi state, temperature and
occupations, basis/magnetic metadata, and existing response radial/angular
metadata. It does not run SCF, search for an independent EF, or reconstruct
radial XC data.

The native service returns the same `lr_ks_susceptibility_result` response shape
used by the reciprocal service. The result is the complete radial/angular
LR-04 canonical matrix, including the existing q/translation convention and
EF/T/integration provenance. Output records:

```text
backend = native_rsgf
native_rsgf_provider = block|chebyshev|auto
native_rsgf_provenance = provider identity, pair workset, q phase, EF/T, eta and quadrature
native_rsgf_route = coefficient-GF -> endpoint augmentation -> LR-04 -> common KXC/Dyson
```

## Evidence

Build and focused checks passed after the integration changes:

```text
cmake --build build -j2                         PASS
UnitTddftProductionDriver                     PASS
UnitLrRsGfRepresentation                      PASS
UnitLrGfEndpointAugmentation                  PASS
UnitLrRsGfSusceptibility                      PASS
TddftProductionDriverFeatureOff               PASS
UnitTddftProductionDriverRejectSoc             PASS (expected rejection)
UnitTddftProductionDriverRejectGeneralizedOverlap PASS (expected rejection)
TddftProductionDriverSmoke                    PASS (reciprocal regression)
TddftProductionDriverNativeSmoke              PASS
```

The native production smoke uses the accepted Fe fixture, `control%recur=block`,
one complete on-site pair, `response_lmax=0` to keep the full 495-point radial
mesh at smoke scale, one q point, one frequency, and three Simpson energy
points. It completed in 38.06 seconds and emitted the native provider identity
and common-route provenance.

The unit reproducibility check calls the public native RSGF service directly
with the finite dense coefficient provider and then calls
`evaluate_tddft_production_sweep` with the same prepared request. The driver’s
native `chiKS` agrees with the direct service result to the asserted `2e-12`
tolerance. This is a service/driver equivalence check, not material validation.

The reciprocal production smoke also passed in the same run, and the feature-off
test passed, establishing selector isolation and reciprocal-route preservation
at the existing regression scope.

## Claim level and next task

This closes **R3 — production-driver registration**. It does not close R4 and
does not promote Fe/Ni physics, provider convergence, SOC, noncollinear,
generalized-overlap, or material agreement claims.

```text
R0 = PASS
R1 = PASS
R2 = PASS
R3 = PASS
R4 = PENDING
```

Next action: rerun **TDVAL-01R** for Fe/Ni material validation. Until that run
passes, the authoritative ledger must retain `R4 = PENDING`.
