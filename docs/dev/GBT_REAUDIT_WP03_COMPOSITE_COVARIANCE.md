# GBT re-audit WP03: composite LMTO covariance

**Date:** 2026-08-27
**Repository:** `fable_v3`
**Validation baseline:** `708b0463` (`test: prove GBT gauge and shifted-k identities`)
**Scope:** fixed-potential, scalar-relativistic bulk `gbt_single_q`

## Verdict

WP03 passes for the composite terms that have an audited implementation:

| Feature | Verdict | Qualified scope |
|---|---|---|
| HOH / second-order LMTO | **supported and exactly tested** | bulk GBT with ordinary `ham_only` metric; fixed-potential RS/k-space operator |
| `Sdot` / two-centre CCOR | **supported and exactly tested** | generated directed `Sdot`; `surface_scalar`, `vmad_scalar`, and `pair_surface` modes |
| generalized overlap | **unsupported and explicitly guarded** | both `generalized_overlap_proxy` and `generalized_overlap_kanpur` |
| local Hubbard-U | **supported and exactly tested** | onsite fixed-potential correction in the collinear rotating frame |

This does not qualify self-consistent Hubbard-U feedback, three-centre or CCD2
CCOR, SOC, intersite Hubbard-V, local-axis GBT, or the complete constrained
SCF workflow. Those remain outside this work package and remain guarded where
the current code has no covariance proof.

## Dependency and test method

The WP02 gate was rerun at the current baseline before this work: the
primitive-link unitarity, bond reversal, closed-loop, reciprocal Hermiticity,
and one-sublattice shifted-k tests pass. The new regression is
`tests/unit/test_gbt_wp03_composite_covariance.py` and is registered as
`UnitGbtWp03CompositeCovariance`.

The test uses an independent dense one-sublattice fixture with complex directed
orbital `S` and `Sdot` blocks, explicit reverse bonds, noncommuting orbital
coefficients, and generic deterministic `(k,q)` points. The ordinary reference
is assembled from independent spin sectors; it never calls a GBT helper and
never phases a completed composite object. The threshold was fixed at
`2e-12`; negative controls must exceed `1e-3`.

The feature sequence is deliberately staged:

1. HOH alone;
2. scalar CCOR alone;
3. pair-surface CCOR alone;
4. HOH+CCOR composition after the individual checks;
5. onsite U locality;
6. generalized-overlap guard/source checks.

The q convention is the repository convention: `q` is Cartesian in units of
`2*pi/alat`, a directed bond uses `alpha=2*pi*q.dot(d)/alat`, and reciprocal
assembly uses `exp(+i*2*pi*k.dot(R))`. Therefore the collinear spin sectors
are compared with ordinary calculations at `k-q/2` and `k+q/2` in the
repository's `(up,down)` ordering.

## HOH derivation and implementation audit

For a directed primitive Hamiltonian bond and local overlap-bar/linearization
objects,

```text
h^G_ij = U_i^H h_ij U_j
o^G_j  = U_j^H o_j U_j
e^G_i  = U_i^H e_i U_i
```

The current GBT builder applies the endpoint link only while lifting raw `S`
and contracting it into `ee`. It then builds onsite `obarm` and `enim` and
forms the directed companion as

```text
eeo_ij = ee_ij obarm_j.
```

The production reciprocal assembler Fourier-sums `ee` and `eeo` separately and
forms `eeo(k) h(k)`. Thus

```text
eeo^G_ij = U_i^H (h o)_ij U_j
h^G o^G h^G = U_i^H (h o h) U_k,
```

with no second GBT phase. `eeo` is a directed intermediate and is not required
to be Hermitian; the contracted second-order Hamiltonian is the Hermitian
object. The stored `eeoee` remains a historical same-bond diagnostic, not the
global HOH operator consumed by the solvers.

The new dense oracle reports the following maxima over q=0 and finite-q
generic points:

| Check | Maximum residual |
|---|---:|
| HOH q=0 / shifted-k matrix identity | `1.931e-16` |
| HOH q=0 / shifted-k eigenvalues | `6.661e-16` |
| HOH reciprocal Hermiticity | `1.653e-16` |
| HOH primitive reverse bond | `3.817e-17` |
| HOH+CCOR composition shifted-k eigenvalues | `6.661e-16` |

The earlier production probes in
`tests/gbt_wp6/cases.json` additionally exercise finite-q bcc-Fe HOH through
the reciprocal and real-space consumers. They pass through the existing WP6
HOH oracle and feature matrix.

## CCOR / Sdot derivation and implementation audit

The GBT builder reads both primitive factors from the same directed neighbor
slot and applies the same endpoint link once:

```text
S_ij    -> S_ij G_ij
Sdot_ij -> Sdot_ij G_ij
D_ij    = W_i (S_ij G_ij) W_j
Ddot_ij = W_i (Sdot_ij G_ij) W_j.
```

The completed scalar or pair-surface CCOR block is then assembled from `D` and
`Ddot`; it is not rephased. For scalar modes the active expression is

```text
Hcc_ij = (VMT-E_lin) [ Ddot_ij + A1_i D_ij + D_ij A1_j ].
```

For `pair_surface`, endpoint `Lambda_i` and `Lambda_j` multiply the same
directed factors on their respective sides. These are onsite rotating-frame
objects, so bond reversal gives the adjoint completed block when the reverse
primitive factors and endpoints are used.

The source audit confirms that `Sdot` uses
`lattice%sdot(j,i,m,ino)`, matching the `sbar(j,i,m,ino)` ordering used for
`S`, and that `build_ccor_pair_block_gbt` calls the lift exactly twice—once for
each primitive factor. It receives the link computed by the main GBT builder;
it does not compute a q phase of its own.

The independent dense results are:

| Check | Scalar CCOR | Pair-surface CCOR |
|---|---:|---:|
| q=0 / shifted-k matrix identity | `9.067e-17` | `1.431e-16` |
| q=0 / shifted-k eigenvalues | `1.110e-16` | `2.220e-16` |
| reciprocal Hermiticity | `1.073e-16` | `8.647e-17` |
| reverse directed bond | `1.646e-17` | `1.512e-17` |

The test includes a negative completed-object phase control. Its residual is
`9.976e-02`, so an accidental second CCOR phase is detectable. The earlier
WP6b production probes cover scalar and pair-surface CCOR, including a
pair-surface run with HOH, and the current WP6b oracle passes.

The following remain explicitly outside the verdict: three-centre CCOR, an
active CCD2 Hamiltonian term, and mixed CCOR-SOC covariance.

## Generalized overlap decision

The current nonorthogonal implementation is not a formal LMTO metric for this
GBT construction. The available proxy assembles `I + eeo(k)`, while the
Kanpur mode is not a completed formal metric path. The reciprocal lifecycle
therefore rejects both generalized-overlap modes for `gbt_single_q` before
eigensolution, with the diagnostic that the available generalized-overlap
representation is not a complete formal GBT metric. The supported GBT
second-order path is `reciprocal_mode='ham_only'`, whose metric is the identity.

No guard was weakened in WP03. The new test checks both the source guard and
the production reciprocal entry point, and the existing negative probes remain
the integration-level rejection coverage.

## Local Hubbard-U decision

`calculate_hubbard_u_potential_general` constructs the current onsite
spin-diagonal correction from `potential%ldm`. In the GBT bulk builder it is
added to `ee(:,:,1,ntype)` at the onsite `m=1` slot after primitive contraction.
That slot has `d=0` and `G_ii=I`, so the correction has no translational q phase
and is covariant as a local rotating-frame operator. The dense test verifies
that adding this correction preserves the shifted-k identity, with maximum
residual `8.653e-17`.

This is fixed-potential operator evidence only. The separate SCF feedback path
that updates `ldm` and regenerates Hubbard-U has not been qualified, so no
self-consistent-U production claim is made.

## Guards retained

The WP03 audit confirms that unsupported combinations continue to stop early:

- SOC in `gbt_single_q`;
- unaudited `local_axis` in `gbt_single_q`;
- intersite Hubbard-V;
- generalized-overlap proxy and Kanpur modes;
- CCOR without generated directed `Sdot` and the required screening inputs.

No double-gauging or empirical energy correction was introduced.

## Reproduction and evidence boundary

Commands run against the current tree:

```text
python3 tests/unit/test_gbt_wp03_composite_covariance.py
python3 tests/unit/test_gbt_wp6_hoh.py
python3 tests/unit/test_gbt_wp6_ccor.py
python3 tests/unit/test_gbt_wp6_matrix.py
python3 tests/unit/test_gbt_oracles.py
ctest --test-dir build-ci-workflow --output-on-failure \
  -R 'UnitGbt(Wp6|Oracles|Structure|Wp01)'
ctest --test-dir /tmp/rslmto_fable_v3_wp03 --output-on-failure \
  -R 'UnitGbt(Wp02GaugeShiftedK|Wp03CompositeCovariance)'
```

All commands pass. The clean WP03 build also compiled the full `rslmto`
library and the WP02 Fortran executable. The existing build configuration used for the CTest run is
GNU Fortran Release with MPI and OpenMP enabled and the lean unit-test targets
already configured. The new WP03 test is Python/NumPy-only and does not require
an SCF fixture or a new compiler-dependent executable.

Not established by WP03:

- commensurate-supercell folding;
- rotating-frame/lab-frame density and moment semantics;
- physical constraint-field energy bookkeeping;
- corrected constrained frozen magnons;
- cone/k-grid convergence, small-q curvature, LKAG closure;
- unrestricted or constrained GBT SCF;
- SOC, intersite V, local-axis, generalized-overlap, or broader CCOR/U modes.

## Completion checklist

- [x] HOH covariance derived from the actual `ee -> eeo -> eeo(k)h(k)` path.
- [x] HOH q=0, shifted-k, reverse-bond, Hermiticity, and eigenvalue checks.
- [x] `Sdot` directed support/order verified.
- [x] Scalar and pair-surface CCOR q=0, shifted-k, reverse-bond, Hermiticity, and eigenvalue checks.
- [x] Generalized-overlap covariance audited; unsupported modes remain guarded.
- [x] Local Hubbard-U fixed-potential status decided from source and dense evidence.
- [x] No unsupported feature silently enabled.
- [x] No double-gauging introduced.
- [x] Per-feature evidence table written.
- [x] Minimal regression test added and registered in CTest.
