# WP6b: CCOR covariance audit

Date: 3 August 2026

## Verdict

The WP6b CCOR slice passes.  GBT CCOR is enabled for the implemented
two-centre modes (`surface_scalar`, `vmad_scalar`, and `pair_surface`) because
each active term now has an independent dense oracle, directed reverse-bond
test, finite-k Hermiticity test, and production smoke evidence.  The old
completed-object rotations have been deleted.

This is a partial G6 PASS for WP6b only.  Final G6 remains open until WP6c and
the WP6 integration feature matrix are complete.

## Primitive-factor derivation

For directed physical bond `i -> j`, use the same endpoint link as the main
Hamiltonian,

```text
G_ij = (U_i^0)^dagger Rz(q.d_ij) U_j^0,
G_ji(-d) = G_ij^dagger.
```

The gauge is owned by both primitive structure-like factors:

```text
S_ij       -> S_ij tensor G_ij
Sdot_ij    -> Sdot_ij tensor G_ij
D_ij       = W_i (S_ij tensor G_ij) W_j
Ddot_ij    = W_i (Sdot_ij tensor G_ij) W_j
```

`W_i` and `W_j` are unchanged diagonal potential factors in the collinear
rotating frame.  `D`, `Ddot`, and the completed CCOR block are composites and
are never phased afterward.

Writing `A0` and `A1` for the orbital-diagonal `ccd0` and `ccd1` coefficients,
the scalar modes are

```text
Kcc_ij = Ddot_ij + A1_i D_ij + D_ij A1_j
Kcc_ii += W_i A0_i W_i
Hcc_ij = lambda Kcc_ij,  lambda = VMT - E_lin.
```

For `pair_surface`, the endpoint spin matrices `Lambda_i` and `Lambda_j` are
onsite diagonal matrices in the same rotating frame:

```text
Hcc_ij = 0.5 (Lambda_i Ddot_ij + Ddot_ij Lambda_j)
        + Lambda_i A1_i D_ij + D_ij A1_j Lambda_j
Hcc_ii += Lambda_i W_i A0_i W_i.
```

Under bond reversal, primitive Hermiticity and the common link give
`D_ji(-d)=D_ij^dagger` and the same identity for `Ddot`.  Swapping all endpoint
onsite factors therefore gives `Hcc_ji(-d)=Hcc_ij^dagger`, which makes the
Bloch sum Hermitian.

For the onsite slot, `d=0`, `alpha=0`, and the shared endpoint reference frame
gives `G_ii=I`.  The onsite `A0` term is evaluated directly in the signed
collinear spin channels.  It receives no translational phase.

## Term map

| Object or coefficient | Class | Frame/gauge ownership |
|---|---|---|
| raw `sbar` / `S_ij` | directed bond | common `G_ij` before contraction |
| normalized `sdot` / `Sdot_ij` | directed bond | same `G_ij`, independently, before contraction |
| `wx0`, `wx1`, `W_i` | onsite | unchanged signed collinear rotating frame |
| screening alpha, alpha-dot, WS radius | onsite/type data | inputs to `ccd`; no GBT phase |
| `ccd0` (`A0`) | onsite | onsite correction only |
| `ccd1` (`A1`) | onsite | left/right endpoint multiplier of a directed composite |
| `ccd2` | onsite diagnostic | computed and logged, not an active Hamiltonian term |
| scalar `VMT-E_lin` | global onsite scalar | multiplies completed scalar CCOR expression |
| pair `Lambda_i`, `Lambda_j` | onsite endpoint matrices | signed collinear rotating frame |
| `D_ij`, `Ddot_ij` | directed composite | derived from already linked primitive factors |
| scalar and pair-surface CCOR summands | directed composite | endpoint-ordered contractions shown above |
| onsite `W_i A0_i W_i` term | onsite composite | alpha=0, common reference axis |
| `eecc`, `hallcc` | directed additive operator | completed CCOR; never rephased |
| `ee+eecc`, `hall+hallcc` work arrays | composite implementation cache | no independent gauge operation |
| `Hcc(k)` | reciprocal composite | ordinary Bloch sum of `eecc` |

Three-centre CCOR, an active `ccd2` Hamiltonian contribution, and mixed
CCOR-SOC terms remain outside the implemented CCOR model.

## Ordinary NC and explicit texture

The ordinary path is q-agnostic.  `periodic_nc` reads per-type endpoint
moments.  `explicit_texture` calls the shared site-indexed endpoint-moment
provider used by the main Hamiltonian.  GBT is rejected if it reaches this
contracted Pauli-component path; it uses the separate primitive-factor builder.

## Deletion list

Deleted from `hamiltonian_ccor.f90`:

1. the `kcomp(:,:,1:2)` full-angle transverse rotation after scalar CCOR
   contraction;
2. the `hcomp(:,:,1:2)` full-angle transverse rotation after pair-surface
   contraction;
3. their `alpha_ss`, bond-vector, and transverse temporary variables;
4. the representation-blind `ccor_apply_spin_spiral` helper and parent
   interface.

No CCOR code now computes `q_ss . bond`, and no contracted CCOR object is
phased.

## Enabled and rejected combinations

Enabled after this audit:

- GBT plus scalar two-centre CCOR (`surface_scalar` or `vmad_scalar`);
- GBT plus `pair_surface` two-centre CCOR;
- the same additive CCOR operator with first-order or audited HOH construction;
- reciprocal `ham_only` and the bulk real-space recursion consumer.

Early rejection remains for:

- a legacy structure backend in GBT;
- missing `Sdot` / `strux_want_sdot`, invalid screening, or invalid VMT mode;
- nonzero-q GBT with SOC (therefore also mixed CCOR-SOC);
- incomplete generalized-overlap modes;
- the still-unimplemented local-cluster GBT assembler and other WP6c terms.

## Evidence

Independent oracle (`tests/unit/test_gbt_wp6_ccor.py`):

| Check | Maximum error |
|---|---:|
| scalar dense expression | `1.153e-16` |
| pair-surface dense expression | `1.417e-16` |
| reverse directed bond | `3.023e-16` |
| onsite alpha=0/common frame | `5.565e-17` |
| scalar k-space Hermiticity | `1.231e-16` |
| pair-surface k-space Hermiticity | `9.743e-17` |

Production probes (`tests/gbt_wp6b/cases.json`):

- `wp6b_ccor_scalar_k4`: PASS; nonzero `maxabs(Hcc)=3.3790e-2`;
  pre-eigensolver `max|H-H^H| <= 2.2276e-16`.
- `wp6b_ccor_pair_surface_k4` with HOH: PASS; nonzero
  `maxabs(Hcc)=1.4334e-1`; pre-eigensolver
  `max|H-H^H| <= 4.4431e-16`.
- `wp6b_ccor_scalar_rs12`: PASS.
- `wp6b_ccor_explicit_texture_q0`: PASS; its ordinary-path q=0
  `maxabs(Hcc)=3.3790e-2` matches the GBT q=0 value.
- strict missing-Sdot negative probe: expected early fatal with
  `gbt_single_q with CCOR requires generated Sdot and strux_want_sdot=.true.`.

Ordinary regression:

- `bccFe_chebyshev_fast_ccor_2c`: PASS against its existing reference.

Complete unit label run: 19/19 PASS.

## Gate recommendation

- WP6b CCOR slice: **PASS**.
- Partial G6 through WP6a+WP6b: **PASS for the audited slices**.
- Final G6: **OPEN / not yet passable**, pending WP6c and integration.
