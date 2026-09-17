# TDVAL-K master blueprint — reciprocal collinear transverse TD-DFT

## 0. Authority and scope

Work only on:

- repository: `https://github.com/rslmtoasa/rslmtoasa`
- branch: `fable_v4`
- pinned starting state for this pack: `795b1e488d`

Before each slice, inspect the live branch.  If the branch has advanced, preserve
the scientific contracts below and report any conflict instead of silently
rewriting them.

This campaign validates only the certified initial response tuple:

- collinear magnetic ground state;
- no SOC;
- Pauli/no-SOC transverse response;
- scalar-relativistic LMTO ground state;
- orthogonal `ham_only`;
- second-order/HOH LMTO;
- `sp` / `spd`;
- full site × response-harmonic × radial response space;
- no site-only replacement.

Do not describe this as an exact scalar-relativistic or fully relativistic
response.

## 1. Deferred work

The following are **out of scope for every Luna slice in this pack**:

- native real-space GF TD-DFT production integration;
- TDRUN-02 / `native_rsgf`;
- recursion/Chebyshev response convergence;
- noncollinear response;
- SOC response;
- Sternheimer/Savrasov implementation;
- exact scalar-relativistic nonspherical response;
- fully relativistic `j,kappa` response;
- BdG;
- new XC physics.

Do not modify the native-RSGF modules merely because they exist.

## 2. Authoritative response physics

Do not rederive these equations in implementation tasks.

### 2.1 Bare reciprocal response

For the certified circular sector,

\[
\chi^{KS}_{IJ}(\mathbf q,\omega)
=
\frac{1}{N_k}\sum_{\mathbf k,nm}
\frac{f_{n\mathbf k}-f_{m,\mathbf k+\mathbf q}}
{\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+\mathbf q}+i\eta}
T_{nm;I}T^*_{nm;J},
\]

with the LR-03 circular normalization already implemented by LR-06.

The endpoint is exactly `k+q`, folded by the existing reciprocal-lattice rule.
No nearest-mesh substitution is allowed.

The production bare-response routes are:

1. `backend='lehmann'` — primary transparent baseline;
2. `backend='reciprocal_gf'` — independent reciprocal GF energy-integral
   cross-check.

They are representations of the same bare KS response, not distinct physical
models.

### 2.2 Reciprocal GF numerical parameters

The reciprocal-GF route has two numerical controls distinct from the physical
response broadening:

- `gf_integration_points`: Simpson quadrature resolution;
- `gf_integration_eta`: spectral-discontinuity width.

If `gf_integration_eta <= 0`, the current service chooses `eta/40`.  The
physical response broadening remains `eta`.

Do not tune `eta` to improve backend agreement.

### 2.3 Response-space representation

Use the LR-04 canonical right-weighted representation already implemented.

No task may insert an additional radial weight.

For `spd`, the formally complete response product cutoff is:

\[
L_{\max}^{response}=4.
\]

A smaller `response_lmax` may be used only as an explicitly labelled
approximation/sensitivity test.  It must never be called exact.

The direct production radial mesh is the reference.  Do not invent radial
decimation in this campaign.

### 2.4 Direct ALSDA route

Use the existing KXC-01 kernel unchanged.

The raw static diagnostic is

\[
b_{xc}=K_{xc}g,\qquad
r=\chi^{KS}(0,0)b_{xc}-g.
\]

Record the LR-04 metric residual and rigid-vector overlap.  Do not rescale
`Kxc` from this residual.

### 2.5 Independent Lounis/LCMM route

Use the existing GSR-01 service unchanged.

It is a separate interaction construction.  It must never be implemented by
rescaling direct ALSDA.

If GSR reports rank deficiency, nonspherical leakage, or `BLOCKED`, preserve
that result and stop for physics review.

### 2.6 Dyson/loss

Use the existing TDDY-01 equation unchanged:

\[
(I-\chi^{KS}K)\chi=\chi^{KS}.
\]

The loss convention is

\[
L=-\frac{\chi-\chi^{\ddagger_W}}{2i\pi}.
\]

Do not negate a circular channel or loss quantity to make a spectrum look
positive.

## 3. Claims discipline

Every result must be labelled as one of:

1. algebraic consistency;
2. independent numerical cross-check;
3. numerical convergence;
4. converged-material validation;
5. literature agreement.

Never promote one class into another.

In particular:

- Lehmann/GF agreement does not independently validate the shared radial vertex;
- a passing GSR route does not validate direct ALSDA;
- a later Goldstone-corrected route cannot rescue a broken raw route;
- a visually plausible magnon curve is not a validation by itself.

## 4. Material order

### First material: bcc Fe

Use the existing tracked bcc-Fe ground-state deck as the starting material
definition.  Do not alter lattice constant, XC choice, magnetic setup, or
moments to improve TD-DFT behavior.

The previously accepted Fe foundation uses the current scalar-relativistic
collinear setup and a production-style `spd` basis.  Re-establish the exact live
provenance before response work.

### Second material: fcc Ni

Do not start Ni until Fe has passed the orchestrator checkpoint after the
Dyson/loss stage.

## 5. Validation order

The required order is:

```text
A. implementation preflight
B. accepted Fe ground-state handoff
C. Fe bare chiKS, Gamma
D. Lehmann ↔ reciprocal-GF closure
E. finite-q / folding / covariance
F. k-mesh and eta convergence
G. static ALSDA + independent GSR diagnostics
H. interacting Fe loss
I. Ni replication
J. literature comparison by the orchestrator
```

Do not jump directly to a magnon dispersion.

## 6. Diagnostic q points

Luna may use direct reciprocal-coordinate q vectors for invariants, but must
not assign high-symmetry labels unless an existing code path or input file
already establishes that label.

Recommended diagnostic sequence:

- `q = (0,0,0)`;
- one mesh-commensurate small q;
- its exact `-q`;
- one arbitrary/off-mesh q;
- later, a small-q sequence chosen by the orchestrator.

These are initially **coordinate diagnostics**, not a claimed Gamma-H,
Gamma-N, or other physical symmetry line.

## 7. Broadening discipline

Use the established diagnostic value `eta=0.02 Ry` for the earliest material
smoke checks because that value already exists in the LR-06 evidence envelope.

Later convergence uses a prescribed ladder, for example:

- `0.02 Ry`;
- `0.01 Ry`;
- `0.005 Ry`.

Do not infer an `eta -> 0` extrapolation and do not subtract eta from a width.

## 8. Backend comparison metric

When comparing Lehmann and reciprocal GF on the **same response space and
accepted state**, compute for every `(q,omega)`:

\[
\Delta=\chi^{KS}_{L}-\chi^{KS}_{GF},
\]

\[
d_F=\sqrt{\sum_{IJ}|\Delta_{IJ}|^2},
\]

\[
r_F=
\frac{d_F}
{\max(\|\chi^{KS}_{L}\|_F,\|\chi^{KS}_{GF}\|_F,\epsilon)},
\]

and

\[
d_\infty=\max_{IJ}|\Delta_{IJ}|.
\]

These are representation-level diagnostics because both matrices use the same
LR-04 canonical representation.

Do not impose an arbitrary material PASS threshold in Luna.  Report the values
and their change with GF quadrature controls.  The orchestrator adjudicates.

## 9. Goldstone correction policy

Do **not** enable GCR-01/BES correction anywhere in the initial reciprocal
campaign.

It may only be reconsidered after:

- Fe Lehmann response is numerically converged;
- the raw direct-ALSDA denominator has an identifiable isolated candidate mode;
- rigid overlap and isolation meet the existing GCR-01 contract;
- at least two convergence samples retain that character.

That decision belongs to the orchestrator, not Luna.

## 10. Literature policy

Luna must not tune calculations against literature and must not decide whether
a numerical dispersion is physically correct.

The orchestrator will perform the literature comparison after receiving the
Fe/Ni evidence pack.

## 11. Required source contracts to read before coding

Read the live versions of:

- `docs/LR_TDDFT_CONVENTIONS.md`
- `docs/LR_RESPONSE_SPACE_ALGEBRA.md`
- `docs/LR_PAULI_TRANSITION_VERTEX.md`
- `docs/LR_KS_SUSCEPTIBILITY.md`
- `docs/LR_GF_SUSCEPTIBILITY_CROSSCHECK.md`
- `docs/KXC_ALSDA_TRANSVERSE_KERNEL.md`
- `docs/GOLDSTONE_SUMRULE_INTERACTION.md`
- `docs/GOLDSTONE_EIGENVALUE_CORRECTION.md`
- `docs/TDDFT_DYSON_AND_LOSS.md`
- `docs/TDDFT_PRODUCTION_DRIVER.md`
- `docs/TDDFT_COLLINEAR_REVALIDATION.md`

Historical `BLOCKED` wording caused by native-RSGF production registration is
not a blocker to this **reciprocal-only** campaign.
