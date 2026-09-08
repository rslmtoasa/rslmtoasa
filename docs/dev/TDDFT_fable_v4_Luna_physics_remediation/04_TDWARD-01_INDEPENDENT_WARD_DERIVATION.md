# TDWARD-01 — Derive an independent LMTO Ward identity and B_xc provenance

## Prerequisites

TDCOV-01 through TDCOV-03 must be green.

## Nature of this task

This is deliberately an **investigative/derivation task**.

Do not redesign the Goldstone writer yet. Do not choose a Ward sign because it makes Fe pass. Do not reconstruct an allegedly independent B_xc from the same `K_xc m` relation being tested.

The output of this task is a physically explicit derivation and a code-level provenance map that TDWARD-02 can implement safely.

---

## Goal

Determine the exact Ward identity realized by the RS-LMTO-ASA conventions, including all signs and factors of two.

The documentation must reconcile:

1. the Hamiltonian spin decomposition used by the code;
2. the definition/sign of the ground-state exchange-correlation spin splitting;
3. the code's `O+ = sigma_x + i sigma_y` and `O- = sigma_x - i sigma_y` convention;
4. any source factor such as 1/2;
5. the definition of `chi_KS`;
6. the definition of pair-potential `Xi`;
7. the relation between an independent ground-state `B_xc` and the Goldstone magnetization vector.

---

## Required source trace

Trace at least:
- the ground-state magnetic LMTO potential parameters;
- Hamiltonian spin decomposition;
- XC response provider;
- legacy transverse kernel construction;
- pair-potential tangent construction;
- circular response conventions;
- static divided-difference chi0;
- existing Ward/Goldstone diagnostic code.

Identify where a ground-state exchange splitting / B_xc can be obtained **without** defining it as `K_xc m` from the same response kernel.

If no independent quantity currently exists in a sufficiently explicit form, state that clearly and propose the narrowest plumbing needed.

---

## Required derivation document

Create:

`docs/TDDFT_WARD_INDEPENDENT_DERIVATION.md`

It must contain:

### 1. Ground-state Hamiltonian convention
Write the local spin-dependent Hamiltonian/potential in the exact convention used by the code, e.g. schematically

\[
H = H_0 I + \mathbf H_1\cdot\boldsymbol{\sigma},
\]

but use the actual signs and factors from the implementation.

### 2. Definition of B_xc
State whether the code-level B_xc corresponds to:
- `V_up - V_down`;
- `(V_up - V_down)/2`;
- the negative of either;
- or another LMTO-specific transformed quantity.

Do not infer this from filenames/comments alone; derive from the Hamiltonian.

### 3. Circular convention
Derive the action of
\[
O^\pm=\sigma_x\pm i\sigma_y
\]
and source normalization.

### 4. Static Ward identity
Derive the exact implemented identity:
\[
\chi_{\rm KS}(0,0) B_{\rm xc} = \pm m
\]
or its matrix/site-resolved analogue.

Fix the sign from definitions.

### 5. Pair Xi identity
Derive how
\[
\Xi m = m
\]
follows from the pair-potential representation and how it transforms under global spin reversal.

### 6. Independent provenance
Name exact variables/routines from which TDWARD-02 should obtain:
- `m`;
- independent `B_xc`;
- `chi_KS`;
- `Xi`.

### 7. Falsification tests
Specify at least:
- +z one-site;
- -z one-site;
- +/-z two-sublattice;
- static direct vs pair-Xi comparison where mathematically justified.

---

## Optional code allowed

You may add a narrow non-production accessor/helper if absolutely necessary to expose the independent ground-state B_xc for tests.

Do not yet replace production Goldstone diagnostics.

---

## Stop condition

If the exact B_xc convention or sign cannot be derived unambiguously from the current ground-state Hamiltonian/potential code, **stop and report the ambiguity**. Do not proceed by fitting a sign/factor to Fe.

That is a successful outcome for this task if documented precisely.

---

## Acceptance checklist

- [x] Ground-state spin Hamiltonian convention derived from code.
- [x] B_xc sign and factor convention derived.
- [x] Circular operator/source normalization reconciled.
- [x] Exact static Ward identity written with justified sign.
- [x] Pair-Xi Goldstone identity derived.
- [x] Independent B_xc provenance identified, or ambiguity documented.
- [x] +z/-z transformation explicitly shown.
- [x] Multisublattice transformation explicitly shown.
- [x] Falsification tests specified.
- [x] No production sign/factor changed merely to make a test pass.

## Execution record

The source trace found an unambiguous half-difference XC convention and an
independent direct XC provenance in `VXC0SP`/`bxc_spin_moment`.  It also found
that the current Ward diagnostic receives `K_perp*m` as its source vector, so
that path is documented as derived rather than independent.  The narrowest
next step is explicit source plumbing for TDWARD-02; no production sign or
factor was changed here.

Focused checks passed:

```text
UnitTddftWardConventions, UnitTddftWard, UnitTddftGoldstone,
UnitTddftDirectXi, UnitLmtoPairPotential, UnitLmtoMagneticTangents
```

`UnitLmtoPairPotential` reported maximum error `8.8558E-10`, centered tangent
error `5.0875E-10`, and one-site `+z/-z` pair Goldstone eigenvalues `1/1`.
`UnitLmtoMagneticTangents` reported maximum error `4.8566E-10`.

## Required one-line commit message

`docs: derive independent TDDFT Ward quantities`
