# TDCOV-01 — Fix the TD-DFT occupation/Fermi-level lifecycle

## Goal

Make the response occupation state a single, explicit, validated object/contract shared by all transverse-response channels and applicable backends.

This task must **not** change pair-potential moment normalization, circular-response definitions, Ward signs, or q-gauge logic.

## Known defect to reproduce first

In the current `fable_v4` lifecycle:
- reverse eigenpair/Lehmann options are copied before the response Fermi level is resolved;
- their Fermi level can remain the type default `0.0_rp`;
- the primary response path can unconditionally recalculate EF even when `auto_find_fermi=.false.`;
- ground-state EF provenance can be unset/misreported.

Before editing production code, add or run a regression that demonstrates at least one of these failures.

---

## Required design

Define the semantic distinction:

### Ground-state Fermi level
The actual chemical potential inherited from the ground-state/input occupation contract.

### Response Fermi level
Resolved **once** for TD-DFT:
- if `auto_find_fermi=.false.`, honor the supplied/ground-state EF exactly;
- if `auto_find_fermi=.true.`, resolve on the response mesh once and record that this occurred.

All response consumers then receive the same resolved state.

Prefer a narrow helper or small state object over duplicated assignments if that improves clarity, e.g. conceptually:

```fortran
type :: response_occupation_state
   real(rp) :: fermi_level = huge(1.0_rp)
   real(rp) :: electronic_temperature = ...
   logical :: fermi_is_resolved = .false.
   ...
end type
```

Do not introduce a large abstraction if a smaller explicit helper fits the code style better.

### Unset-state rule

`0.0_rp` must not mean "unset".

Use an unmistakable sentinel/validity flag consistent with the codebase. Every backend initializer that requires EF must fail loudly if handed an unresolved value.

---

## Required implementation

1. Trace where the SCF/input EF becomes available.
2. Resolve the response occupation state **before** channel-specific backend options are finalized.
3. Remove unconditional `find_fermi=.true.` semantics from TD-DFT.
4. Propagate the resolved EF and electronic temperature consistently to:
   - primary eigenpair chi0;
   - reverse eigenpair chi0;
   - primary k-space Lehmann/Green options;
   - reverse k-space Lehmann/Green options;
   - real-space options where applicable.
5. Ensure primary and reverse channels cannot diverge in occupation state unless a future API explicitly requests that behavior.
6. Correct provenance:
   - `ground_state_fermi_level_Ry`;
   - `response_fermi_level_Ry`;
   - source/policy describing whether response EF was inherited or recomputed.
7. Keep any internal `huge()` sentinels internal; do not print them as if they were physical metadata.

---

## Mandatory regression tests

### Test A — conspicuous fixed EF
Use a deterministic Hamiltonian-backed response fixture with:
```text
auto_find_fermi = .false.
fermi = -0.123456 Ry
```

Assert:
- reported ground-state EF is the intended ground-state value;
- response EF is exactly `-0.123456` within roundoff;
- primary and reverse backend options use the same EF;
- no canonical response-mesh EF solve overrides it.

### Test B — auto-resolved EF
With `auto_find_fermi=.true.`:
- resolve EF once;
- confirm every consumer sees the same result;
- record source as response-mesh recomputed.

### Test C — unset EF is fatal
Construct backend options with unresolved EF and assert initialization fails with a diagnostic that names the missing occupation state.

### Test D — circular channel equality of thermodynamic metadata
For a `circular_channel='both'` run, assert identical:
- EF;
- temperature;
- occupation tolerance;
- band-window policy;
for primary and reverse channels.

---

## Forbidden shortcuts

- Do not simply add two late reverse-channel EF assignments and leave lifecycle ambiguity intact if the primary user-fixed-EF bug remains.
- Do not retain `0.0_rp` as an unset default.
- Do not change moment signs, kernel signs, q signs, Lehmann denominators, or mode finding in this task.
- Do not weaken electron-count checks to make the new lifecycle pass.

---

## Acceptance checklist

- [x] Reproduced at least one old EF lifecycle failure before patching.
- [x] `auto_find_fermi=.false.` is honored by TD-DFT.
- [x] Auto-resolved EF is computed once.
- [x] Primary and reverse eigenpair backends receive identical resolved occupation state.
- [x] Primary and reverse Lehmann/Green backends receive identical resolved occupation state.
- [x] Real-space path remains consistent.
- [x] Unresolved EF is impossible to consume silently.
- [x] Ground-state EF provenance is real, not zero/default.
- [x] Response EF provenance states inherited vs recomputed policy.
- [x] Fixed-EF regression passes.
- [x] Auto-EF regression passes.
- [x] Unset-EF negative test passes.
- [x] Existing relevant TD-DFT tests pass.

## Required one-line commit message

`td-dft: fix response occupation-state lifecycle`
