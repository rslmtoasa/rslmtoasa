You are working on the current HEAD of the fable_v3 branch of:

https://github.com/rslmtoasa/rslmtoasa

This task belongs to Phase II — Scientific Feature Establishment.

Before editing anything, inspect the CURRENT repository state, including:

- docs/dev/PHASE_I_STABILIZATION.md
- docs/DEVELOPER_MAP.md
- tests/KNOWN_ISSUES.md
- the relevant source modules/submodules
- existing unit/functional/validation tests
- current CTest labels and fixtures

Phase I structural stabilization is complete. Do NOT start another broad
refactoring campaign.

============================================================
GLOBAL PHASE-II RULES
============================================================

1. PHYSICS FIRST

The purpose of Phase II is to establish existing scientific capabilities.

A feature is not Validated merely because:
- it compiles;
- it produces output;
- a smoke test passes;
- a generated reference matches itself.

Validation requires a meaningful mathematical, physical, cross-representation,
or trusted-reference contract.

2. PRESERVE THE RS-LMTO-ASA ARCHITECTURE

RS-LMTO-ASA remains fundamentally a real-space LMTO / Green-function magnetic
code, with reciprocal space as a complementary representation.

Do not create a parallel architecture for the feature being validated.

3. PROTECT THE LEGACY ATOMIC LMTO CORE

Do not edit mature atomic/radial/potential routines in self.f90 and related
files unless the task proves that the root cause lies there and modification
is unavoidable.

Refactor/fix around trusted atomic machinery before changing it.

4. ROOT CAUSES, NOT SYMPTOMS

If validation exposes a discrepancy:
- reproduce it minimally;
- identify the first violated invariant;
- trace backwards to the producing layer;
- fix the owning cause.

Do not:
- relax tolerances until green;
- disable a valid feature combination;
- silently fall back;
- zero/clamp problematic contributions;
- replace NaN/Inf with zero;
- regenerate references without understanding the change.

5. BASELINE FIRST

Before editing:
- run the smallest relevant unit/functional groups;
- record existing observables;
- record current failures separately;
- preserve the same baseline through the implementation.

6. TEST ORACLE POLICY

Prefer, in order:
- analytic identities and limits;
- physical invariants;
- an independent production route already in RS-LMTO-ASA;
- a trusted archived/reference calculation;
- mature external results when genuinely comparable.

Do NOT create a duplicate Python implementation of the target physics as an
oracle.

Python may parse, compare, fit, and report production outputs.

7. VALIDATION TESTS ARE NOT QUICK TESTS

Material/convergence campaigns should normally be labelled `validation` and
kept out of the normal quick developer gate.

Only add a new permanent unit/functional test when it protects a compact,
stable contract exposed by this task.

8. FEATURE STATUS IS SCOPED

Do not promote an entire feature family because one limited case passes.

State exactly what has been established, e.g.:

    "collinear LDA+U with this double-counting convention: Validated"

rather than:

    "LDA+U: Production"

9. NO PERFORMANCE/GPU WORK

Phase II establishes CPU/reference science.

Do not introduce GPU/HIP optimization or performance architecture unless the
task explicitly says otherwise.

10. NO UNRELATED FUTURE-SCOPE DESIGN

Stay within existing scientific capabilities being established here.

11. ONE FOCUSED COMMIT

At completion:
- tick the supplied checklist;
- list files changed;
- list tests/campaigns run;
- report numerical evidence;
- state any remaining limitation honestly;
- make one focused commit with the supplied one-line message.
