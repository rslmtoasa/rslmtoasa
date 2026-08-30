You are working on the current HEAD of the fable_v3 branch of:

https://github.com/rslmtoasa/rslmtoasa

This task belongs to Phase I — Structural Stabilization of RS-LMTO-ASA.

Before editing anything, inspect the CURRENT repository state. Do not assume
that file structure, tests, line numbers, or implementation details are
identical to earlier plans.

============================================================
GLOBAL DEVELOPMENT RULES
============================================================

1. PRESERVE THE SCIENTIFIC ARCHITECTURE

RS-LMTO-ASA is primarily a real-space LMTO / Green-function magnetic
electronic-structure code.

The newer reciprocal-space functionality is an organic complementary
representation of the same underlying LMTO physics.

Do not reorganize the code around reciprocal space, accelerators, TD-DFT,
or generic software abstractions.

2. PREFER THE EXISTING HOUSE STYLE

Prefer:
- physical derived types;
- constructors/destructors;
- parent modules with implementation submodules;
- concrete reusable helpers;
- reusable workspaces;
- small dependency-light mathematical kernels.

Do NOT introduce generic factories, provider hierarchies, dependency-injection
systems, or new solver abstraction layers unless this specific task has a
demonstrated requirement for one.

3. PROTECT LEGACY ATOMIC LMTO PHYSICS

Mature atomic/radial/potential LMTO routines in self.f90 and related files are
protected scientific infrastructure.

Do not edit them for:
- stylistic cleanup;
- modernization;
- formatting;
- abstraction;
- test convenience.

If this task genuinely requires modifying such a routine, stop first and
establish:
- why the change is unavoidable;
- its numerical contract;
- a pre-change baseline;
- the minimum possible patch.

Prefer refactoring AROUND trusted legacy physics rather than THROUGH it.

4. ROOT CAUSES, NOT SYMPTOMS

Do not make a failing valid pathway green by:
- disabling it;
- returning early;
- silently falling back to another algorithm;
- zeroing a problematic contribution;
- clipping an unphysical value;
- suppressing NaN/Inf;
- relaxing tolerances without justification;
- removing assertions;
- updating reference numbers simply to match new output.

If a feature combination is physically valid and intended to work, a failure
is a defect.

Find the first violated invariant, trace it to the producing layer, and fix
the owning cause.

A guard is acceptable only for a genuine physical/mathematical incompatibility,
a deliberately unavailable implementation, or corrupted internal state.

5. BASELINE FIRST

Before modifying production code:

a. identify the smallest relevant CTest labels/tests;
b. run them;
c. record pass/fail state;
d. record relevant numerical observables;
e. identify known failures separately.

After the change:

a. run exactly the same focused tests;
b. compare with the recorded baseline;
c. run the broader lean relevant group;
d. run the full lean unit suite where practical.

Do not regenerate references until every numerical change is understood.

6. TEST QUALITY

Use the existing production implementations.

Do not create a Python reimplementation of the target physics as an oracle.

Good oracles are:
- analytic invariants;
- mathematical limits;
- another independent production implementation;
- an already-trusted recorded baseline;
- existing canonical unit/functional tests.

Python may be used for orchestration, parsing, and comparison only.

Keep new tests small and focused.

7. STRUX POLICY

Normal test fixtures now use strux_lib unless a test explicitly exists to
exercise or compare the legacy structure-constant backend.

Do not change the application's production default as part of STAB work.

8. FEATURE COMBINATIONS

Do not reduce the supported feature space merely because a combination exposes
a bug.

Classify any incompatibility as:
- physically/mathematically invalid;
- meaningful but not implemented;
- intended to work but broken.

Only the first two justify a permanent rejection.

9. SCOPE DISCIPLINE

Perform only the task described below.

Do not opportunistically:
- rename unrelated routines;
- reformat large files;
- modernize unrelated legacy code;
- change numerical algorithms during a pure structural task;
- combine optimization and refactoring unless explicitly requested.

10. COMPLETION

When done:
- tick completed checklist items;
- leave genuinely incomplete items unticked;
- summarize the root cause or structural change;
- list files changed;
- list tests run and results;
- state whether any numerical output changed and why;
- make ONE focused commit using the supplied commit message.

If the requested change cannot be completed safely, do not invent a workaround.
Document precisely what blocked it.

