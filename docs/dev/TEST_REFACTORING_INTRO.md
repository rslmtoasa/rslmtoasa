You are working on the current fable_v3 branch of:

https://github.com/rslmtoasa/rslmtoasa

This task is part of the RS-LMTO-ASA test-harness stabilization program.

GENERAL RULES

1. Work from the current repository state. Inspect the relevant source, tests,
   CMake registration, runners, manifests, fixtures, and existing references
   before editing anything.

2. This task is primarily about tests/test infrastructure. Do not modify
   production physics code unless the task explicitly requires it or a genuine
   production defect makes the requested test impossible.

3. Treat root causes, not symptoms.
   - Do not disable valid functionality to make tests pass.
   - Do not add guards around broken valid combinations and call that a fix.
   - Do not relax tolerances merely to obtain green tests.
   - Do not silently fall back to another implementation.
   - Do not update reference values merely because the new run differs.
   If a production defect is exposed, identify it clearly and preserve the
   evidence rather than hiding it.

4. Protect the legacy atomic LMTO routines, especially mature atomic/radial/
   potential machinery in self.f90 and related files. Do not edit those routines
   for test convenience, cleanup, or modernization.

5. Testing policy:
   - Prefer small Fortran unit tests and small functional tests of real
     production routines.
   - Prefer analytic invariants, independently implemented production routes,
     and trustworthy recorded baselines.
   - Do NOT create a Python reimplementation of a physics algorithm merely to
     serve as an oracle.
   - Python is fine for test orchestration, parsing, comparison, and harness
     mechanics.
   - Keep tests lean and modular.

6. Baseline-first workflow:
   a. Identify the existing tests affected by this task.
   b. Run them before editing.
   c. Record their pass/fail status and relevant numerical outputs.
   d. Implement the change.
   e. Run the same tests again.
   f. Run the smallest relevant broader test group.
   g. Explain every changed numerical reference before accepting it.

7. Do not perform unrelated cleanup.

8. Preserve the current modular/class-like Fortran style and existing
   parent-module/submodule architecture.

9. When the task is complete:
   - tick every completed checkbox in the task checklist;
   - leave unchecked anything genuinely incomplete;
   - summarize files changed, tests run, and numerical changes;
   - make one focused commit using the supplied one-line commit message.

