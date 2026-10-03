# RS-LMTO-ASA — orientation for a fresh AI session

Start here, then read:

- **[`docs/DEVELOPER_MAP.md`](docs/DEVELOPER_MAP.md)** — entry points, class
  chains, kernel inventory, testing map. The fastest way to find where a
  workflow lives.
- **[`docs/ROADMAP.md`](docs/ROADMAP.md)** — planned feature work and current
  status.
- **[`docs/DECISIONS.md`](docs/DECISIONS.md)** — closed campaign conclusions
  and archive locations.
- **[`tests/README.md`](tests/README.md)** — suite index, CI trigger/matrix
  strategy, GPU coverage tiers.
- **[`tests/KNOWN_ISSUES.md`](tests/KNOWN_ISSUES.md)** — bugs found via
  coverage work, deliberately left unfixed pending a dedicated task.

When a task names a spec file, that spec governs the task's technical
content; this file governs how you work.

**Removed campaign code.** The linear-response campaign's source, tests
and docs were removed from `fable_v4b`; they remain at tag
`lr-campaign-archive-2026-10`. Do not restore, copy or consult them for
new work: their numbers, conventions and tests are not references.

## The rules that matter most

1. **Bit-level behavior is the contract.** `tests/regression` and the
   `tests/{scf,postproc}` example suites must pass at the same tolerances
   before and after every change. Run them before starting and after every
   task.
2. **One task, one commit.** Don't batch unrelated changes into one commit. Only allow single line commit messages and do not add "Co-authored by".  If too complicated, ask.
3. **KISS.** Write the smallest code that does the task. See "Simplicity".
4. **Class-based architecture stays.** Derived types, constructors,
   `restore_to_default`, type-bound procedures. New code follows the same
   pattern.
5. **`self.f90` and `symbolic_atom.f90` are off-limits.** Physics-dense
   legacy code scheduled for a separate audit — edit around the edges only.
6. **Stay inside the defined task.** See "Scope discipline" — this one is
   about token budget, and it is not optional. Do not expand the task.
7. **Sprint, no marathons.** If a check still fails after two focused
   attempts, stop and report what you measured. Never get around a failure by
   adding options, routes or special cases, or by loosening a check.
8. **You do not grade your own work.** See "Verification integrity".

## Verification integrity — no self-oracling

A check is evidence only if it could have failed and its expected value came
from somewhere other than the code under test. The archived campaign failed
on exactly this point: its checks confirmed their own arithmetic in cases
where the answer was fixed by construction.

1. **Expected values come from outside the code under test:** the
   developer-owned oracle directories (currently `tests/spin_response/oracles/`), committed regression
   references, a closed-form expression written independently in the test, or
   a cited publication. Never from running the code under test, calling the
   routine being tested, or re-deriving its formula inside the test.
2. **Never move the yardstick.** Do not edit, regenerate or re-baseline
   reference data, expected values or tolerances to make a check pass, and do
   not delete, skip or weaken a failing test. Developer-owned oracle
   directories are read-only. If a reference looks wrong, stop and report it.
3. **Do not write the oracle for your own work.** If a task needs a reference
   that does not exist, say so and stop. Unit tests of your own internal
   plumbing are fine but never count as validation. When the developer asks
   you to draft an oracle for review, that is the whole task: do not also
   implement the code it checks.
4. **Every new test must be able to fail.** For each test, name in one line
   the bug it would catch. Where practical, show it failing once with that bug
   injected locally, not committed.
5. **Label tautologies.** A quantity defined by the identity being checked, a
   1x1 eigenvector, a residual that aliases another, or any other check that
   holds by construction is reported as `by construction` and never counted as
   a pass.
6. **Trivial cases are not coverage.** One-site, q = 0, omega = 0, U = 0 or
   1x1 cases may be included but never as the only evidence; include at least
   one case whose answer is not fixed by symmetry or definition.
7. **One code path.** No special-case branch or solver for test inputs,
   q = 0, omega = 0 or one-site cells that production does not take.
8. **Report numbers, not status.** Never write "certified", "validated",
   "closed", "exact", "passes the gate" or similar in code, comments, docs or
   commit messages. Report measured values; the developer assigns status.

## Simplicity — no overcomplicated code

1. **Nothing that was not asked for.** No new namelist keys, options, modes,
   backends, selectors, providers, registries, capability gates, contract
   types or abstraction layers unless the task requests them by name.
2. **Rule of two.** Generalize only when a second real use exists in the
   same task. One implementation of each thing.
3. **Reuse, don't duplicate.** Use the existing routine (occupations, Fermi
   solve, constants, eigenpairs). If it is private or awkward, propose
   exposing it; never copy it into a new module.
4. **Plain over clever.** A subroutine before a new type; a new type only for
   state with a lifetime. Numerical kernels take plain arrays.
5. **No dead or scaffolding code.** No `if (.false.)` blocks, commented-out
   alternatives, "kept for reference" copies or unused arguments.
6. **Test code lives in `tests/`.** No fixtures, mock providers or test hooks
   exported from production modules.
7. **Plain output.** Columns with a one-line header naming units. No
   provenance essays in data files. Never print an uninitialized or sentinel
   value; print `n/a`.
8. **Plan before large changes.** If a change will add more than about 200
   lines or touch more than three files, first reply with a plan (files,
   procedures, estimated lines) and wait. If the task states a size budget,
   stop before exceeding it.

## Reporting — how every task ends

- A table: check | measured | expected or tolerance | command.
- A separate list of what holds by construction, and what was not checked.
- `git diff --stat` and `wc -l` of new or changed files.
- At most five lines of notes. No narrative summary.

## Documentation and research-ticket policy

- Agent prompts and prompt packs are never committed.
- A task's record is the report above. Do not write campaign reports,
  certification records or status documents.
- A research or validation campaign ends with one line in
  `docs/DECISIONS.md`; its reports are deleted at campaign close (they stay
  at the campaign's tag).
- Every research ticket states at creation whether its result will be
  promoted into a named module/submodule or archived at a tag.
- Campaigns add scripts and data under `tests/validation/` only — never new
  source modules, namelist values or `backend` strings.

## Scope discipline — do only the task that was asked for

The user defines the task. Work on that, finish it, report, and **stop**. Do not
extend the scope because an adjacent problem looks tractable or interesting.
This is a budget constraint, not a limit on capability: exploratory work burns
tokens that were not authorized for it.

**Green-light rule.** Anything beyond the stated task needs the user's explicit
approval *before* you start it. That includes: fixing a bug you found on the
side, adding test cases or examples that were not requested, refactoring code
you happened to read, widening a validation study, and "while I'm here"
cleanups. Propose it in one or two sentences and wait.

**When you find something real mid-task** — and you will, this codebase has
live bugs — the correct move is:

1. Record it in `tests/KNOWN_ISSUES.md` with what you actually verified,
   clearly separated from what you are guessing.
2. Leave a visible comment at the code site if it is load-bearing.
3. Mention it in your report as a *proposal*, with a size estimate.
4. Carry on with the original task.

Do **not** silently fix it, and do **not** start triaging it. If it blocks the
task you were given, say so and stop rather than working around it.

**Honesty about untriaged findings.** If you did not diagnose something, say
so plainly and name your own test setup as a suspect where that is a real
possibility. A hand-built input deck is a live suspect until a known-good deck
reproduces the symptom. Never present an unverified guess with the same
confidence as a measured result.

**Corollary — don't invent studies either.** Running the regression suite
before and after is rule 1 and always in scope. Building convergence studies,
cross-route comparisons or oracles on your own initiative is not. If the
checks available look too weak to support what the task claims, say so in the
report; that is never a reason to present a weak check as validation.
