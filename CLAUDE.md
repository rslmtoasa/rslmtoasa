# RS-LMTO-ASA — orientation for a fresh AI session

Start here, then read:

- **[`docs/DEVELOPER_MAP.md`](docs/DEVELOPER_MAP.md)** — entry points, class
  chains, kernel inventory, testing map. The fastest way to find where a
  workflow lives.
- **[`docs/dev/REFACTORING_PLAN.md`](docs/dev/REFACTORING_PLAN.md)** — Phase 1
  ground rules and task history (structural refactor, now complete).
- **[`REFACTORING_PHASE2.md`](REFACTORING_PHASE2.md)** — current phase: test
  coverage, CI, documentation. Check its progress checklist before starting
  work.
- **[`tests/README.md`](tests/README.md)** — suite index, CI trigger/matrix
  strategy, GPU coverage tiers.
- **[`tests/KNOWN_ISSUES.md`](tests/KNOWN_ISSUES.md)** — bugs found via
  coverage work, deliberately left unfixed pending a dedicated task.

## The six rules that matter most

1. **Bit-level behavior is the contract.** `tests/regression` and the
   `tests/{scf,postproc}` example suites must pass at the same tolerances
   before and after every change. Run them before starting and after every
   task.
2. **One task, one commit.** Don't batch unrelated changes into one commit.
3. **KISS.** No argument-validation ladders, wrapper layers, or "just in
   case" guards on internal routines — checks belong at true boundaries
   (namelist parsing, file I/O, plugin availability).
4. **Class-based architecture stays.** Derived types, constructors,
   `restore_to_default`, type-bound procedures. New code follows the same
   pattern.
5. **`self.f90` and `symbolic_atom.f90` are off-limits.** Physics-dense
   legacy code scheduled for a separate audit — edit around the edges only.
6. **Stay inside the task Anders defined.** See below — this one is about
   token budget, and it is not optional.

## Scope discipline — do only the task that was asked for

Anders defines the task. Work on that, finish it, report, and **stop**. Do not
extend the scope because an adjacent problem looks tractable or interesting.
This is a budget constraint, not a limit on capability: exploratory work burns
tokens that were not authorized for it.

**Green-light rule.** Anything beyond the stated task needs Anders' explicit
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

**Corollary — don't over-verify either.** Running the regression suite before
and after is rule 1 and always in scope. Building fresh convergence studies,
cross-route comparisons, or new oracles to satisfy your own curiosity is not,
unless the task asked for validation.
