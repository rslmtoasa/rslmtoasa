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

## The five rules that matter most

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
