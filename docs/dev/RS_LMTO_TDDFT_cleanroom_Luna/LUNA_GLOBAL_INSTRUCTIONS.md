# Global Luna instructions

## Prime directive

Your goal is not to make TD-DFT work. Your goal is to implement only physics fixed by literature contracts and exact RS-LMTO basis mappings.

A BLOCKED task can be more successful than working code.

## Before each task

Read `00_MASTER_BLUEPRINT.md`, `00A_LITERATURE_CONTRACTS.md`, `00B_BLOCKER_POLICY.md`, and the task file.

Record:
- `git rev-parse HEAD`
- `git status --porcelain`

Do not start physics work from an unexplained dirty tree.

## Forbidden reasoning shortcuts

Do not use:
- "should be equivalent";
- "reasonable approximation";
- "probably reusable";
- "standard convention";
- "good enough";
- "for now use".

Provide a derivation/citation or declare BLOCKED.

## Old TD-DFT

After LR-00, old TD-DFT code is not an authority for signs, normalization, response basis, kernel, Goldstone or finite-q physics.

Do not copy it from Git history unless a task explicitly requests archaeology.

## Layout/naming

New files stay directly under `source/`.

Do not use author/institution names in modules, types, public options or filenames.

## No coding-agent physics choices

You may not decide:
- response projection;
- orbital/angular truncation;
- radial approximation;
- factors 2/4pi/muB;
- signs;
- Goldstone method;
- q convention;
- hybridization of published schemes.

Stop if not fixed by the task or completed mapping document.

## Tests

Tests must falsify physics contracts. Do not test one implementation only against another implementation of the same formula and call it independent.

Never weaken a deterministic tolerance merely to pass.

## Final task report

Return:
1. `PASS`, `BLOCKED`, or `FAIL`;
2. contracts touched;
3. source provenance;
4. equations/normalization;
5. changes;
6. tests and commands;
7. numerical evidence;
8. blockers/risks;
9. checklist;
10. commit hash/message if committed.
