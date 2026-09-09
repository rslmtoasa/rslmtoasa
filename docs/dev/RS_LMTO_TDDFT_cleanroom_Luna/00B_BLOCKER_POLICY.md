# Blocker policy — points for stopping, not enabling

Every report statement must be tagged:
`[LITERATURE]`, `[BASIS MAPPING]`, `[IMPLEMENTATION]`, `[HYPOTHESIS]`, or `[BLOCKER]`.

Only the first three may enter production physics.

## Automatic blockers

- **B01** radial `f(r)` semantics/normalization not proved.
- **B02** Bxc available only through compressed LMTO potential parameters.
- **B03** ground-state XC and response-kernel XC provenance mismatch.
- **B04** unresolved definition of Bxc (`Vup-Vdn`, half splitting, sign, muB).
- **B05** unresolved radial measure/Jacobian.
- **B06** response-basis mapping requires undocumented orbital/site truncation.
- **B07** proposed replacement of radial/angular response by site scalar.
- **B08** old pair-Xi/kernel/Goldstone semantics reused without re-derivation.
- **B09** finite-q Fourier phase unresolved.
- **B10** unresolved Ry/Ha/muB units.
- **B11** Goldstone correction chosen empirically from Fe/Ni.
- **B12** published workflows mixed without exact equivalence proof.
- **B13** needed radial data discarded before response and not exactly recoverable.
- **B14** non-collinear transverse-only shortcut proposed in place of four-component response.

## On BLOCKED

1. Stop production edits.
2. Write `docs/<TASK>_BLOCKER_REPORT.md`.
3. Name the missing object/derivation precisely.
4. Name the nearest code/data evidence.
5. State what evidence would unblock it.
6. Do not implement a fallback.
