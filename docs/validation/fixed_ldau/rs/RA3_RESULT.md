# RA3 result: fixture and spin guards

**Measured: pass.** The nsp=1 bcc Fe bulk and B2FeCo impurity fixed-U decks
completed one-step runs. The retained nsp=2 bulk and impurity decks still
fail with the existing fixed-U SOC/noncollinear guard. Adding
`hubbard_radial_moments='local'` to those nsp=2 decks separately reaches and
fails the required `requires nsp=1` radial-mode guard before the SOC check.

This preserves both guard checks without claiming nsp=2 is a supported local
fixed-U fixture. Inputs, logs, and hashes are in
[`RA3_PROVENANCE.json`](RA3_PROVENANCE.json); raw outputs remain in
`/tmp/rslmto-ra3-verified/`.

The VAL-04 driver was also corrected to preserve its fixture's `nsp=1` instead
of rewriting it to `nsp=2`. Its complete five-run k-space sweep passed against
the final-source executable without rebuilding. The exact command, executable
and runner hashes, and per-case output hashes are recorded in
[`RA3_PROVENANCE.json`](RA3_PROVENANCE.json); raw runs remain in
`/tmp/val04-ra3-nsp1/`. This driver-level check complements the one-step RA3
fixture and guard cases above.
