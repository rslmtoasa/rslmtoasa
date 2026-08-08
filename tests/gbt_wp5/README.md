# GBT WP5 route-unification probes

These fixed-potential bcc-Fe cases use the same `gbt_single_q` primitive-bond
builder for both solver routes. The reciprocal cases vary the full uniform
k mesh; the real-space cases vary block-recursion depth. Every run writes
`gbt_bonds.out`, which must be byte-identical across routes for the same q.

The convergence table and the exact commands/results accepted at G5 are
recorded in `docs/dev/GBT_WP5_G5_REPORT.md`.
