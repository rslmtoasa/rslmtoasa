# C4 A4.2 — frozen native coefficient cross-swap

Both original contexts use `recur='block'`; swaps preserved each receiver's mesh, precision, window, input, and native route and stopped before any potential update. The two directions exchanged full production A, exact post-`zsqr` B, and all four native terminator outputs (`a_inf`, `b_inf`, `a_inf0`, `b_inf0`); readback was exact across 6,804 coefficient entries and 324 terminator matrix rows. Both swaps reproduce the opposite native d-LDOS and complex P_OUT/QL exactly at serialized precision, so the A4.2 evaluator-inconsistency condition is not met.

At A level k=0, the first off-diagonal difference is roundoff-sized (5.5e−17); at k=1, A(2,2) is −0.0176877658814 Ry (bulk) versus −0.0176877714694 Ry (impurity). The first k=1 B² and exact-B diagonal values differ by about 1.49e−8 Ry² and 1.53e−8 Ry, respectively. The full-array maxima are:

| Production coefficient | Compared extent | Maximum absolute bulk/impurity difference |
|---|---:|---:|
| A | 18×18×21 (k=0…20) | 2.7791e−6 Ry |
| B² | 18×18×21 (k=0…20) | 1.7285e−6 Ry² |
| exact post-`zsqr` B | 18×18×21 (k=0…20) | 2.1315e−6 Ry |

The historical C3 ~1e−16 d-only comparison used a diagnostic coefficient capture; it is distinct from the frozen production `PROD_COEFF` capture. Production uses the `cheb_backend=fast` FP32 recurrence, but these records do not isolate precision as the cause of the observed full-coefficient differences. Full paths, commands, input/source/executable and output hashes, capture bounds, per-run environments, and swap validation are in [C4_A42_PROVENANCE.json](C4_A42_PROVENANCE.json). No A4.1 ratio or overall Part A gate verdict is inferred here.

| Run | Donor coefficients and terminator | Result |
|---|---|---|
| `native_bulk_final` | Bulk | Frozen bulk reference |
| `native_imp_final` | Impurity | Frozen impurity reference |
| `swap_bulk_from_imp_final` | Impurity | Exact match to impurity native output |
| `swap_imp_from_bulk_final` | Bulk | Exact match to bulk native output |
