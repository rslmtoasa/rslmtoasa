# LR-CLOSE-03 archived baselines

These tests close the vacuous `lr_baseline` label by running the current
executable and comparing the emitted physics rows with references generated
from the isolated pre-refactor archive at commit
`70eb73b4483d1099960fcfc014e23b8ee8576252`.

The archived decks are the predecessor forms of the checked-in current cases:

* `rotation`: old `post_processing='exchange_q'`, `rotation_dynamics=.true.`,
  native Turek reference, second-order `ham_only`, physical bcc Fe.
* `tddft_lehmann`: old `&tddft backend='product_lehmann'`, the bounded
  TDVK-02R2 compact Lehmann smoke seam.  The full radial-point Lehmann deck
  remains the existing `LinearResponseTddftSmoke` integration test; its
  complete Fe response space is 4,455 dense unknowns and is not a practical
  closeout smoke test.
* `native_rsgf`: old `&tddft backend='native_rsgf'`, block provider, three
  Green-function integration points.

Only numeric data rows in the requested output artifacts are compared.
Comments, build identity, timing, paths, log order, rotation channel text,
and status words are ignored.  The references are never regenerated from the
current executable.

Configure the regression tests and run the closeout set with:

```sh
cmake -S . -B build-p11-main -DRUN_REG_TESTS=ON
ctest --test-dir build-p11-main -N -L lr_baseline
ctest --test-dir build-p11-main --output-on-failure -L lr_baseline
```
