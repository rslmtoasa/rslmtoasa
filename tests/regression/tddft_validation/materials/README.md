# TDDFT-11 material campaign decks

These are the reproducible current-branch smoke decks used to compare the
three transverse `chi0` providers at selected `(q, omega)` points. They use
the committed Fe and Ni restart potentials from `results/validation/VAL-18_bccFe`
and `results/validation/VAL-19_fccNi`; no ground-state file is modified.

The three decks for each material differ only in `chi0_backend` and output
prefix, except that the native real-space deck emits bare `chi0` only because
the current provider has no validated exact static kernel for Dyson output.
All q coordinates are direct reciprocal coordinates and all energies are Ry.

Run from the material directory with the built executable, passing each deck
explicitly:

```sh
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_eigenpairs.nml
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_kspace_lehmann.nml
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_realspace_gf.nml
```

The release report records the exact commands and resulting backend/evidence
status. The checked-in `.dat` files are the raw outputs used by the machine
readable collector.
