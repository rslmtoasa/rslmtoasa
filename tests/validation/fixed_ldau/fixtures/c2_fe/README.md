# C2 Fe bulk/impurity integration fixtures

These small input files let `c2_operator_mapping_regression.py` exercise the
production bulk and local Hubbard insertion routines with one accepted Fe
starting potential. The bulk and impurity decks are copied from the matched
one-step C1 Fe validation decks; the runner changes only their `&hamiltonian`
groups to create the listed U/J controls. It writes fresh potential inputs in
each scratch case and never modifies these fixture files.

`Fe_plus_u_seed.nml` is the accepted C1 Fe bulk U=3.8 eV/J=1.1 eV endpoint
potential, copied byte-for-byte from
`/private/tmp/rslmto-ra4-frozen/c1_qlmode_continuations/fe_bulk_u38j11_beta0015/Fe_out.nml`.
Its SHA256 is `d4abdc888fa8ee2d64deb71ac9feb81fd0c960fe7749b7cb5cb2f3ff754bbe38`.
The runner derives a U=J=0 seed for the D2 controls by changing only the stored
`hubbard_u` and `hubbard_j` fields; those transformed seeds are recorded in
each run's `PRELAUNCH_SHA256SUMS.txt`.

Fixture SHA256 values:

| File | SHA256 |
| --- | --- |
| `input_bulk.nml` | `3aeb8b5e5570506095fd1508147c66b6167a1209454f55ed5edbc86ed83ffc4b` |
| `input_impurity.nml` | `cbf421d32d61ef143a01eb56d7bbcc53a816f688c660634fef0269de62e5597e` |
| `Fe_plus_u_seed.nml` | `d4abdc888fa8ee2d64deb71ac9feb81fd0c960fe7749b7cb5cb2f3ff754bbe38` |

From the repository root, run with the frozen audit-enabled executable:

```sh
python3 tests/validation/fixed_ldau/c2_operator_mapping_regression.py \
  --binary /path/to/rslmto.x \
  --workdir /tmp/c2-operator-mapping-run
```

The runner requires a fresh work directory and preserves all generated inputs,
logs, traces, audit records, and hashes there. It checks D1 same-material U/J
block and P_OUT agreement; D2(a) J=0; D2(b) host-only U; and D2(c) impurity-only
U through both general and impurity-override namelist fields. The C2 audit
environment switch is enabled, and the safe moments-only C1 trace is used to
capture P_OUT without invoking the unrelated legacy bond-trace path.

The same impurity deck also exercises the production residual-mapping bounds
regression. Use a Debug executable built with Fortran bounds checking:

```sh
python3 tests/validation/fixed_ldau/c2_residual_mapping_bounds_smoke.py \
  --binary /path/to/rslmto-debug-bounds.x \
  --workdir /tmp/c2-residual-mapping-bounds
```

This isolated two-step diagnostic asserts that the generated `ntype=2` cluster
has a physical representative index above `ntype`, that the production helper
completes twice, and that recursion site 1 is reported as impurity type 2. It
does not test convergence or serve as a physics acceptance record.
