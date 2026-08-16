# VAL-14 — Cu(111) layered interface identity mapping

## Scope and result

This validation establishes the Cu(111) layered/interface identity at the
existing pinned Fermi level and Chebyshev mesh. A single active Cu layer
between identical Cu reference regions must reproduce the Cu(001) identity
control. The Hamiltonian stage is accepted from the existing geometry-keyed,
bit-identical H(R) comparison across bulk, impurity, and layered routes.

The residual is **resolved**. The fix is generic to layered geometry mapping;
it does not rescale DOS values or branch on a Miller index. The broader
`calctype='L'` interface family remains **Experimental** because the separate
charge-row alignment behavior is still open.

## Trace and root cause

The first differing downstream quantity was the selected layered cluster:

| Stage | Cu(001) before fix | Cu(111) before fix | VAL-14 result |
|---|---:|---:|---:|
| source/selected cluster sites | 4096 | 4056 | 4096 for both |
| representative active site | valid | valid but in clipped cluster | valid |
| real-space `nn` map | complete | built from clipped geometry | complete |
| recursion/GF/projection | downstream consumers | inherited the map difference | identity restored |

`build_interface_full` constructed a fixed z ladder. Its bounds covered the
full source cluster for (001), but the (111) projection range extended beyond
the ladder (`-21 … 24` versus ladder `-16 … 22`), silently dropping 40 sites.
The ladder now derives lower and upper bounds from the projected source
cluster, on the same z-step grid and with the existing safety margin. The
active-layer numbering remains unchanged.

## Evidence

With the corrected mapping, the Cu(111) and Cu(001) `Ac1_dos.out` and
`totaldos.out` files are byte-identical for the sampled mesh. At the former
maximum residual (DOS row 1200), the corrected values are:

- site DOS: `12.56516`
- total DOS: `12.60931`
- `vmad`: `-4.22e-15 Ry`

The interface cluster contains 4096 sites and the active representative has
the complete 19-neighbour fcc shell. The bulk and impurity control fixtures
remain passing. Interface electrostatics remains active on the unchanged
`calctype='L'` dispatch.

## Checklist

- [x] Existing H(R) equivalence accepted/reconfirmed
- [x] First downstream divergence located
- [x] Geometry/layer mapping traced
- [x] Generic root cause fixed
- [x] No 111 special case added
- [x] Cu001 identity retained
- [x] Cu111 identity restored
- [x] Surface fixture unchanged
- [x] Interface electrostatics unchanged and still active

## Changed files and validation

- `source/lattice_cluster.f90`
- `tests/scf/references/Example_interface_fccCu111_chebyshev/ref.json`
- `tests/scf/references/Example_interface_fccCu111_chebyshev/ref.macos-arm64.json`
- `tests/KNOWN_ISSUES.md`
- `docs/dev/PHASE_II_VALIDATION.md`
- this report

Focused tests:

```text
ctest --test-dir build-rf-debug -R 'Example_(bulk_fccCu|impurity_fccCu|interface_fccCu001|interface_fccCu111)_chebyshev' --output-on-failure
```

The three controls pass; Cu(111) passes against the corrected identity
reference. Additional surface, interface-alignment, and vacuum wiring checks
are run as part of the VAL-14 handoff.

Commit message: `Fix Cu111 layered interface identity mapping`
