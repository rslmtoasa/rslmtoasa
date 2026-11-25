# Hamiltonian export validation

Build with CMake and `-DRUN_UNIT_TESTS=ON`. Run `ctest -R UnitEffectiveHamiltonian -V` for the deterministic mixed local/bulk algebra and wrapped-image tests. The run directory for `HamiltonianActionProbe` must contain a prepared frozen deck. Set `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` for the recorded comparisons.

In a deck with `post_processing='hamiltonian_export'`, run `rslmto.x` first to produce the canonical files. Then run `HamiltonianActionProbe native` and `HamiltonianActionProbe import` in the same directory. The test probe uses three deterministic complex full-support vectors and the unchanged native unscaled H and V actions. Run the independent checker as follows.

```shell
python3 tests/hamiltonian_export/check_exports.py <deck> --action
python3 tests/hamiltonian_export/check_exports.py <cluster-deck> --cluster-file
python3 tests/hamiltonian_export/check_exports.py <periodic-deck> --kspace <reciprocal-deck>
python3 tests/hamiltonian_export/check_exports.py <special-kspace> --native-eigenpairs <native-fermi-surface>
python3 tests/hamiltonian_export/check_exports.py <import-deck> --transport <native-deck> --baseline <pinned-baseline-deck>
```

Run the FeO installer with the attached archive. It checks the SHA-256, saves the original archive contents separately and creates four derived decks. The native structure builder regenerates screened constants. The four archived restart potentials and the two DM sidecars are copied byte for byte. Their embedded backend is `rs`, and the input explicitly selects local radial moments. The current `hubbard_backend` selector is for SCF, so the frozen builder reads the backend provenance in the saved potentials.

```shell
python3 tests/hamiltonian_export/install_feo.py 'feou4hoh (1).tar' feo-validation
```

The checker requires NumPy and SciPy. The FeO installer requires f90nml.

The standard CTest cases `EffectiveHamiltonian-round-trip` and
`EffectiveHamiltonian-block-dos` compare the native and imported results on
a mixed-boundary Cu replica. They cover canonical and seven-column replay,
fixed Fermi references, damaged payload and metadata rejection, open-boundary
translations, and block DOS. They exercise export under block, scalar Lanczos
and Chebyshev inputs and compare the resulting matrix bytes.

`EffectiveHamiltonian-magnetic` compares bcc Fe exchange with onsite and
nearest-neighbor pairs through native and imported Chebyshev and block
recursion. The same HOH/SOC export is used for two Kubo operator slots with
spin torque and SOC torque. `EffectiveHamiltonian-surface` compares orbital
accumulation in a Cu surface with frozen and active atom types. Configure
with `-DRUN_EXAMPLE_TESTS=ON` to register these CLI tests.
