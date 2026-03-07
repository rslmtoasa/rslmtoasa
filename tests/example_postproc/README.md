Example post-processing tests (exchange, conductivity)
===============================================

- Purpose: These are fast smoke-style checks for post-processing workflows (exchange, conductivity). They run a small, deterministic input and confirm the program runs without fatal errors.
- Location: test cases live under `tests/example_postproc/cases`.

Smoke vs regression
-------------------
- By default these are smoke tests: they are registered in CTest and will run quickly to detect runtime failures. They do not perform full numerical regression comparisons unless you explicitly generate and enable references.

Generating reference data (regression)
-------------------------------------
1. Use the included generator to create reference `data.nml` outputs from a known-good binary. From the repository root run:

```bash
tests/example_postproc/generate_reference_data.sh /path/to/known-good/rslmto.x [--include-cheb] [--scratch-root <path>]
```

- This will create `tests/example_postproc/references/<TestName>/ref.nml` and a `meta.txt` with the run parameters.

2. Enable reference-based tests when configuring CMake by setting:

```bash
-DRUN_EXAMPLE_TESTS=ON -DRUN_EXAMPLE_REF_TESTS=ON -DEXAMPLE_REF_PYTHON_EXECUTABLE=/path/to/python
```

- The CTest entries for the postproc tests will then include `_ref` variants that run the comparator using `f90nml` (make sure `f90nml` is installed in the Python you point to).

Notes
-----
- Postproc reference outputs are separate from `tests/example_scf/references` and are stored under `tests/example_postproc/references`.
- Chebyshev-based runs are registered but disabled by default; use `--include-cheb` when generating references and set `-DRUN_EXAMPLE_CHEB_TESTS=ON` to enable them in CTest.
