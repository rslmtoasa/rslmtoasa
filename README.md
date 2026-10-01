# RS-LMTO-ASA

**Electronic structure, magnetism, and transport in real space.**

[![Build](https://github.com/rslmtoasa/rslmtoasa/actions/workflows/binaries.yml/badge.svg)](https://github.com/rslmtoasa/rslmtoasa/actions/workflows/binaries.yml)
[![Example tests](https://github.com/rslmtoasa/rslmtoasa/actions/workflows/tests.yml/badge.svg)](https://github.com/rslmtoasa/rslmtoasa/actions/workflows/tests.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

RS-LMTO-ASA is an open-source electronic-structure code based on the **real-space linearized muffin-tin orbital method within the atomic sphere approximation**. It combines self-consistent electronic-structure calculations with recursion and Green-function techniques to study bulk materials, surfaces, and embedded impurities.

The code is written primarily in Fortran and supports OpenMP and MPI parallelization. Its real-space formulation connects local electronic structure with magnetic interactions and transport responses, including effects of spin–orbit coupling and noncollinear magnetism.

## Contents

- [Capabilities](#capabilities)
- [Build from source](#build-from-source)
- [Run your first calculation](#run-your-first-calculation)
- [Inputs and calculation modes](#inputs-and-calculation-modes)
- [Outputs and restarts](#outputs-and-restarts)
- [Parallel execution](#parallel-execution)
- [Examples](#examples)
- [Tests](#tests)
- [Contributing and support](#contributing-and-support)
- [Citing the code](#citing-the-code)
- [License](#license)

## Capabilities

- **Self-consistent electronic structure:** charge densities, total energies, densities of states, and local spin and orbital moments.
- **Magnetism:** collinear and noncollinear configurations, with scalar-relativistic or spin–orbit-coupled Hamiltonians.
- **Magnetic interactions:** exchange couplings, Dzyaloshinskii–Moriya interactions, and anisotropic exchange.
- **Transport:** charge, spin, and orbital conductivity calculations using a real-space Chebyshev formulation.
- **Geometry:** Bravais lattices, more general bulk structures, surfaces, and embedded impurity calculations.
- **Recursion:** Lanczos, block recursion, and Chebyshev methods, with availability depending on the calculation workflow.
- **Parallel computing:** OpenMP threads, MPI processes, and combined MPI/OpenMP execution.

The ASA represents the material using atomic spheres. Sphere radii, empty spheres where needed, and the real-space environment are part of the physical model and should be converged alongside the numerical parameters.

## Build from source

### Requirements

| Dependency | Purpose |
| --- | --- |
| CMake 3.19 or newer | Build configuration |
| Fortran compiler | Compile the electronic-structure code |
| C/C++ toolchain | Configure the top-level CMake project |
| BLAS and LAPACK | Linear algebra |
| Make or Ninja | Build execution |
| MPI implementation, optional | Distributed-memory parallelization |
| Python 3 and `requirements.txt`, for tests | Test runners and namelist handling |

GNU Fortran is a practical starting point. Compiler flags also include support for Intel compilers. Use a consistent compiler, MPI implementation, and linear-algebra library throughout a build.

### Linux: GNU Fortran and OpenBLAS

On Ubuntu or Debian, install the basic dependencies:

```bash
sudo apt update
sudo apt install git build-essential cmake gfortran libopenblas-dev liblapack-dev
```

Clone and build:

```bash
git clone https://github.com/rslmtoasa/rslmtoasa.git
cd rslmtoasa

cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_OPENMP=ON \
  -DENABLE_MPI=OFF

cmake --build build --parallel
```

The executable is **`build/bin/rslmto.x`**. This build uses one process and can use multiple OpenMP threads.

### MPI build

For GNU Fortran with Open MPI on Ubuntu or Debian:

```bash
sudo apt install openmpi-bin libopenmpi-dev

cmake -S . -B build-mpi \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_MPI=ON \
  -DENABLE_OPENMP=ON

cmake --build build-mpi --parallel
```

For Intel oneAPI, load the compiler, MPI, and math-library environment first, and select the matching MPI Fortran wrapper, such as `mpiifx`. BLAS/LAPACK selection may require a `BLA_VENDOR` setting appropriate to the installed libraries. Check CMake's detected compiler and library paths before building.

Use a separate build directory when changing compilers or MPI implementations.

### Build options

| Option | Default | Description |
| --- | --- | --- |
| `CMAKE_BUILD_TYPE` | `Release` | `Release`, `Debug`, or `Testing` |
| `ENABLE_OPENMP` | `ON` | Enable OpenMP when detected |
| `ENABLE_MPI` | `OFF` | Enable MPI; requires an MPI toolchain |
| `ENABLE_MARCH_NATIVE` | `ON` | Request host CPU optimization |
| `COLOR` | `ON` | Enable colored terminal output |
| `ENABLE_FLUSH` | `OFF` | Flush printouts |
| `RUN_EXAMPLE_TESTS` | `OFF` | Register SCF and post-processing example tests |
| `RUN_REG_TESTS` | `OFF` | Register the legacy recursion regression tests |

For GNU builds targeting another CPU, set `-DENABLE_MARCH_NATIVE=OFF`. CPU flag selection is compiler-dependent; check the resolved flags when preparing portable binaries.

For installation after building:

```bash
cmake --install build --prefix "$HOME/.local"
```

## Run your first calculation

From the repository root, copy the supplied bcc Fe example into a working directory:

```bash
mkdir -p runs
cp -r example/bulk/bccFe runs/bccFe
cd runs/bccFe

export OMP_NUM_THREADS=4
../../build/bin/rslmto.x input.nml > run.log 2>&1
```

Run the executable **inside the calculation directory**, where its atom files and other inputs are available. If no argument is supplied, the program reads `input.nml` by default.

The bcc Fe example sets `pre_processing = 'bravais'`, block recursion, and up to 100 SCF steps. It also includes saved atom states, making it useful for checking a build with an existing electronic-structure starting point.

Inspect `run.log`, `report.out`, and the resulting atom files. Reaching the end of a run does not by itself establish SCF convergence: check the reported convergence behavior and your chosen tolerances.

## Inputs and calculation modes

Calculations use **Fortran namelists**. Start from the closest supplied example and adapt its geometry, atom labels, magnetic configuration, and numerical settings.

| Input | Role |
| --- | --- |
| `input.nml` | Main calculation settings |
| `lattice.nml`, when required | Geometry for workflows using a separate lattice definition |
| `<label>.nml` | Atom or symbolic-atom input associated with the labels in `&atoms` |
| `<label>_out.nml`, when present | Saved atom state used by the restart loader |
| Additional structure files | Geometry-dependent inputs, such as an existing cluster for an impurity calculation |

Labels identify symbolic atom types and may distinguish inequivalent sites of the same element.

### Main namelists

| Namelist | Controls |
| --- | --- |
| `&calculation` | Workflow selection and verbosity |
| `&lattice` | Lattice, real-space environment, and exchange-pair settings |
| `&atoms` | Atom labels and input database location |
| `&control` | Calculation type, spin treatment, and recursion settings |
| `&self` | SCF iterations, convergence settings, and SOC scaling |
| `&energy` | Fermi energy, energy range, and energy sampling |
| `&mix` | Mixing algorithm and parameters |
| `&hamiltonian` | Hamiltonian settings and current operators for transport |

The accepted variables are defined in [`source/include_codes/namelists/`](source/include_codes/namelists/). Defaults and additional behavior are implemented in the corresponding source modules.

### Geometry and workflow selection

Common `&calculation` selectors are:

| Setting | Workflow |
| --- | --- |
| `pre_processing = 'bravais'` | Bravais-lattice setup and SCF |
| `pre_processing = 'newclubulk'` | General bulk/cluster setup and SCF |
| `pre_processing = 'buildsurf'` | Surface construction and SCF |
| `pre_processing = 'newclusurf'` | Surface/cluster setup and SCF |
| `post_processing = 'exchange'` | Magnetic interaction calculations |
| `post_processing = 'conductivity'` | Conductivity calculations |

In `&control`, `calctype = 'B'`, `'S'`, and `'I'` select bulk, surface, and impurity treatments, respectively. Geometry and workflow selectors must be consistent with the supporting input files.

For exchange calculations, the supplied example defines `njij` and `ijpair` in `&lattice`. Pair indices refer to the generated real-space cluster; verify them after changing the structure.

### Spin treatment

The `nsp` parameter in `&control` selects the magnetic and relativistic treatment:

| `nsp` | Treatment |
| --- | --- |
| `1` | Collinear, scalar relativistic |
| `2` | Collinear, with spin–orbit coupling |
| `3` | Noncollinear, scalar relativistic |
| `4` | Noncollinear, with spin–orbit coupling |

**`nsp = 1` does not mean nonmagnetic.** It selects the collinear scalar-relativistic treatment.

### Numerical convergence

Check convergence with respect to the real-space environment, recursion depth, energy range and sampling, and SCF mixing. For Chebyshev calculations, the energy interval must cover the relevant Hamiltonian spectrum. Transport and exchange observables can require stricter convergence than total energies or local moments.

## Outputs and restarts

Files produced depend on the selected workflow and output settings.

| Output | Contents |
| --- | --- |
| `run.log` | Terminal output captured by the shell command above |
| `report.out` | Summary including energies, moments, and magnetic forces |
| `totaldos.out` | Density-of-states output |
| `<label>_scf.nml` | Atom state saved during SCF |
| `<label>_out.nml` | Saved atom state for subsequent calculations |
| `jij.out`, `jijso.out`, `jijfo.out` | Exchange outputs, depending on the evaluation routine |
| `dijso.out`, `dijfo.out` | DMI outputs from second- and first-order evaluations |
| `cond_total.out` | Total conductivity output |
| `<label>_cond.out` | Conductivity contributions for per-type evaluation |

Column definitions, normalization, and units depend on the output and calculation path. Consult the writing routine before converting results to a spin Hamiltonian or comparing transport responses across geometries. The `so` and `fo` suffixes in the exchange files identify perturbation order.

### Restart from an SCF checkpoint

The restart loader reads `<label>_out.nml` files when they exist. To resume from saved `_scf.nml` states, first back up any existing `_out.nml` files, then run this in the calculation directory:

```bash
for file in *_scf.nml; do
  [ -f "$file" ] || continue
  cp -- "$file" "${file%_scf.nml}_out.nml"
done
```

Run the calculation again with the same geometry and compatible atom labels. These files provide saved atom states; they do not replace the other required inputs.

## Parallel execution

For an OpenMP build:

```bash
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close
/path/to/rslmto.x input.nml > run.log 2>&1
```

For an MPI/OpenMP build:

```bash
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close
mpirun -np 10 /path/to/rslmto.x input.nml > run.log 2>&1
```

The second example requests 10 MPI processes with 4 threads each, requiring 40 allocated CPU cores. Match the process/thread layout to the scheduler allocation and use the MPI runtime associated with the build. Benchmark representative calculations to choose the layout for your system.

## Examples

| Directory | Examples |
| --- | --- |
| [`example/bulk/`](example/bulk/) | bcc Fe, Si, graphene, and Mn₃Sn |
| [`example/surface/`](example/surface/) | Cu(001) surface |
| [`example/impurity/`](example/impurity/) | B2 FeCo impurity setup |
| [`example/exchange/`](example/exchange/) | Exchange calculations for bcc Fe |
| [`example/conductivity/`](example/conductivity/) | Charge, spin, and orbital transport setups for Fe, Cu, Pt, and PtSe₂ |

Some examples include converged or partially converged atom states. Copy the whole example directory to preserve its supporting files, and review its settings before using it as a production calculation.

## Tests

From the repository root, install the Python dependencies and register the example tests:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

cmake -S . -B build \
  -DRUN_EXAMPLE_TESTS=ON \
  -DEXAMPLE_PYTHON_EXECUTABLE="$VIRTUAL_ENV/bin/python"
cmake --build build --parallel

ctest --test-dir build -L example --output-on-failure
```

To run one suite:

```bash
ctest --test-dir build -L scf --output-on-failure
ctest --test-dir build -L postproc --output-on-failure
```

Example tests run in isolated directories under `build/Testing/` and compare selected outputs with stored references where available. Cases without numerical references run as smoke tests. Short SCF tests check implementation behavior; they do not establish production convergence.

See the [SCF test guide](tests/scf/README.md) and [post-processing test guide](tests/postproc/README.md) for case definitions, tolerances, and reference generation.

## Contributing and support

Bug reports, improvements, examples, and documentation contributions are welcome through [GitHub Issues](https://github.com/rslmtoasa/rslmtoasa/issues) and [pull requests](https://github.com/rslmtoasa/rslmtoasa/pulls).

For a reproducible bug report, include:

- The commit or branch used.
- Compiler, CMake, BLAS/LAPACK, and MPI versions.
- Build options and the command used to launch the calculation.
- A minimal input set, relevant logs, and the expected behavior.

For numerical changes, describe the physical quantity affected and validate the change against a suitable example or reference calculation.

## Citing the code

When publishing results, identify RS-LMTO-ASA, link to this repository, and record the release or commit used. Also cite the methodological publications relevant to your calculation, including the electronic-structure approach and the exchange or transport formalism employed.

For an up-to-date citation recommendation, contact the maintainers through the repository.

## License

RS-LMTO-ASA is distributed under the **GNU General Public License, version 3**. See [`LICENSE`](LICENSE) for the full terms.
