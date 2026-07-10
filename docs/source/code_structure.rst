.. _code_structure:

=======================================
Code Structure
=======================================

.. note::

   For a task-oriented developer cheat-sheet (entry points, per-workflow call
   chains, kernel inventory, "where to add X" recipes, and the testing map),
   see ``docs/DEVELOPER_MAP.md`` in the repository root. This page gives the
   higher-level directory and subsystem overview; the developer map is the
   fastest way to find where a specific workflow lives.

Directory Layout
================

.. code-block:: text

   rslmto_fable/
   ├── source/                       # Fortran source code
   │   ├── main.f90                   # Program entry point
   │   ├── calculation.f90            # Top-level driver: pre/processing/post select-cases
   │   ├── control.f90                # Control parameters and settings
   │   ├── self.f90                   # Self-consistent field (SCF) loop core (legacy, off-limits)
   │   │
   │   ├── lattice.f90                # Crystal structure & geometry — parent module (type + interfaces)
   │   ├── lattice_lifecycle.f90      #   submodule: constructor / restore_to_default / build_from_file
   │   ├── lattice_cluster.f90        #   submodule: real-space cluster construction (newclu, surfaces)
   │   ├── lattice_strux.f90          #   submodule: structure-constant backend dispatch
   │   ├── lattice_print.f90          #   submodule: I/O and reporting
   │   │
   │   ├── hamiltonian.f90            # TB-LMTO Hamiltonian — parent module (type + interfaces)
   │   ├── hamiltonian_build.f90      #   submodule: build_bulkham / build_locham
   │   ├── hamiltonian_ccor.f90       #   submodule: combined-correction terms
   │   ├── hamiltonian_hubbard.f90    #   submodule: +U / +V (Hubbard) corrections
   │   ├── hamiltonian_paoflow_io.f90 #   submodule: PAOFLOW-format import/export
   │   │
   │   ├── recursion.f90              # Recursion method — parent module (type + interfaces)
   │   ├── recursion_core.f90         #   submodule: Lanczos (recur / crecal)
   │   ├── recursion_haydock.f90      #   submodule: block-Lanczos (recur_b / crecal_b)
   │   ├── recursion_chebyshev.f90    #   submodule: Chebyshev moment recursion
   │   ├── recursion_transport.f90    #   submodule: intersite / transport variants
   │   ├── chebyshev_fast.f90         # Fast Chebyshev CPU kernels (fast/batched/MKL)
   │   ├── haydock_fast.f90           # Fast block-Lanczos CPU kernels
   │   ├── rsrec_cuda_plugin.f90      # CUDA plugin wrapper (iso_c_binding); kernels in source/cuda/
   │   │
   │   ├── green.f90                  # Green's function calculations
   │   ├── energy.f90                 # Energy mesh and Fermi level
   │   ├── charge.f90                 # Charge density and electrostatics
   │   ├── mix.f90                    # Density mixing strategies
   │   ├── bands.f90                  # Band moments / magnetic torques
   │   ├── density_of_states.f90      # Real-space DOS calculations
   │   │
   │   ├── reciprocal.f90             # k-space (reciprocal) route — parent module
   │   ├── reciprocal_lifecycle.f90   #   submodule: construction / lifecycle
   │   ├── reciprocal_fourier.f90     #   submodule: H(k) via Fourier transform of hoppings
   │   ├── reciprocal_bands.f90       #   submodule: diagonalization, band structure
   │   ├── reciprocal_dos.f90         #   submodule: tetrahedron / gaussian DOS integration
   │   ├── reciprocal_projection.f90  #   submodule: orbital-projected DOS, band moments
   │   │
   │   ├── exchange.f90               # Heisenberg exchange (J_ij / D_ij)
   │   ├── conductivity.f90           # Real-space transport / conductivity tensor
   │   ├── spin_dynamics.f90          # Atomistic spin dynamics (LLG integrator)
   │   ├── xc.f90                     # Exchange-correlation functionals
   │   ├── potential.f90              # Potential parameters
   │   ├── element.f90                # Element properties
   │   ├── symbolic_atom.f90          # Atomic site information (legacy, off-limits)
   │   ├── basis.f90                  # Basis-set handling
   │   ├── symmetry.f90               # Symmetry operations
   │   ├── spglib_interface.f90       # spglib bindings
   │   ├── sparse.f90                 # Sparse-matrix storage/utilities
   │   ├── namelist_generator.f90     # Input namelist utilities
   │   ├── report.f90                 # Run report / summary output
   │   │
   │   ├── precision.f90              # Floating-point precision definitions
   │   ├── math.f90                   # Mathematical utilities
   │   ├── string.f90                 # String manipulation
   │   ├── logger.f90                 # Logging and debugging
   │   ├── timer.f90                  # Performance timing
   │   ├── mpi.f90                    # MPI wrappers
   │   ├── safe_alloc.f90             # Memory tracking
   │   ├── globals.f90                # Global constants
   │   ├── os.f90                     # OS utilities, command-line argument parsing
   │   ├── lists.f90                  # List and array utilities
   │   ├── array.f90                  # Advanced array operations
   │   ├── bounds.f90                 # Bounds / range helpers
   │   ├── face.F90                   # Preprocessor-facing helper macros
   │   ├── cuda/                      # CUDA device kernels + C API (rsrec_gpu.cu, rsrec_cuda.h)
   │   └── include_codes/             # Vendored/included code (strux_lib, abspinlib, namelists, ...)
   │
   ├── cmake/                         # CMake configuration modules
   │   ├── SetFortranFlags.cmake      # Compiler-specific flags
   │   ├── SetCompileFlag.cmake       # Compilation settings
   │   └── git-watcher.cmake          # Version tracking
   │
   ├── tests/                         # Test suites
   │   ├── regression/                # Bit-level regression suite (backend matrix)
   │   ├── scf/                       # Bulk/surface/impurity SCF examples
   │   ├── postproc/                  # Exchange, conductivity, bands, DOS, PAOFLOW import
   │   ├── run_regression_tests.sh
   │   └── generate_references.py
   │
   ├── docs/                          # Sphinx documentation (this)
   │   ├── source/
   │   │   ├── index.rst
   │   │   ├── theory/                # Theoretical foundations
   │   │   ├── user_guide/            # User documentation
   │   │   ├── keywords/              # Input parameter reference
   │   │   └── reference/             # Code reference
   │   ├── conf.py
   │   └── DEVELOPER_MAP.md           # Task-oriented developer cheat-sheet
   │
   ├── CMakeLists.txt                 # Root CMake configuration
   ├── source/CMakeLists.txt          # Source file compilation rules
   ├── CLAUDE.md                      # Orientation for AI-assisted sessions
   └── README.md                      # Quick start guide

.. note::

   Several formerly monolithic modules were split into a **parent module +
   Fortran submodules** during the structural refactoring. The parent file
   (e.g. ``recursion.f90``) declares the derived type and the interfaces of
   its type-bound procedures; the ``<parent>_<name>.f90`` files
   (e.g. ``recursion_chebyshev.f90``) are ``submodule (recursion_mod)`` units
   that carry the implementations. This splits *implementation* across files
   without flattening the class — the public type and its interface stay in
   one place. The families are ``lattice_*``, ``hamiltonian_*``,
   ``recursion_*``, and ``reciprocal_*``.

Build System (CMake)
====================

**Main Files:**

- ``CMakeLists.txt`` - Root configuration; defines project (``rslmto``),
  compiler, options, and the ``rslmto.x`` executable target
- ``source/CMakeLists.txt`` - Source file list and per-target options
  (MKL kernels, CUDA plugin, strux_lib, abspinlib)
- ``cmake/SetFortranFlags.cmake`` - Compiler detection and flags (gfortran, ifort, etc.)
- ``cmake/SetCompileFlag.cmake`` - Additional compilation rules

The core sources compile into a static library (``rslmto``) which is linked
into the ``rslmto.x`` executable (built from ``source/main.f90``). Optional
components are gated by CMake options: ``ENABLE_MKL_KERNELS`` (MKL Chebyshev
kernels) and the CUDA plugin (``USE_CUDA_PLUGIN``).

**Key Build Targets:**

.. code-block:: bash

   make              # Compile the code (produces rslmto.x)
   make html         # Build HTML documentation (from docs/)
   make install      # Install executable to <prefix>/bin
   make clean        # Remove build artifacts

Module Dependency Graph
=======================

**Core Dependencies:**

.. code-block:: text

   main.f90
   └── calculation.f90 (top-level driver)
       │   pre_processing  select-case -> pre_processing_*  (bravais/buildsurf/newclu*)
       │   processing      select-case -> processing_sd
       │   post_processing select-case -> post_processing_* (exchange/conductivity/
       │                                   paoflow2rs/orbital_modern/bands/dos/...)
       ├── control.f90
       ├── lattice.f90 (+ lattice_* submodules)
       ├── charge.f90
       ├── self.f90 (SCF loop core)
       │   ├── hamiltonian.f90 (+ hamiltonian_* submodules)
       │   ├── recursion.f90 (+ recursion_* submodules, chebyshev_fast, haydock_fast, cuda plugin)
       │   ├── green.f90
       │   ├── bands.f90
       │   ├── density_of_states.f90
       │   └── mix.f90
       ├── reciprocal.f90 (+ reciprocal_* submodules)   # k-space route (bands / DOS)
       ├── exchange.f90
       ├── conductivity.f90
       ├── spin_dynamics.f90
       └── ... (output / report modules)

**Utility Modules (used everywhere):**

- ``precision.f90`` - Real/complex number precision
- ``math.f90`` - Mathematical constants and functions
- ``string.f90`` - String utilities
- ``logger.f90`` - Logging and output
- ``timer.f90`` - Performance timing
- ``mpi.f90`` - MPI wrappers
- ``safe_alloc.f90`` - Allocation tracking
- ``globals.f90`` - Global constants and parameters

Major Subsystems
================

**1. Initialization & Control** (``control.f90``, ``os.f90``, ``namelist_generator.f90``)

Reads and validates input parameters from Fortran namelist files. ``calculation.f90``
dispatches calculation modes through three independent select-cases
(``pre_processing`` / ``processing`` / ``post_processing``).

**2. Structure & Geometry** (``lattice.f90`` + ``lattice_*`` submodules, ``symbolic_atom.f90``, ``element.f90``, ``symmetry.f90``, ``spglib_interface.f90``)

Handles crystal structure definition, atomic positions, Bravais lattice,
real-space cluster construction, structure-constant backend dispatch, and
element-specific parameters.

**3. Basis & Hamiltonian** (``potential.f90``, ``basis.f90``, ``hamiltonian.f90`` + ``hamiltonian_*`` submodules)

Constructs the tight-binding Hamiltonian matrix from TB-LMTO parameters and
potential functions. Supports spin-orbit coupling, magnetic moments,
combined-correction terms (``hamiltonian_ccor``), Hubbard +U/+V
(``hamiltonian_hubbard``), and PAOFLOW import/export
(``hamiltonian_paoflow_io``).

**4. Electronic Structure — recursion route** (``recursion.f90`` + ``recursion_*`` submodules, ``green.f90``)

Core real-space computational engine:

- **Recursion method**: Lanczos (``recursion_core``), block-Lanczos
  (``recursion_haydock``), and Chebyshev (``recursion_chebyshev``) recursion.
  Hot kernels live in ``chebyshev_fast.f90`` / ``haydock_fast.f90``, with an
  optional CUDA path (``rsrec_cuda_plugin.f90`` + ``source/cuda/``).
- **Green's functions**: on-site and inter-site Green's functions at complex
  energies, using recursion coefficients to avoid explicit diagonalization.

**5. Electronic Structure — k-space route** (``reciprocal.f90`` + ``reciprocal_*`` submodules)

A parallel, diagonalization-based path for band structure and DOS. Builds
``H(k)`` by Fourier transform of the real-space hoppings
(``reciprocal_fourier``), diagonalizes (``reciprocal_bands``), and integrates
the DOS via tetrahedron or gaussian methods (``reciprocal_dos``). This route
does **not** go through ``recursion``/``green``.

**6. Self-Consistent Field** (``self.f90``, ``mix.f90``, ``charge.f90``)

Iterative SCF loop: computes charge density from Green's functions, updates
effective potentials, applies density mixing (linear, Broyden, ...), and
checks convergence.

**7. Properties Calculation** (``density_of_states.f90``, ``bands.f90``, ``exchange.f90``, ``conductivity.f90``, ``spin_dynamics.f90``)

Post-processing / driven modules:

- **DOS / Bands**: local and total density of states, band moments
- **Exchange**: Heisenberg J_ij / D_ij parameters (for magnets)
- **Conductivity**: real-space transport coefficients
- **Spin dynamics**: atomistic LLG time integration

**8. Exchange-Correlation** (``xc.f90``)

LDA and GGA exchange-correlation functionals (Barth-Hedin, PBE, etc.), via libXC.

**9. Utilities** (``math.f90``, ``array.f90``, ``lists.f90``, ``sparse.f90``, ``logger.f90``, ``safe_alloc.f90``, ``report.f90``, etc.)

Mathematical functions, array/sparse operations, logging, memory tracking,
run reporting.

Calculation Workflow
====================

A typical RS-LMTO-ASA calculation follows:

1. **Initialization** (``main.f90``)

   - Initialize MPI and OpenMP
   - Parse optional command-line argument (input filename, default: ``input.nml``)
   - Construct the ``calculation`` object (``restore_to_default`` + ``build_from_file``)

2. **Dispatch** (``calculation.f90::process``)

   Three independent select-cases run in sequence — ``pre_processing``,
   ``processing``, ``post_processing`` — each defaulting to ``'none'``. Each
   ``*_processing_*`` subroutine is a self-contained pipeline that constructs
   its own ``control``/``lattice``/``charge``/``hamiltonian``/... objects from
   the input file.

3. **SCF Loop** (``self.f90::run``, via ``pre_processing_bravais`` / ``_buildsurf`` / ``_newclu*``)

   For each iteration:

   a. Construct Hamiltonian (``hamiltonian_build.f90``)
   b. Compute recursion coefficients (``recursion_*``, dispatched on ``control%recur``)
   c. Calculate Green's functions (``green.f90``)
   d. Integrate DOS to get charge density (``density_of_states.f90``)
   e. Update effective potential (``charge.f90``)
   f. Apply mixing (``mix.f90``)
   g. Check convergence

4. **Properties** (``calculation.f90::post_processing_*``)

   - Band structure / DOS via the k-space route (``reciprocal_*``)
   - Exchange parameters (``exchange.f90``)
   - Conductivity (``conductivity.f90``)
   - Orbital moments (``recursion%chebyshev_orbital_mod``)

5. **Output**

   - Write results to files
   - Print final report (``report.f90``)

Code Organization Principles
=============================

**Modularity**

Each physical quantity or calculation method has a dedicated Fortran module with:

- Derived type definition (e.g., ``type :: lattice``)
- Constructor, ``destructor``, and ``restore_to_default``
- Core procedures (type-bound methods)
- I/O routines (``build_from_file``, ``print_state``)

Larger modules split their *implementation* across Fortran submodules
(``submodule (parent_mod) name``) while keeping the type and its interface in
the parent module — the class is never flattened into bare module procedures.

**Encapsulation**

- Private data members where practical
- Public interface procedures
- Pointer-based communication between modules

**Parallelization**

- OpenMP for shared-memory loops (e.g. atom loops, Chebyshev sweeps)
- MPI for inter-node communication
- Hybrid OpenMP/MPI support; optional CUDA plugin for the recursion kernels

**Precision Control**

All floating-point calculations use ``precision_mod`` to ensure consistent precision:

.. code-block:: fortran

   use precision_mod, only: rp  ! rp = kind(1.0d0) for double precision

Key Files for Common Tasks
===========================

**Add a new input parameter:**

- ``control.f90`` (add field to ``type control``)
- ``source/include_codes/namelist_generator`` (namelist generation)
- ``docs/source/keywords/*.rst`` (document the parameter)

**Modify Hamiltonian construction:**

- ``hamiltonian_build.f90`` (main build routines)
- ``potential.f90`` (TB parameters)
- ``hamiltonian_ccor.f90`` / ``hamiltonian_hubbard.f90`` (corrections)
- ``recursion_*`` (if recursion coefficients affected)

**Add a new post-processing property:**

- ``calculation.f90`` (add the string to ``check_post_processing`` and a case
  to the ``post_processing`` select-case; reuse
  ``prepare_post_processing_stack`` for the common setup)
- Relevant module (e.g. ``exchange.f90``, ``conductivity.f90``)

**Optimize critical loops:**

- ``chebyshev_fast.f90`` / ``haydock_fast.f90`` (fast recursion kernels)
- ``recursion_*`` submodules (dispatch)
- Marked with ``!$omp parallel`` directives

Provenance
==========

This documentation is derived from:

- File structure: directory walk of ``source/``
- Main entry point: ``source/main.f90``, ``source/os.f90`` (argument parser)
- Build system: ``CMakeLists.txt``, ``source/CMakeLists.txt``, ``cmake/SetFortranFlags.cmake``
- Module / submodule structure: Fortran ``module`` / ``submodule`` declarations in each source file
- Calculation workflow: ``calculation.f90``, ``self.f90``, ``docs/DEVELOPER_MAP.md``

See Also
========

- ``docs/DEVELOPER_MAP.md`` - Task-oriented developer cheat-sheet (call chains, kernels, testing map)
- :doc:`getting_started` - Installation and running
- :ref:`theory/lmto_asa_overview` - Theoretical foundations
- :ref:`reference/module_overview` - Module and type documentation
