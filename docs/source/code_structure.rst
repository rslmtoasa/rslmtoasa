.. _code_structure:

=======================================
Code Structure
=======================================

Directory Layout
================

.. code-block:: text

   rslmto_devel/
   ├── source/                    # Fortran source code
   │   ├── main.f90               # Program entry point
   │   ├── calculation.f90        # Main calculation driver
   │   ├── control.f90            # Control parameters and settings
   │   ├── self.f90               # Self-consistent field (SCF) loop
   │   ├── lattice.f90            # Crystal structure and geometry
   │   ├── hamiltonian.f90        # Hamiltonian construction
   │   ├── recursion.f90          # Recursion method for Lanczos/Chebyshev
   │   ├── green.f90              # Green's function calculations
   │   ├── energy.f90             # Energy mesh and Fermi level
   │   ├── charge.f90             # Charge density and electrostatics
   │   ├── exchange.f90           # Exchange interaction calculations
   │   ├── conductivity.f90       # Electronic transport properties
   │   ├── bands.f90              # Band structure calculations
   │   ├── density_of_states.f90  # DOS calculations
   │   ├── xc.f90                 # Exchange-correlation functionals
   │   ├── potential.f90          # Potential parameters
   │   ├── element.f90            # Element properties
   │   ├── symbolic_atom.f90      # Atomic site information
   │   ├── namelist_generator.f90 # Input namelist utilities
   │   ├── mix.f90                # Density mixing strategies
   │   ├── face.F90               # Face-centered cubic (FCC) utilities
   │   ├── precision.f90          # Floating-point precision definitions
   │   ├── math.f90               # Mathematical utilities
   │   ├── string.f90             # String manipulation
   │   ├── logger.f90             # Logging and debugging
   │   ├── timer.f90              # Performance timing
   │   ├── mpi.f90                # MPI wrappers
   │   ├── safe_alloc.f90         # Memory tracking
   │   ├── globals.f90            # Global constants
   │   ├── os.f90                 # Operating system utilities
   │   ├── lists.f90              # List and array utilities
   │   ├── array.f90              # Advanced array operations
   │   └── include_codes/         # Include files and macros
   │
   ├── cmake/                     # CMake configuration modules
   │   ├── SetFortranFlags.cmake  # Compiler-specific flags
   │   ├── SetCompileFlag.cmake   # Compilation settings
   │   └── git-watcher.cmake      # Version tracking
   │
   ├── example/                   # Example calculations
   │   ├── bulk/                  # Bulk materials (Si, bccFe, etc.)
   │   ├── surface/               # Surface calculations
   │   ├── impurity/              # Magnetic impurity calculations
   │   ├── exchange/              # Exchange interaction examples
   │   └── conductivity/          # Transport property examples
   │
   ├── tests/                     # Regression tests
   │   ├── regression/            # Test suite
   │   └── run_regression_tests.sh
   │
   ├── build/                     # Build artifacts (CMake generated)
   │   ├── bin/rslmto.x           # Compiled executable
   │   └── modules/               # Compiled Fortran modules
   │
   ├── docs/                      # Sphinx documentation (this)
   │   ├── source/
   │   │   ├── index.rst
   │   │   ├── theory/            # Theoretical foundations
   │   │   ├── user_guide/        # User documentation
   │   │   ├── keywords/          # Input parameter reference
   │   │   └── reference/         # Code reference
   │   └── conf.py
   │
   ├── CMakeLists.txt             # Root CMake configuration
   ├── README.md                  # Quick start guide
   └── requirements.txt           # Python dependencies (if any)

Build System (CMake)
====================

**Main Files:**

- ``CMakeLists.txt`` - Root configuration; defines project, compiler, and options
- ``source/CMakeLists.txt`` - Source file compilation rules
- ``cmake/SetFortranFlags.cmake`` - Compiler detection and flags (gfortran, ifort, etc.)
- ``cmake/SetCompileFlag.cmake`` - Additional compilation rules

**Key Build Targets:**

.. code-block:: bash

   make              # Compile the code
   make html         # Build HTML documentation (from docs/)
   make install      # Install executable and modules
   make clean        # Remove build artifacts
   make distclean    # Remove all CMake files

Module Dependency Graph
=======================

**Core Dependencies:**

.. code-block:: text

   main.f90
   └── calculation.f90 (main driver)
       ├── control.f90
       ├── lattice.f90
       ├── charge.f90
       ├── self.f90 (SCF loop)
       │   ├── hamiltonian.f90
       │   ├── recursion.f90
       │   ├── green.f90
       │   ├── bands.f90
       │   ├── density_of_states.f90
       │   ├── mix.f90
       │   └── exchange.f90
       ├── exchange.f90
       ├── conductivity.f90
       └── ... (output modules)

**Utility Modules (used everywhere):**

- ``precision.f90`` - Real/complex number precision
- ``math.f90`` - Mathematical constants and functions
- ``string.f90`` - String utilities
- ``logger.f90`` - Logging and output
- ``timer.f90`` - Performance timing
- ``mpi.f90`` - MPI wrappers
- ``globals.f90`` - Global constants and parameters

Major Subsystems
================

**1. Initialization & Control** (``control.f90``, ``namelist_generator.f90``)

Reads and validates input parameters from Fortran namelist files. Manages calculation modes 
(SCF, band structure, DOS, transport, etc.).

**2. Structure & Geometry** (``lattice.f90``, ``symbolic_atom.f90``, ``element.f90``)

Handles crystal structure definition, atomic positions, Bravais lattice, and element-specific 
parameters.

**3. Basis & Hamiltonian** (``potential.f90``, ``hamiltonian.f90``)

Constructs the tight-binding Hamiltonian matrix from TB-LMTO parameters and potential functions. 
Supports spin-orbit coupling, magnetic moments, and relativistic corrections.

**4. Electronic Structure** (``recursion.f90``, ``green.f90``)

Core computational engine:

- **Recursion method**: Implements Lanczos and Chebyshev recursion for efficient Green's function 
  and moment calculations
- **Green's functions**: Computes on-site and inter-site Green's functions at complex energies
- Uses recursion coefficients to avoid explicit diagonalization

**5. Self-Consistent Field** (``self.f90``, ``mix.f90``, ``charge.f90``)

Iterative SCF loop:

- Computes charge density from Green's functions
- Updates effective potentials
- Applies density mixing (linear, Broyden, etc.)
- Checks convergence

**6. Properties Calculation** (``density_of_states.f90``, ``bands.f90``, ``exchange.f90``, ``conductivity.f90``)

Post-processing modules:

- **DOS**: Local and total density of states
- **Bands**: Energy band structures
- **Exchange**: Heisenberg exchange parameters (for magnets)
- **Conductivity**: Transport coefficients

**7. Exchange-Correlation** (``xc.f90``)

LDA and GGA exchange-correlation functionals (Barth-Hedin, PBE, etc.).

**8. Utilities** (``math.f90``, ``arrays.f90``, ``logger.f90``, etc.)

Mathematical functions, array operations, logging, memory management.

Calculation Workflow
====================

A typical RS-LMTO-ASA calculation follows:

1. **Initialization** (``main.f90``)
   
   - Initialize MPI and OpenMP
   - Parse optional command-line argument (input filename, default: input.nml)
   - Read input namelist

2. **Setup** (``calculation.f90::build_from_file``)
   
   - Load crystal structure
   - Initialize element parameters
   - Setup atomic basis and potentials

3. **SCF Loop** (``self.f90::process``)
   
   For each iteration:
   
   a. Construct Hamiltonian (``hamiltonian.f90``)
   b. Compute recursion coefficients (``recursion.f90``)
   c. Calculate Green's functions (``green.f90``)
   d. Integrate DOS to get charge density (``density_of_states.f90``)
   e. Update effective potential (``charge.f90``)
   f. Apply mixing (``mix.f90``)
   g. Check convergence

4. **Properties** (``calculation.f90::post_processing_*``)
   
   - Calculate band structure (``bands.f90``)
   - Compute exchange parameters (``exchange.f90``)
   - Calculate transport (``conductivity.f90``)

5. **Output**
   
   - Write results to files
   - Print final report

Code Organization Principles
=============================

**Modularity**

Each physical quantity or calculation method has a dedicated Fortran module with:

- Derived type definition (e.g., ``type :: lattice``)
- Constructor and destructor
- Core procedures (methods)
- I/O routines (``build_from_file``, ``print_state``)

**Encapsulation**

- Private data members
- Public interface functions
- Pointer-based communication between modules

**Parallelization**

- OpenMP for shared-memory loops (e.g., atom loops)
- MPI for inter-node communication
- Hybrid OpenMP/MPI support

**Precision Control**

All floating-point calculations use ``precision_mod`` to ensure consistent precision:

.. code-block:: fortran

   use precision_mod, only: rp  ! rp = kind(1.0d0) for double precision

Key Files for Common Tasks
===========================

**Add a new input parameter:**

- ``control.f90`` (add field to ``type control``)
- ``source/CMakeLists.txt`` (add to namelist generation)
- ``docs/keywords/*.rst`` (document the parameter)

**Modify Hamiltonian construction:**

- ``hamiltonian.f90`` (main routines)
- ``potential.f90`` (TB parameters)
- ``recursion.f90`` (if recursion coefficients affected)

**Add a new output property:**

- ``calculation.f90`` (add to ``post_processing_*``)
- Relevant module (e.g., ``exchange.f90``)

**Optimize critical loops:**

- ``recursion.f90`` (largest computation)
- Marked with ``!$omp parallel`` directives

Provenance
==========

This documentation is derived from:

- File structure: ``ls -la source/``, directory walk
- Main entry point: ``source/main.f90``
- Build system: ``CMakeLists.txt``, ``cmake/SetFortranFlags.cmake``
- Module listing: Fortran ``module`` declarations in each ``.f90`` file
- Calculation workflow: ``calculation.f90``, ``self.f90``, call graphs

See Also
========

- :doc:`getting_started` - Installation and running
- :ref:`theory/lmto_asa_overview` - Theoretical foundations
- :ref:`reference/module_overview` - Module and type documentation
