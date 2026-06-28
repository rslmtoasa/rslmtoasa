.. _reference_index:

=======================================
Reference
=======================================

Code documentation and module overview.

.. toctree::
   :maxdepth: 2

   module_overview
   data_structures
   algorithms
   report_module
   exchange_module
   conductivity_module
   bands_module

Overview
========

This section provides reference documentation for the RS-LMTO-ASA codebase:

- **Module Overview** - Main Fortran modules and their purposes
- **Data Structures** - Type definitions and data organization
- **Algorithms** - Key computational algorithms used

While RS-LMTO-ASA is primarily a Fortran codebase, this section provides a curated
overview rather than exhaustive API documentation (which is better suited for 
Doxygen or similar tools).

For detailed code documentation, see:

- Source code comments (inline documentation)
- Module descriptions in .f90 files
- Generated API docs (if available)

Module Organization
===================

**Core Calculation Modules:**

- **calculation.f90** - Main calculation driver
- **control.f90** - Control parameters
- **lattice.f90** - Crystal structure and geometry
- **hamiltonian.f90** - Hamiltonian matrix construction
- **self.f90** - Self-consistent field (SCF) loop
- **recursion.f90** - Lanczos/Chebyshev recursion
- **green.f90** - Green's function calculations

**Property Modules:**

- **density_of_states.f90** - DOS calculation and integration
- **bands.f90** - Band structure post-processing
- **exchange.f90** - Exchange interactions
- **conductivity.f90** - Transport properties (conductivity tensor)

**Utility Modules:**

- precision.f90 - Floating-point precision
- math.f90 - Mathematical constants and functions
- string.f90 - String utilities
- logger.f90 - Logging and output
- timer.f90 - Performance timing
- mpi.f90 - MPI wrappers
- And others...

Key Procedures by Topic
========================

**Initialization:**

- ``calculation::build_from_file()`` - Read input and initialize
- ``lattice::build_from_file()`` - Load structure
- ``control::build_from_file()`` - Read control parameters

**SCF Loop:**

- ``self::process()`` - Main SCF iteration loop
- ``self::iterate()`` - Single SCF iteration
- ``self::converge_scf()`` - Check convergence

**Hamiltonian Construction:**

- ``hamiltonian::build_bulkham()`` - Bulk Hamiltonian
- ``hamiltonian::build_locham()`` - Local Hamiltonian
- ``hamiltonian::build_lsham()`` - Spin-orbit coupling

**Electronic Structure:**

- ``recursion::recur()`` - Lanczos recursion
- ``recursion::hop()`` - Single recursion step
- ``recursion::chebyshev_recur_ll()`` - Chebyshev moments
- ``green::sgreen()`` - On-site Green's function
- ``green::bgreen()`` - Block Green's function

**Properties:**

- ``density_of_states::density()`` - Integrate DOS
- ``bands::calculate()`` - Band structure
- ``exchange::calculate()`` - Exchange interactions
- ``conductivity::calculate()`` - Conductivity tensor

Quick Code Navigation
=====================

**To understand SCF:**

Start with ``source/self.f90::process()`` -> calls ``calculation.f90``, ``hamiltonian.f90``, 
``recursion.f90``, ``green.f90``, ``density_of_states.f90``, ``mix.f90``.

**To understand recursion:**

Start with ``source/recursion.f90::recur()`` -> calls ``hop()``, ``crecal()``.

**To understand Green's functions:**

Start with ``source/green.f90::sgreen()`` -> uses recursion coefficients from recursion module.

**To understand magnetic properties:**

Look for ``nsp > 1`` handling in:

- ``control.f90`` - Collinearity/relativistic settings
- ``hamiltonian.f90`` - Separate Hamiltonians for spins
- ``green.f90`` - Spin-up/down Green's functions
- ``exchange.f90`` - Inter-site coupling

See Also
========

- :doc:`../code_structure` - Code directory layout
- :doc:`../theory/index` - Theoretical background
- Individual keyword pages for algorithm details
