.. _reference/module_overview:

=======================================
Module Overview
=======================================

Core Modules
============

**calculation.f90**

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Type**
     - calculation
   * - **Purpose**
     - Main calculation driver

Main procedures:

- ``constructor(input_file)`` - Initialize from namelist
- ``process()`` - Execute full calculation (pre-, main, post-processing)
- ``build_from_file(fname)`` - Parse input

Key members:

- ``pre_processing`` - Setup type ('none', 'bravais', 'buildsurf', etc.)
- ``processing`` - Main calculation ('none', 'sd', etc.)
- ``post_processing`` - Output type ('none', 'dos', 'bands', etc.)

**control.f90**

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Type**
     - control
   * - **Purpose**
     - Control parameters

Major parameters:

- ``nsp`` - Relativistic type (1=SR, 2=FR, 3=NC-SR, 4=NC-FR)
- ``llsp, lld`` - Recursion cutoffs
- ``max_iterations`` - SCF limit
- ``dq_tol`` - Convergence tolerance
- ``idos`` - LDOS output flag
- ``verbose`` - Verbosity flag

See :ref:`keywords/control_parameters` for complete listing.

**lattice.f90**

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Type**
     - lattice
   * - **Purpose**
     - Crystal structure

Major procedures:

- ``build_from_file(fname)`` - Load structure
- ``setup_cluster()`` - Generate cluster geometry
- ``print_state()`` - Display structure info

Key members:

- ``alat`` - Lattice constant
- ``nbulk`` - Atoms per unit cell
- ``nx, ny, nz`` - Cluster dimensions
- ``r2`` - Hopping cutoff radius²

**hamiltonian.f90**

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Type**
     - hamiltonian
   * - **Purpose**
     - Hamiltonian matrix construction

Major procedures:

- ``build_bulkham()`` - Construct bulk Hamiltonian
- ``build_locham()`` - Local Hamiltonian
- ``build_lsham()`` - Spin-orbit coupling terms
- ``build_from_paoflow()`` - Interface with PAOflow

Key members:

- ``ee, eeo`` - Bulk Hamiltonian blocks
- ``hall, hallo`` - Local Hamiltonian blocks
- ``lsham`` - Spin-orbit Hamiltonian
- ``hmag, hhmag`` - Magnetic Hamiltonian blocks

**reciprocal.f90**

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Type**
     - reciprocal
   * - **Purpose**
     - Reciprocal-space k-mesh/k-path workflows (bands, DOS, k-space diagnostics)

Major procedures:

- ``generate_mp_mesh()`` - Build Monkhorst-Pack mesh
- ``build_kspace_hamiltonian()`` - Fourier transform real-space Hamiltonian blocks
- ``diagonalize_hamiltonian()`` - Solve k-space eigenproblems
- ``calculate_band_structure()`` / ``calculate_density_of_states()`` - High-level post-processing drivers

Mode behavior:

- ``kspace_ham_order='second'``: second-order ASA ``H(k)`` with optional Bloch-summed CCOR.
- ``kspace_ham_order='first'``: first-order ``h(k) + H_LS`` with optional Bloch-summed CCOR.
- ``kspace_ham_order='proper'``: deprecated alias for ``second``.

**self.f90**

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Type**
     - self
   * - **Purpose**
     - Self-consistent field (SCF) loop

Major procedures:

- ``process()`` - Main SCF loop
- ``iterate()`` - Single iteration
- ``converge_scf()`` - Check convergence

Key members:

- Pointers to lattice, charge, control, hamiltonian, recursion, green, dos, bands, exchange

**recursion.f90**

.. table::
   :align: left

   +--------+-------------------+
   | Type   | recursion         |
   +--------+-------------------+
   | Purpose| Recursion method  |
   +--------+-------------------+

Major procedures:

- ``recur()`` - Full Lanczos recursion
- ``hop()`` - Single Lanczos step
- ``crecal()`` - Recursion coefficient calculation
- ``chebyshev_recur_ll()`` - Chebyshev moments
- ``block_green()`` - Block recursion for Green's function

Key members:

- ``a, b2`` - Recursion coefficients (scalar)
- ``a_b, b2_b`` - Recursion coefficients (block)
- ``psi, pmn`` - Lanczos wavefunctions
- ``mu_n`` - Chebyshev moments

**green.f90**

.. table::
   :align: left

   +--------+-------------------+
   | Type   | green             |
   +--------+-------------------+
   | Purpose| Green's functions |
   +--------+-------------------+

Major procedures:

- ``sgreen()`` - On-site Green's function
- ``bgreen()`` - Block Green's function
- ``block_green()`` - Full block computation
- ``chebyshev_green()`` - Chebyshev expansion
- ``calculate_intersite_gf()`` - Inter-site Green's functions

Key members:

- ``g0`` - On-site Green's function
- ``gij, gji`` - Inter-site Green's functions
- ``gij_eta, gji_eta`` - Energy-dependent broadening version

Property Modules
================

**density_of_states.f90**

.. table::
   :align: left

   +--------+-------------------+
   | Type   | dos               |
   +--------+-------------------+
   | Purpose| Density of states |
   +--------+-------------------+

Procedures:

- ``density()`` - Integrate Green's function for DOS
- ``chebyshev_dos()`` - DOS from Chebyshev moments
- ``chebyshev_dos_full()`` - Per-orbital DOS

Members:

- ``doscheb`` - DOS array
- Pointers to recursion, symbolic_atom, lattice, control, energy

**bands.f90**

.. table::
   :align: left

   +--------+-------------------+
   | Type   | bands             |
   +--------+-------------------+
   | Purpose| Band structure    |
   +--------+-------------------+

Procedures:

- ``density()`` - Calculate band structure
- ``write_bands()`` - Output bands

**exchange.f90**

.. table::
   :align: left

   +--------+-------------------+
   | Type   | exchange          |
   +--------+-------------------+
   | Purpose| Exchange coupling |
   +--------+-------------------+

Procedures:

- ``calculate()`` - Compute exchange interactions
- ``write_exchange()`` - Output J values

Key members:

- J matrix elements between atoms

**conductivity.f90**

.. table::
   :align: left

   +--------+-------------------+
   | Type   | conductivity      |
   +--------+-------------------+
   | Purpose| Transport         |
   +--------+-------------------+

Procedures:

- ``calculate()`` - Compute conductivity tensor
- ``write_conductivity()`` - Output σ values

Utility Modules
===============

**precision.f90**

- Defines ``rp`` (real precision)
- Standard: double precision (kind=8)

**math.f90**

- Constants: π, √π, etc.
- Functions: distance, angle, cross_product, normalization, rotation

**string.f90**

- Utilities: int2str, real2str, fmt, path_join, freplace

**logger.f90**

- Logging system (info, warning, error, debug levels)
- Output control

**timer.f90**

- Performance timing
- ``g_timer`` global timer instance

**mpi.f90**

- MPI wrappers for Fortran
- Handles rank, size, barriers, reductions

**charge.f90**

- Charge density management
- Potential updates
- Madelung corrections

**mix.f90**

- Density mixing strategies
- Linear and Broyden implementations

**element.f90**

- Element properties
- Database interface

**potential.f90**

- TB potential parameters
- Database loading

**energy.f90**

- Energy mesh management
- Fermi level calculation

**symbolic_atom.f90**

- Symbolic atom representation
- Atomic site information

Data Flow Overview
==================

.. code-block:: text

   Input (namelist)
      ↓
   Control + Lattice + Element + Potential
      ↓
   Hamiltonian Construction
      ↓
   SCF Loop:
      ├-> Recursion Coefficients
      ├-> Green's Functions
      ├-> DOS/Charge Density
      ├-> Potential Update
      └-> Convergence Check
      ↓
   Properties (Exchange, DOS, Bands, Conductivity)
      ↓
   Output Files

Provenance
==========

All major modules in ``source/`` directory:

- Each module in its own ``.f90`` file
- Module name matches file name (lowercase)
- Type definition and procedures together

Typical module structure:

.. code-block:: fortran

   module my_module_name
      use other_modules
      implicit none
      private
      
      type, public :: my_type
         ! Members
      contains
         procedure :: method1
         procedure :: method2
      end type my_type
      
      interface my_type
         procedure :: constructor
      end interface
      
   contains
      ! Implementations
   end module my_module_name

See Also
========

- :doc:`../code_structure` - Directory layout
- :doc:`data_structures` - Type definitions
- :doc:`algorithms` - Key algorithms
