.. _reference/module_overview:

=======================================
Module Overview
=======================================

.. note::

   The modules ``lattice``, ``hamiltonian``, ``recursion``, and ``reciprocal``
   declare their derived type and procedure interfaces in a parent file
   (e.g. ``hamiltonian.f90``) and carry the implementations in Fortran
   submodules (``hamiltonian_build.f90``, ``hamiltonian_ccor.f90``, ...). The
   type and its public interface are documented here against the parent
   module; see :ref:`code_structure` for the full submodule breakdown.

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

**reciprocal_green.f90** (submodule of ``reciprocal``, B2)

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Type**
     - submodule of ``reciprocal``
   * - **Purpose**
     - k-space Green's-function engine: fills the same ``green%gij``/``gji``,
       ``gij_eta`` ladder, and torque arrays the real-space recursion route
       fills, from two backends, so downstream consumers (exchange, damping,
       conductivity) are route-agnostic.

Major procedures:

- ``fill_green(green_obj, mesh, pairs, sigma_provider)`` - populate the
  Green's-function arrays via ``green_backend='lehmann'`` or ``'dyson'``
  (see :ref:`keywords/output_options` for ``gf_route``/``green_backend``/``green_eta``)

Related files:

- ``sigma_provider.f90`` - abstract self-energy provider interface
  (``sigma_zero`` is the only concrete implementation so far; the
  interface exists so a future CPA/DMFT Σ(z) provider slots in without
  changing ``fill_green``'s call sites)
- ``lehmann_kernel.f90`` - strict Lehmann representation kernel (Σ=0, one
  eigensolve per k, amortized over energies and pairs)
- ``dyson_kernel.f90`` - direct Dyson inversion kernel
  (:math:`[zS-H(k)-\Sigma(z)]^{-1}`, batched per (k, z)); also underlies
  ``reciprocal_bsf.f90`` below
- ``moment_kernel.f90`` / **reciprocal_moments.f90** (submodule of
  ``reciprocal``, B5.1) - exact k-space Chebyshev moment generator
  (``fill_moments``) feeding the existing KPM conductivity pipeline with
  moments computed exactly from eigenpairs instead of recursion

See ``docs/dev/plans/B2_reciprocal_green.md`` for the convention checklist
(phase/bond convention, local spin frames, normalization) and
``docs/dev/plans/B5_route_agnostic_postprocessing.md`` for the moment
generator and the Jij/damping/conductivity triad tests.

**self.f90**

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Type**
     - self
   * - **Purpose**
     - Self-consistent field (SCF) loop

Major procedures:

- ``run()`` - Main SCF loop
- ``is_converged()`` - Check convergence
- ``report()`` - Write iteration/convergence report

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

**reciprocal_bsf.f90** (submodule of ``reciprocal``, B3)

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Type**
     - submodule of ``reciprocal``
   * - **Purpose**
     - Bloch spectral function :math:`A(k,E) = -(1/\pi)\,\text{Im}\,\text{Tr}\,G(k,E+i\eta)`
       along the canonical high-symmetry k-path, using the same
       ``dyson_kspace_inverse`` primitive as ``reciprocal_green.f90``'s
       backend D.

Major procedures:

- ``calculate_bsf(output_file)`` - drives the k-path sweep and writes
  total/spin-up/spin-down traces (``post_processing='bsf'``, see
  :ref:`keywords/output_options`)

Related files:

- ``bsf_kernel.f90`` - the trace kernel (``bsf_spectral_trace``)

**frozen_magnon.f90** (B1)

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Type**
     - frozen_magnon
   * - **Purpose**
     - Namelist type + lifecycle for the frozen-magnon E(q)/:math:`\omega`\ (q)
       sweep (``post_processing='frozen_magnon'``); the sweep driver itself
       lives in ``calculation.f90``
       (``post_processing_frozen_magnon[_acoustic|_auto]``).

See :ref:`keywords/frozen_magnon` for the full parameter reference and
``docs/dev/plans/B1_gbt_frozen_magnons.md`` for the underlying
generalized-Bloch-theorem fix.

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
- B7: ``interfacepot``/``align_regions`` — two-sided deviation-variable
  electrostatics for interface (``calctype='L'``) clusters, reusing
  ``surfmat``'s Madelung matrices rather than a separate interface matrix
  (see ``docs/dev/plans/B7_interfaces_and_vacuum_leads.md``)

**region_registry.f90** (B7)

- Explicit per-site region bookkeeping (``region_registry``,
  ``region_descriptor`` types: kinds vacuum/bulk/active/lead_a/lead_b),
  replacing magic-offset arithmetic previously in ``charge%surfpot``
- Alignment solver: ``alignment_gauge_anchor``, ``alignment_initial_guess``,
  ``alignment_update``, ``alignment_consistency_check`` — fixed-point
  solve for each frozen region's offset V_r, gauge-anchored to the first
  frozen non-vacuum region (gate G-B7-2, ``docs/dev/CONTRACT_FROZEN_REGION.md``)

**vacuum_lead.f90** (B7)

- Generates frozen potential parameters for a semi-infinite vacuum
  (empty-sphere) lattice by running the code's own radial solver at
  constant V(r), replacing hand-set empty-sphere stacks
- Validated against an analytic spherical-Bessel oracle
  (``tests/unit/test_vacuum_lead.f90``)
- Consumed via ``&lattice region_b_kind='vacuum'``
  (see :ref:`keywords/lattice_geometry`)

**electrostatics_multipole.f90** (B6)

- ``compute_dipole_moments`` — l=1 (z) dipole moment Q₁₀ per atom from
  on-site cross-orbital density-matrix elements (s-p_z, p_z-d_z²) and
  exported radial partial-wave amplitudes, with SCF mixing
- Feeds ``charge.f90``'s Madelung potential shift via the pre-existing
  ``dsz`` dipole-monopole matrix (no new lattice-sum machinery — narrower
  in scope than the original blueprint; see
  ``docs/dev/plans/B6_surface_electrostatics.md``)
- Gated by ``&control dipole_electrostatics`` (see :ref:`keywords/control_parameters`)

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
