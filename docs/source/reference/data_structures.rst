.. _reference/data_structures:

=======================================
Data Structures
=======================================

Key Fortran Derived Types
=========================

This section documents the major data structures (Fortran types) in RS-LMTO-ASA.

calculation
-----------

**File:** source/calculation.f90

**Purpose:** Main calculation workflow

**Key members:**

.. code-block:: fortran

   type :: calculation
      character(len=sl) :: pre_processing    ! 'none', 'bravais', 'buildsurf', etc.
      character(len=sl) :: processing        ! 'none', 'sd', etc.
      character(len=sl) :: post_processing   ! 'none', 'dos', 'bands', 'exchange', etc.
      logical :: verbose                     ! Enable verbose output
      character(len=sl) :: fname             ! Input namelist filename
   end type calculation

control
-------

**File:** source/control.f90

**Purpose:** Calculation control parameters

**Key members:**

.. code-block:: fortran

   type :: control
      integer :: nsp                         ! 1=SR, 2=FR, 3=NC-SR, 4=NC-FR
      integer :: llsp, lld                   ! Recursion cutoffs
      integer :: max_iterations              ! SCF iteration limit
      integer :: idos                        ! LDOS output type
      integer :: nlim, npold                 ! Cluster/history settings
      integer :: mext                        ! Spin acceleration
      logical :: lrot, incorb                ! Rotation settings
      real(rp) :: dq_tol                     ! Convergence tolerance
      logical :: verbose                     ! Verbose output
   end type control

lattice
-------

**File:** source/lattice.f90

**Purpose:** Crystal structure and geometry

**Key members:**

.. code-block:: fortran

   type :: lattice
      real(rp) :: alat                       ! Lattice constant (Å)
      real(rp) :: r2                         ! Hopping cutoff (Å²)
      integer :: nbulk                       ! Atoms per unit cell
      integer :: nx, ny, nz                  ! Cluster dimensions
      ! ... atomic positions, lattice vectors, etc.
   end type lattice

hamiltonian
-----------

**File:** source/hamiltonian.f90

**Purpose:** Electronic Hamiltonian matrix

**Key members:**

.. code-block:: fortran

   type :: hamiltonian
      complex(rp), allocatable :: ee(:,:,:,:)     ! Bulk Hamiltonian
      complex(rp), allocatable :: hall(:,:,:,:)   ! Local Hamiltonian
      complex(rp), allocatable :: lsham(:,:,:)    ! Spin-orbit coupling
      complex(rp), allocatable :: obarm(:,:,:)    ! Overlap matrix
      logical :: hoh, local_axis, orb_pol         ! Options
   end type hamiltonian

recursion
---------

**File:** source/recursion.f90

**Purpose:** Recursion method coefficients and workspace

**Key members:**

.. code-block:: fortran

   type :: recursion
      real(rp), allocatable :: a(:,:,:,:)         ! Lanczos coeff a_n
      real(rp), allocatable :: b2(:,:,:,:)        ! Lanczos coeff b_n²
      complex(rp), allocatable :: a_b(:,:,:,:)    ! Block form
      complex(rp), allocatable :: b2_b(:,:,:,:)   ! Block form
      complex(rp), allocatable :: psi(:,:,:)      ! Lanczos vectors
      complex(rp), allocatable :: mu_n(:,:,:,:)   ! Chebyshev moments
   end type recursion

green
-----

**File:** source/green.f90

**Purpose:** Single-particle Green's functions

**Key members:**

.. code-block:: fortran

   type :: green
      complex(rp), allocatable :: g0(:,:,:,:)      ! On-site GF
      complex(rp), allocatable :: gij(:,:,:,:)     ! Inter-site GF
      complex(rp), allocatable :: gij_eta(:,:,:,:) ! G with broadening
      ! ... more components for different couplings
   end type green

dos (density_of_states)
-----------------------

**File:** source/density_of_states.f90

**Purpose:** Density of states calculation and storage

**Key members:**

.. code-block:: fortran

   type :: dos
      real(rp), allocatable :: doscheb(:,:,:)     ! DOS array
      ! Pointers to: recursion, symbolic_atom, lattice, control, energy
   end type dos

charge
------

**File:** source/charge.f90

**Purpose:** Charge density and electrostatic potential

**Key members:**

.. code-block:: fortran

   type :: charge
      real(rp), allocatable :: dq(:)              ! Charge neutrality corrections
      real(rp) :: cht, trq                        ! Charge transfer, torque
      ! ... monopole, dipole, higher multipole moments
      ! ... Madelung/Coulomb potential components
   end type charge

self
----

**File:** source/self.f90

**Purpose:** Self-consistent field iteration

**Key members:**

.. code-block:: fortran

   type :: self
      ! Pointers to other modules:
      class(lattice), pointer :: lattice
      class(charge), pointer :: charge
      class(control), pointer :: control
      class(hamiltonian), pointer :: hamiltonian
      class(recursion), pointer :: recursion
      class(green), pointer :: green
      class(dos), pointer :: dos
      class(mix), pointer :: mix
      ! ... others
   end type self

mix
---

**File:** source/mix.f90

**Purpose:** Density mixing for SCF

**Key members:**

.. code-block:: fortran

   type :: mix
      character(len=sl) :: mixing_type            ! 'linear' or 'broyden'
      real(rp) :: alpha                           ! Mixing parameter
      integer :: n_history                        ! Broyden history length
      ! ... arrays for old densities and residuals
   end type mix

element
-------

**File:** source/element.f90

**Purpose:** Atomic element properties

**Key members:**

.. code-block:: fortran

   type :: element
      character(len=10) :: symbol                 ! Element name
      integer :: atomic_number, f_core, core
      real(rp) :: valence
      integer :: num_quant_s, num_quant_p, num_quant_d
   end type element

potential
---------

**File:** source/potential.f90

**Purpose:** Tight-binding potential parameters

**Key members:**

.. code-block:: fortran

   type :: potential
      integer :: lmax                             ! Max orbital angular momentum
      real(rp), allocatable :: center_band(:,:)   ! Orbital center energies
      real(rp), allocatable :: width_band(:,:)    ! Bandwidth parameters
      ! ... ORT parameters (c, srdel, qpar, ppar, etc.)
   end type potential

energy
------

**File:** source/energy.f90

**Purpose:** Energy mesh and Fermi level

**Key members:**

.. code-block:: fortran

   type :: energy
      integer :: channels_ldos                    ! Number of energy points
      real(rp) :: fermi, energy_min, energy_max   ! Energy range
      real(rp) :: edel                            ! Energy mesh spacing
      real(rp), allocatable :: ene(:)             ! Energy array
      logical :: fix_fermi                        ! Fix Fermi level?
   end type energy

symbolic_atom
--------------

**File:** source/symbolic_atom.f90

**Purpose:** Symbolic (representative) atom for cluster

**Key members:**

.. code-block:: fortran

   type :: symbolic_atom
      integer :: atom_index                       ! Position in cluster
      integer :: type                             ! Atom type (element)
      real(rp) :: position(3)                     ! Cartesian position
      integer :: neighbors                        ! Number of neighbors
      ! ... magnetic moment, charge, etc.
   end type symbolic_atom

xc
--

**File:** source/xc.f90

**Purpose:** Exchange-correlation functional

**Key members:**

.. code-block:: fortran

   type :: xc
      character(len=3) :: txch                    ! Functional type
      integer :: txc, nss                         ! Functional ID, spin settings
      real(rp) :: exchf, alpm, aa, bb, cc, ...    ! Functional parameters
   end type xc

bands
-----

**File:** source/bands.f90

**Purpose:** Band structure calculation

**Key members:**

.. code-block:: fortran

   type :: bands
      ! Pointers to: green, lattice, energy, recursion, dos, etc.
      ! Band energy array
   end type bands

exchange
--------

**File:** source/exchange.f90

**Purpose:** Exchange interaction calculations

**Key members:**

.. code-block:: fortran

   type :: exchange
      ! Pointers to: green, lattice, recursion, control, etc.
      ! Exchange coupling matrix J(i,j)
   end type exchange

conductivity
------------

**File:** source/conductivity.f90

**Purpose:** Electron transport properties

**Key members:**

.. code-block:: fortran

   type :: conductivity
      ! Pointers to: control, lattice, hamiltonian, green, etc.
      ! Conductivity tensor σ(α,β,E)
   end type conductivity

Memory Layout
=============

**Array dimensions (typical):**

- ``ee(norb, norb, natom, natom)`` - Inter-atomic Hamiltonian blocks
- ``hall(norb, norb, nclust, 1)`` - Local Hamiltonian
- ``g0(norb, norb, 1, nclust)`` - On-site Green's function
- ``gij(norb, norb, nclust, nclust)`` - Inter-site Green's function
- ``a(nclust, nclust, norb, ll_max)`` - Recursion coefficients
- ``doscheb(nclust, norb, channels_ldos)`` - DOS

where:

- ``norb`` ≈ 3×lmax (typically 9 for lmax=2)
- ``nclust`` ≈ nbulk × nx × ny × nz
- ``ll_max`` = llsp or lld

Complex vs. Real
================

- **Real arrays:** `a`, `b2` (recursion coefficients are real)
- **Complex arrays:** `ee`, `hall`, `gij`, `lsham`, `psi` (quantum amplitudes)
- **Real arrays:** all energies, positions, mesh points
- **Integer arrays:** atom indices, type indices

Precision Control
=================

All floating-point variables use ``rp`` from precision.f90:

.. code-block:: fortran

   use precision_mod, only: rp
   
   real(rp) :: fermi, charge
   complex(rp) :: green_function

Default: ``rp = kind(1.0d0)`` (double precision)

Allocation Patterns
===================

**Typical allocation (in constructor):**

.. code-block:: fortran

   allocate(this%a(natom, natom, norb, ll_max))
   allocate(this%g0(norb, norb, 1, natom))

**Deallocation (in destructor):**

.. code-block:: fortran

   if (allocated(this%a)) deallocate(this%a)
   if (allocated(this%g0)) deallocate(this%g0)

Safe allocation wrapper available: ``safe_alloc_mod`` (optional, controlled by ``USE_SAFE_ALLOC``)

Provenance
==========

All type definitions in source code:

- Type declaration section: first part of each module
- Constructor: function named ``constructor()``
- Destructor: final subroutine (Fortran 2008)
- Methods: type-bound procedures

See Also
========

- :doc:`module_overview` - Module organization
- :doc:`algorithms` - Computational algorithms
- :doc:`../code_structure` - Directory layout
