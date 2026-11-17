.. _theory/green_functions:

=======================================
Green's Functions
=======================================

Introduction
============

Green's functions are central to the RS-LMTO-ASA approach. Instead of diagonalizing the 
Hamiltonian, the code computes the **single-particle Green's function** at various energies,
from which all spectral and thermodynamic properties can be derived.

Green's Function Definition
============================

The **retarded single-particle Green's function** is defined as:

.. math::

   G_{ij}(E) = \langle i | (E - H + i\eta)^{-1} | j \rangle

where:

- :math:`E` is the energy (real)
- :math:`\eta > 0` is an infinitesimal broadening (convergence parameter)
- :math:`H` is the effective one-electron Hamiltonian
- :math:`i, j` index spatial positions or orbital states

**Compact Matrix Form:**

.. math::

   \mathbf{G}(E) = [E \mathbf{S} - \mathbf{H} + i\eta]^{-1}

where :math:`\mathbf{S}` is the overlap matrix.

Physical Interpretation
=======================

**Spectral Function and Density of States:**

The imaginary part of the Green's function gives the density of states (DOS):

.. math::

   \rho(E) = -\frac{1}{\pi} \text{Im} \sum_i G_{ii}(E)

**Local density of states (LDOS):**

.. math::

   \rho_i(E) = -\frac{1}{\pi} \text{Im} G_{ii}(E)

Shows the contribution to DOS from orbital :math:`i`.

**Charge Density:**

Integrate LDOS below the Fermi energy:

.. math::

   n_i = \int_{-\infty}^{E_F} dE \, \rho_i(E) \, f(E)

where :math:`f(E)` is the Fermi-Dirac distribution.

**Response Functions:**

Magnetic susceptibilities, linear response, etc. are derived from Green's function derivatives.

Types of Green's Functions in RS-LMTO
=====================================

**On-site Green's Function:**

.. math::

   G_{ii}(E) = \langle i | (E - H + i\eta)^{-1} | i \rangle

Contains eigenvalue spectrum of local site.

**Inter-site Green's Function:**

.. math::

   G_{ij}(E) = \langle i | (E - H + i\eta)^{-1} | j \rangle \quad (i \neq j)

Describes hopping between sites :math:`i` and :math:`j`.

**Block Green's Function:**

For multi-orbital systems, the full matrix :math:`G` with dimension (# atoms) :math:`\times` (# orbitals per atom).

Real-Space Calculation
======================

**Advantage of Real-Space Approach:**

Instead of computing eigenvalues and eigenvectors, RS-LMTO uses the **recursion method**
to compute :math:`G_{ij}(E)` at desired energies directly, without diagonalization.

**Key Benefits:**

- No need for k-point mesh sampling
- Direct access to real-space quantities
- Efficient for non-periodic systems (surfaces, clusters, defects)
- Linear scaling (approximately) with cluster size

Recursive Construction
~~~~~~~~~~~~~~~~~~~~~~~

The Green's function can be recursively built using a continued-fraction representation:

.. math::

   G_{00}(E) = \frac{1}{E - a_0 - \frac{b_1^2}{E - a_1 - \frac{b_2^2}{E - a_2 - \cdots}}}

where :math:`a_n` and :math:`b_n` are **recursion coefficients** computed via the Lanczos algorithm.

This avoids explicit matrix inversion; see :ref:`theory/recursion_method` for details.

Density of States from Green's Function
========================================

**Total DOS:**

.. math::

   \rho(E) = -\frac{1}{\pi} \text{Im} \text{Tr}[G(E)]

Computed by integrating Chebyshev polynomial expansions (more stable than direct Im).

**Chebyshev Expansion:**

Green's function spectral function expanded in Chebyshev polynomials:

.. math::

   \rho(E) = \frac{1}{\pi \sqrt{a^2 - (E-b)^2}} \left[ g_0 \mu_0 + 2 \sum_{n \ge 1} g_n \mu_n T_n\left(\frac{E-b}{a}\right) \right]

where :math:`\mu_n` are Chebyshev moments (computed recursively), :math:`g_n`
is the Jackson kernel, :math:`a = (E_{\max} - E_{\min})/1.7` and
:math:`b = (E_{\max} + E_{\min})/2` (see :ref:`theory/recursion_method`).

**Advantages:**

- Numerically stable
- Efficient convergence
- No need to evaluate at many energy points

Implementation in RS-LMTO
=========================

**Key Routines in Code:**

- ``green.f90::sgreen()`` - Scalar (on-site) Green's function
- ``green.f90::bgreen()`` - Block Green's function for all sites
- ``green.f90::block_green_eta()`` - Deprecated wrapper of ``block_green_core``; on-site G at the Fermi point plus an imaginary shift
- ``green.f90::chebyshev_green_eta()``, ``chebyshev_green_ij()`` - Chebyshev on-site (imaginary shift) and intersite G
- ``green.f90::chebyshev_dos_dispatch()`` - Chebyshev on-site reconstruction for the DOS
- ``recursion.f90::recur()`` - Lanczos recursion coefficients
- ``density_of_states.f90::chebyshev_dos()`` - DOS from Chebyshev moments

**Data Structures in Code:**

.. code-block:: fortran

   ! From green.f90
   type :: green
       complex(rp), allocatable :: g0(:,:,:,:)      ! On-site Green's function
       complex(rp), allocatable :: gij(:,:,:,:)     ! Inter-site Green's function
       complex(rp), allocatable :: gij_eta(:,:,:,:) ! Intersite G vs imaginary eta
       ! ... more fields for different interactions
   contains
       procedure :: sgreen          ! Compute on-site
       procedure :: bgreen          ! Compute full block
       procedure :: chebyshev_green_ij ! Chebyshev intersite G
   end type green

Imaginary Energy Shift
======================

The ``*_eta`` routines evaluate :math:`G` at the Fermi point with a complex
shift :math:`\eta = i\,r` passed in by the caller; ``bands.f90`` scans
:math:`r` over a quadrature grid. There is no input key for :math:`\eta` and
no energy-dependent :math:`\eta(E)`.

Fermi Level Determination
=========================

The **Fermi level** :math:`E_F` is found by integrating the DOS:

.. math::

   N_e = \int_{-\infty}^{E_F} dE \, \rho(E)

where :math:`N_e` is the target number of electrons.

**Options:**

1. **Fixed Fermi level:** Specify ``fix_fermi = .true.`` in input
2. **Self-consistent:** Iteratively adjust :math:`E_F` to match electron count

The real-space Fermi level is found at zero temperature.

Magnetic Response and Spin-Resolved Functions
==============================================

For magnetic systems:

- Compute **separate Green's functions** for spin-up and spin-down
- Extract **spin magnetization**: :math:`m_z = n_{\uparrow} - n_{\downarrow}`
- Calculate **magnetic moments** per atom from integrated DOS

For **non-collinear** magnetism:

- Full matrix Green's function in spin space
- Spin density as 3D vector quantity

Computing Transport Properties
==============================

The conductivity uses the Kubo-Bastin formula, evaluated from Chebyshev
double moments of the velocity operators (``source/recursion_transport.f90``,
``source/conductivity.f90``).

Advantages and Limitations
==========================

**Advantages:**

[OK] No diagonalization needed
[OK] Direct access to local properties
[OK] Efficient for large clusters
[OK] Natural way to include broadening
[OK] Enables disorder calculations

**Limitations:**

[X] More complex than eigenvalue approach for simple systems
[X] Requires careful energy mesh choices
[X] Convergence depends on broadening parameter :math:`\eta`

Provenance
==========

This section documents:

- **Main implementation:**
  
  - ``source/green.f90`` - All Green's function types and routines
  - ``source/recursion.f90`` - Recursion coefficients for Green's function
  - ``source/density_of_states.f90`` - DOS integration and Chebyshev moments

- **Usage context:**
  
  - ``source/self.f90::run()`` - Calls Green's function in SCF loop
  - ``source/bands.f90`` - Fermi level, DOS and moments from Green's functions
  - ``source/conductivity.f90`` - Transport from Green's functions

See Also
========

- :ref:`theory/recursion_method` - Details on recursive computation
- :ref:`theory/scf_cycle` - Integration into SCF loop
- ``density_of_states.f90`` - DOS calculation module
