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

- $E$ is the energy (real)
- $\eta > 0$ is an infinitesimal broadening (convergence parameter)
- $H$ is the effective one-electron Hamiltonian
- $i, j$ index spatial positions or orbital states

**Compact Matrix Form:**

.. math::

   \mathbf{G}(E) = [E \mathbf{S} - \mathbf{H} + i\eta]^{-1}

where $\mathbf{S}$ is the overlap matrix.

Physical Interpretation
=======================

**Spectral Function and Density of States:**

The imaginary part of the Green's function gives the density of states (DOS):

.. math::

   \rho(E) = -\frac{1}{\pi} \text{Im} \sum_i G_{ii}(E)

**Local density of states (LDOS):**

.. math::

   \rho_i(E) = -\frac{1}{\pi} \text{Im} G_{ii}(E)

Shows the contribution to DOS from orbital $i$.

**Charge Density:**

Integrate LDOS below the Fermi energy:

.. math::

   n_i = \int_{-\infty}^{E_F} dE \, \rho_i(E) \, f(E)

where $f(E)$ is the Fermi-Dirac distribution.

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

Describes hopping between sites $i$ and $j$.

**Block Green's Function:**

For multi-orbital systems, the full matrix $G$ with dimension (# atoms) $\times$ (# orbitals per atom).

Real-Space Calculation
======================

**Advantage of Real-Space Approach:**

Instead of computing eigenvalues and eigenvectors, RS-LMTO uses the **recursion method**
to compute $G_{ij}(E)$ at desired energies directly, without diagonalization.

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

where $a_n$ and $b_n$ are **recursion coefficients** computed via the Lanczos algorithm.

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

   \rho(E) = \frac{1}{\sqrt{(E_{\max} - E)(E - E_{\min})}} \sum_n \mu_n T_n\left(\frac{2E - E_{\max} - E_{\min}}{E_{\max} - E_{\min}}\right)

where $\mu_n$ are Chebyshev moments (computed recursively).

**Advantages:**

- Numerically stable
- Efficient convergence
- No need to evaluate at many energy points

Implementation in RS-LMTO
=========================

**Key Routines in Code:**

- ``green.f90::sgreen()`` - Scalar (on-site) Green's function
- ``green.f90::bgreen()`` - Block Green's function for all sites
- ``green.f90::block_green_eta()`` - Green's function with energy-dependent broadening
- ``green.f90::chebyshev_green()`` - Chebyshev moment representation
- ``recursion.f90::recur()`` - Lanczos recursion coefficients
- ``density_of_states.f90::chebyshev_dos()`` - DOS from Chebyshev moments

**Data Structures in Code:**

.. code-block:: fortran

   ! From green.f90
   type :: green
       complex(rp), allocatable :: g0(:,:,:,:)      ! On-site Green's function
       complex(rp), allocatable :: gij(:,:,:,:)     ! Inter-site Green's function
       complex(rp), allocatable :: gij_eta(:,:,:,:) ! G with energy-dependent eta
       ! ... more fields for different interactions
   contains
       procedure :: sgreen          ! Compute on-site
       procedure :: bgreen          ! Compute full block
       procedure :: chebyshev_green ! Chebyshev expansion
   end type green

Energy Dependent Broadening
===========================

The code supports **energy-dependent broadening** via parameter $\eta(E)$:

.. math::

   G_{ij}(E) = \langle i | (E - H + i\eta(E))^{-1} | j \rangle

**Typical forms:**

1. **Constant:** :math:`\eta = \eta_0` (broadening parameter)
2. **Energy-dependent:** :math:`\eta(E) = \eta_0 (1 + |E - E_F|)` (sharper near Fermi level)

Controlled by keyword ``ieta`` in input namelists (see :ref:`keywords/energy_mesh`).

Fermi Level Determination
=========================

The **Fermi level** $E_F$ is found by integrating the DOS:

.. math::

   N_e = \int_{-\infty}^{E_F} dE \, \rho(E) \, f(E)

where $N_e$ is the target number of electrons.

**Options:**

1. **Fixed Fermi level:** Specify ``fix_fermi = .true.`` in input
2. **Self-consistent:** Iteratively adjust $E_F$ to match electron count
3. **Temperature effects:** Include Fermi-Dirac distribution with temperature $T$

Magnetic Response and Spin-Resolved Functions
==============================================

For magnetic systems:

- Compute **separate Green's functions** for spin-up and spin-down
- Extract **spin magnetization**: $m_z = n_{\uparrow} - n_{\downarrow}$
- Calculate **magnetic moments** per atom from integrated DOS

For **non-collinear** magnetism:

- Full matrix Green's function in spin space
- Spin density as 3D vector quantity

Computing Transport Properties
==============================

**Conductivity tensor:**

.. math::

   \sigma_{\alpha\beta}(E) = \frac{e^2}{2\pi} \text{Im} \text{Tr}[v_\alpha G(E) v_\beta G(E)]

where $v_\alpha$ are velocity operators, derived from Hamiltonian.

See ``conductivity.f90`` for implementation.

**Kubo-Greenwood Formula:**

Relates conductivity to Green's function:

.. math::

   \sigma(E) = \frac{e^2}{V} \text{Im} \sum_{i,j} G_{ij}(E) v_{ij} G_{ji}(E) v_{ji}

Advantages and Limitations
==========================

**Advantages:**

✓ No diagonalization needed
✓ Direct access to local properties
✓ Efficient for large clusters
✓ Natural way to include broadening
✓ Enables disorder calculations

**Limitations:**

✗ More complex than eigenvalue approach for simple systems
✗ Requires careful energy mesh choices
✗ Convergence depends on broadening parameter $\eta$

Provenance
==========

This section documents:

- **Main implementation:**
  
  - ``source/green.f90`` - All Green's function types and routines
  - ``source/recursion.f90`` - Recursion coefficients for Green's function
  - ``source/density_of_states.f90`` - DOS integration and Chebyshev moments

- **Usage context:**
  
  - ``source/self.f90::scf_iteration()`` - Calls Green's function in SCF loop
  - ``source/bands.f90`` - Band structure from Green's functions
  - ``source/conductivity.f90`` - Transport from Green's functions

See Also
========

- :ref:`theory/recursion_method` - Details on recursive computation
- :ref:`theory/scf_cycle` - Integration into SCF loop
- ``density_of_states.f90`` - DOS calculation module
