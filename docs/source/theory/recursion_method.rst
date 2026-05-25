.. _theory/recursion_method:

=======================================
Recursion Method
=======================================

Introduction
============

The **recursion method** (Lanczos or Chebyshev recursion) is the computational heart of RS-LMTO.
It enables efficient calculation of Green's functions and spectral quantities without explicit
matrix diagonalization, providing a significant computational advantage for large clusters.

Lanczos Recursion
=================

The Lanczos algorithm transforms the Hamiltonian into a **tridiagonal matrix** via orthogonal
transformations using a set of iteratively-generated vectors:

**Algorithm:**

Starting with an initial vector :math:`|0\rangle`:

1. :math:`|1\rangle = |0\rangle / \||0\||`
2. For :math:`n = 1, 2, \ldots, N_{\text{rec}}`:
   
   .. math::

      a_n = \langle n | H | n \rangle
      |n+1\rangle = (H - a_n) |n\rangle - b_n^2 |n-1\rangle
      b_{n+1} = \||n+1\||
      |n+1\rangle \rightarrow |n+1\rangle / b_{n+1}

**Result:** Recursion coefficients :math:`\{a_n, b_n\}` that capture the spectral properties.

**Advantages:**

- Produces tridiagonal Hamiltonian in recursion basis
- Requires only Hamiltonian-vector products (no matrix storage or inversion)
- $O(N)$ memory and $O(N)$ per step (for sparse matrices)
- Converges quickly for physical properties (first few moments are important)

Green's Function from Recursion Coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The on-site Green's function is represented as a **continued fraction**:

.. math::

   G_{00}(E) = \frac{1}{E - a_0 - \frac{b_1^2}{E - a_1 - \frac{b_2^2}{E - a_2 - \cdots}}}

This is evaluated numerically without explicit inversion, hence "recursion."

**Scalar Form (single orbital):**

.. math::

   G(E) = (E - H + i\eta)^{-1}

expanded using coefficients $a_0, a_1, \ldots, a_N$ and $b_1, b_2, \ldots, b_N$.

**Block Form (multi-orbital):**

Each orbital/atom combination gets its own recursion, or blocked algorithms handle multiple orbitals simultaneously.

Implementation: Lanczos vs. Chebyshev
=====================================

**Lanczos Recursion:**

- Produces orthogonal Lanczos vectors
- Tridiagonal matrix structure
- Best for Green's function evaluation
- More numerically sensitive near band edges

**Chebyshev Recursion:**

- Uses Chebyshev polynomial basis
- Requires energy rescaling to $[-1, 1]$ interval
- Better numerical stability
- More suitable for DOS and moment calculations

**RS-LMTO uses both:**

- **Lanczos** for direct Green's function computation (``recursion.f90::hop()``)
- **Chebyshev** for DOS and moment integrations (``recursion.f90::chebyshev_*``)

Chebyshev Moments and DOS
=========================

Instead of computing $\rho(E)$ directly from Green's function, RS-LMTO computes **Chebyshev moments**:

.. math::

   \mu_n = \text{Tr}[H^n] / N

These are then expanded:

.. math::

   \rho(E) = \frac{W(\sqrt{E})}{\pi \sqrt{1 - (E/W)^2}} \sum_{n=0}^{N_{\text{cheb}}} \mu_n T_n(E/W)

where:

- $T_n$ are Chebyshev polynomials
- $W$ is a bandwidth parameter for scaling
- $N_{\text{cheb}}$ is the number of moments (typically 100-500)

**Advantages:**

- Very stable numerically
- Fast convergence (Chebyshev polynomials are optimal)
- No explicit density matrix needed

Implementation in RS-LMTO
=========================

**Main Routines:**

.. code-block:: fortran

   ! From recursion.f90
   type :: recursion
      real(rp), allocatable :: a(:,:,:,:)       ! Lanczos coeff a_n
      real(rp), allocatable :: b2(:,:,:,:)      ! Lanczos coeff b_n^2
      complex(rp), allocatable :: a_b(:,:,:,:)  ! Block form of a
      complex(rp), allocatable :: b2_b(:,:,:,:) ! Block form of b^2
      complex(rp), allocatable :: psi(:,:,:)    ! Lanczos wavefunctions
      complex(rp), allocatable :: mu_n(:,:,:,:) ! Chebyshev moments
   contains
      procedure :: hop         ! Lanczos hop (single recursion step)
      procedure :: hop_b       ! Block Lanczos
      procedure :: recur       ! Full Lanczos recursion
      procedure :: recur_b     ! Block recursion
      procedure :: chebyshev_recur_ll ! Chebyshev recursion
      ! ... more procedures
   end type recursion

**Key Procedures:**

1. **``hop()`` and ``hop_b()``**
   
   Single recursion step: compute $a_n$, $b_n$, and next Lanczos vector.

2. **``recur()`` and ``recur_b()``**
   
   Full recursion loop: iterate to cutoff ``llsp`` or ``lld`` (recursion depth).

3. **``chebyshev_recur_ll()``**
   
   Chebyshev moment calculation (stochastic or exact traces).

4. **``chebyshev_dos()``**
   
   Compute DOS from Chebyshev moments via polynomial expansion.

Recursion Cutoff Parameters
===========================

**Control parameters** (from ``control.f90``):

- **``llsp``** - Recursion cutoff for s and p electrons (typical: 50-100)
- **``lld``** - Recursion cutoff for d electrons (typical: 50-100)

These determine the **depth** of recursion (number of coefficients computed).

**Trade-off:**

- Higher cutoff → better accuracy but more computation
- Lower cutoff → faster but coarser spectral resolution

Typical values balance accuracy and speed; see examples in ``example/`` directory.

Stochastic Trace Evaluation
============================

For large systems, computing exact Chebyshev moments $\mu_n = \text{Tr}[H^n]$ is expensive.

**Stochastic approach:**

Use random vectors $|r\rangle$ to estimate the trace:

.. math::

   \mu_n \approx \frac{1}{N_{\text{random}}} \sum_i \langle r_i | H^n | r_i \rangle

**Advantages:**

- Scales as $O(N_{\text{random}})$ independent of system size
- Unbiased estimate with statistical error

**Control:**

Keyword ``random_vec_num`` (in ``control.f90``) specifies number of random vectors.

Numerical Aspects
=================

**Broadening Parameter (eta):**

Added as $i\eta$ to avoid singularities:

.. math::

   G(E) = (E - H + i\eta)^{-1}

Controls spectral resolution. Typical: $\eta \sim 10^{-3}$ to $10^{-2}$ Ry.

**Energy Rescaling for Chebyshev:**

To use Chebyshev polynomials ($T_n : [-1,1] \to \mathbb{R}$), energies are rescaled:

.. math::

   E_{\text{scaled}} = \frac{2E - (E_{\max} + E_{\min})}{E_{\max} - E_{\min}}

where $E_{\min}$, $E_{\max}$ define the spectral window.

**Convergence Criteria:**

- Check when Chebyshev moments drop below numerical noise
- Typical convergence: $\mu_n < 10^{-12}$ for machine precision

Hopping Region and Cluster Truncation
=====================================

In real-space formulation, the recursion operates within a **hopping region** defined by:

1. **Central atom(s)** where properties are computed
2. **Neighboring atoms** at distances up to ``r2`` (input parameter)

Beyond the hopping region, contributions are neglected (Hamiltonian is zero).

**Determines:**

- Accuracy for cluster properties (need to reach bulk limit)
- Computational cost (scales with cluster size)

Larger ``r2`` → includes more neighbors → more accurate but slower.

Block Recursion for Multi-Orbital Systems
==========================================

For systems with multiple orbitals per atom (s, p, d), use **block recursion:**

Instead of scalar recursion for single orbital, perform recursion on full matrix blocks:

.. math::

   \mathbf{a}_n \in \mathbb{C}^{(n_{\text{orb}} \times n_{\text{orb}})}

Each $a_n$ is now a matrix, $b_n^2$ is a matrix product.

Computational cost: $O(n_{\text{orb}}^3)$ higher per step, but fewer iterations needed.

**Implementation:**

- ``recursion.f90::hop_b()`` - Block version of ``hop()``
- ``recursion.f90::recur_b()`` - Full block recursion
- ``recursion_mod::a_b``, ``b2_b`` - Block coefficient arrays

Comparison: Direct Diagonalization vs. Recursion
=================================================

.. table::
   :align: left

   +-----------------------+------------------------+------------------+
   | Aspect                | Direct Diagonalization | Recursion Method |
   +=======================+========================+==================+
   | Memory                | $O(N^2)$               | $O(N \times L)$  |
   +-----------------------+------------------------+------------------+
   | Time (per energy)     | $O(N^3)$               | $O(N \times L)$  |
   +-----------------------+------------------------+------------------+
   | Spectral resolution   | Fixed by eigenvalues   | Flexible         |
   +-----------------------+------------------------+------------------+
   | Works for disorder    | ✓ (harder)             | ✓ (natural)      |
   +-----------------------+------------------------+------------------+
   | Requires k-mesh       | ✓                      | ✗                |
   +-----------------------+------------------------+------------------+
   | Code complexity       | Low                    | High             |
   +-----------------------+------------------------+------------------+

where $L$ is recursion cutoff (``llsp`` or ``lld``).

Provenance
==========

Main implementation locations:

- **Lanczos recursion:**
  
  - ``source/recursion.f90::hop()`` - Single step
  - ``source/recursion.f90::recur()`` - Full loop
  - ``source/recursion.f90::crecal()`` - Coefficient calculation

- **Chebyshev moments:**
  
  - ``source/recursion.f90::chebyshev_recur_ll()`` - Moment calculation
  - ``source/density_of_states.f90::chebyshev_dos()`` - DOS from moments

- **Green's function:**
  
  - ``source/green.f90::sgreen()`` - Uses recursion coefficients
  - ``source/green.f90::block_green()`` - Block version

References
==========

- Pettifor & Weaire, "The Electronic Structure of Complex Systems" - Recursion method overview
- Lanczos, "Solution of Systems of Linear Equations by Minimized Iterations", J. Res. Natl. Bur. Stand. 49, 33 (1952)
- Weaire et al., "The use of continued fractions in the calculation of electronic structure", Comm. Phys. 1, 25 (1976)

See Also
========

- :ref:`theory/green_functions` - Green's function formalism
- :ref:`theory/scf_cycle` - Integration into SCF procedure
- :ref:`code_structure` - Module organization
