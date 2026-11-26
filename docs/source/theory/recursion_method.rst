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
- :math:`O(N)` memory and :math:`O(N)` per step (for sparse matrices)
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

expanded using coefficients :math:`a_0, a_1, \ldots, a_N` and :math:`b_1, b_2, \ldots, b_N`.

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
- Requires energy rescaling to :math:`[-1, 1]` interval
- Better numerical stability
- More suitable for DOS and moment calculations

**RS-LMTO uses both:**

- **Lanczos** for direct Green's function computation (``recursion.f90::hop()``)
- **Chebyshev** for DOS and moment integrations (``recursion.f90::chebyshev_*``)

Chebyshev Moments and DOS
=========================

Instead of computing :math:`\rho(E)` directly from Green's function, RS-LMTO computes **Chebyshev moments**:

.. math::

   \mu_n = \langle 0 | T_n(\tilde H) | 0 \rangle, \qquad \tilde H = (H - b)/a

for each recursion site and orbital block :math:`|0\rangle` (not a trace over
the cluster). The DOS is the kernel-damped sum

.. math::

   \rho(E) = \frac{1}{\pi \sqrt{a^2 - (E-b)^2}}
   \left[ g_0 \mu_0 + 2 \sum_{n \ge 1} g_n \mu_n T_n\!\left(\frac{E-b}{a}\right) \right]

where:

- :math:`T_n` are Chebyshev polynomials
- :math:`g_n` is the Jackson kernel
- :math:`a`, :math:`b` are the scaling parameters given under "Energy Rescaling" below
- the number of moments is :math:`2\,\texttt{lld} + 2`

**Advantages:**

- Very stable numerically
- Fast convergence (Chebyshev polynomials are optimal)
- No explicit density matrix needed

Spectral positivity
~~~~~~~~~~~~~~~~~~~

The real-axis on-site Chebyshev reconstruction is regularized before DOS,
Fermi-level, and radial charge/moment integrations. Each orbital's 2x2 spin
spectral-density block is projected to be positive semidefinite: negative
diagonal spectral weight is set to zero and spin coherence is reduced only
when it exceeds the physical bound. This removes small negative-DOS artifacts
from finite-order or reduced-precision paths while preserving the dispersive
part of the Green function. The correction is reported in the run log.

This regularization applies only to the local real-energy DOS route; complex-
energy and intersite Green functions are not changed. A large or frequent
correction indicates that the Chebyshev bounds/order should be checked. CUDA
runs can also request ``gpu_precision = 'fp64'`` in ``&control``.

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
      complex(rp), allocatable :: psi(:,:)      ! Lanczos wavefunctions
      complex(rp), allocatable :: mu_n(:,:,:,:) ! Chebyshev moments
   contains
      procedure :: hop         ! Lanczos hop (single recursion step)
      procedure :: hop_b       ! Block Lanczos
      procedure :: recur       ! Full Lanczos recursion
      procedure :: recur_b     ! Block recursion
      procedure :: chebyshev_recur_ll ! One Chebyshev recursion step
      ! ... more procedures
   end type recursion

**Key Procedures:**

1. **``hop()`` and ``hop_b()``**
   
   Single recursion step: compute :math:`a_n`, :math:`b_n`, and next Lanczos vector.

2. **``recur()`` and ``recur_b()``**
   
   Full recursion loop: iterate to the recursion depth ``lld``.

3. **``chebyshev_recur()``**
   
   Chebyshev moment calculation, per site (``chebyshev_recur_ll()`` is one
   recursion step).

4. **``chebyshev_dos()``**
   
   Compute DOS from Chebyshev moments via polynomial expansion.

Recursion Cutoff Parameters
===========================

**Control parameters** (from ``control.f90``):

- **``lld``** - Recursion depth for all orbitals (typical: 50-100)
- **``llsp``** - Enters only array dimensions (``max(llsp, lld)``); it does
  not change the recursion

``lld`` sets the number of coefficients computed.

**Trade-off:**

- Higher cutoff -> better accuracy but more computation
- Lower cutoff -> faster but coarser spectral resolution

Typical values balance accuracy and speed; see examples in ``example/`` directory.

Random Vectors (Conductivity Only)
==================================

The DOS/SCF Chebyshev moments are computed per site, without random vectors.
The Kubo-Bastin conductivity moments (``compute_moments_stochastic`` in
``source/recursion_transport.f90``) use random starting vectors when
``cond_calctype = 'random_vec'`` (default ``'per_type'``); ``random_vec_num``
in ``&control`` sets their number.

Numerical Aspects
=================

**Broadening Parameter (eta):**

Added as :math:`i\eta` to avoid singularities:

.. math::

   G(E) = (E - H + i\eta)^{-1}

Controls spectral resolution. Typical: :math:`\eta \sim 10^{-3}` to :math:`10^{-2}` Ry.

**Energy Rescaling for Chebyshev:**

To use Chebyshev polynomials (:math:`T_n : [-1,1] \to \mathbb{R}`), energies are rescaled:

.. math::

   E_{\text{scaled}} = \frac{E - b}{a}, \qquad
   a = \frac{E_{\max} - E_{\min}}{1.7}, \qquad
   b = \frac{E_{\max} + E_{\min}}{2}

where :math:`E_{\min}`, :math:`E_{\max}` are the Hamiltonian bounds when
they are set (``bounds%e_max > bounds%e_min``), otherwise ``energy_min`` and
``energy_max`` from ``&energy`` (``resolve_chebyshev_window``).

Chebyshev Bounds Workflow
=========================

Chebyshev scaling in production recursion is driven by Hamiltonian bounds and implemented through ``source/bounds.f90``.

- ``bounds_algorithm='gershgorin'``: use real-space Gershgorin bounds from assembled ``ee``.
- ``bounds_algorithm='hgamma'`` or ``'sturm'``: use ``H(Gamma)`` spectrum estimate (via ``compute_spectrum_bounds``), with Gershgorin fallback.
- ``bounds_algorithm='both'`` (or ``'hybrid'``): conservative envelope of both estimators.
- ``bounds_scaling``: symmetric expansion factor applied to the selected interval.
- Default is ``bounds_algorithm='none'``: keep legacy behavior and use ``energy_min``/``energy_max`` from ``&energy``.

This is a recursion/Chebyshev feature, not a k-space-only feature.

Example (``&hamiltonian``):

.. code-block:: fortran

   &hamiltonian
      bounds_algorithm = 'both'
      bounds_scaling = 1.05
   /

**Convergence Criteria:**

- Check when Chebyshev moments drop below numerical noise
- Typical convergence: :math:`\mu_n < 10^{-12}` for machine precision

Hopping Region and Cluster Truncation
=====================================

In real-space formulation, the recursion operates within a **hopping region** defined by:

1. **Central atom(s)** where properties are computed
2. **Neighboring atoms** within radius :math:`\sqrt{r2}` (``r2`` is the squared cutoff radius, ``r2 = ct**2``)

Beyond the hopping region, contributions are neglected (Hamiltonian is zero).

**Determines:**

- Accuracy for cluster properties (need to reach bulk limit)
- Computational cost (scales with cluster size)

Larger ``r2`` -> includes more neighbors -> more accurate but slower.

Block Recursion for Multi-Orbital Systems
==========================================

For systems with multiple orbitals per atom (s, p, d), use **block recursion:**

Instead of scalar recursion for single orbital, perform recursion on full matrix blocks:

.. math::

   \mathbf{a}_n \in \mathbb{C}^{(n_{\text{orb}} \times n_{\text{orb}})}

Each :math:`a_n` is now a matrix, :math:`b_n^2` is a matrix product.

Computational cost: :math:`O(n_{\text{orb}}^3)` higher per step, but fewer iterations needed.

**Implementation:**

- ``recursion.f90::hop_b()`` - Block version of ``hop()``
- ``recursion.f90::recur_b()`` - Full block recursion
- ``recursion_mod::a_b``, ``b2_b`` - Block coefficient arrays

Inter-site Green functions: pairs and one-shot
==============================================

Exchange, damping and inertia need the inter-site blocks
:math:`G_{ij}(z) = P_i^\dagger (z - H)^{-1} P_j`, where :math:`P_i` embeds the
:math:`n_b \times n_b` orbital-spin basis of site :math:`i`. The control key
``intersite_mode`` selects one of two constructions; both fill the same
``green%gij`` and ``green%gji`` arrays, so the consumers are unchanged.

**Pairs route** (``intersite_mode = 'pairs'``). For each pair, four recursions
start from :math:`|\pm\rangle = (|i\rangle \pm |j\rangle)/\sqrt{2}` and
:math:`|\pm i\rangle = (|i\rangle \pm i|j\rangle)/\sqrt{2}`. Their on-site Green
functions :math:`G_1, \ldots, G_4` combine to

.. math::

   G_{ij} = \tfrac12\left[(G_1 - G_2) - i(G_3 - G_4)\right], \qquad
   G_{ji} = \tfrac12\left[(G_1 - G_2) + i(G_3 - G_4)\right].

The complex seeds supply the s-p and p-d elements, which cancel in the real
combination alone (Inoue, Okada and Ohta 1993). A self-pair uses one recursion
from :math:`|i\rangle`.

**One-shot route** (``intersite_mode = 'one_shot'``). One recursion starts from
:math:`P_i` alone (the unit block at site :math:`i`). The pair list is grouped
by its first atom (``source/intersite_request.f90``); during the recursion the
block of every target :math:`j` of that seed is read off, so all pairs with the
same first atom share one recursion. This gives :math:`G_{ji}` directly,
including the s-p and p-d elements, and the reverse block follows from
Hermiticity of :math:`H`:

.. math::

   G_{ij}(z) = \left[G_{ji}(\bar z)\right]^\dagger .

*Chebyshev.* With :math:`V_0 = P_i`, :math:`V_1 = \tilde H V_0` and
:math:`V_{m+1} = 2\tilde H V_m - V_{m-1}`, the moments are
:math:`\mu_m(j,i) = P_j^\dagger V_m`, the :math:`j` block of :math:`V_m`, and
:math:`\mu_m(i,j) = \mu_m(j,i)^\dagger` exactly. Bra and ket differ, so the
product identity that doubles the moments of a single seed does not apply:
:math:`2\,\mathrm{lld}+2` moments need :math:`2\,\mathrm{lld}+1` Hamiltonian
applications (the pairs route needs :math:`4(\mathrm{lld}+1)` per pair). The
moments are the same as in the pairs route, so the kernel, the reconstruction
and the results are the same up to rounding. With ``cheb_lp`` the linear
prediction is fitted once to each block :math:`\mu_m(j,i)`; its real shared
coefficients keep :math:`\mu_m(i,j) = \mu_m(j,i)^\dagger` for the extended
moments.

*Block Haydock.* Block Lanczos from :math:`Q_1 = P_i` gives normalized blocks
:math:`Q_m` and coefficients :math:`A_m` and Hermitian :math:`B_m` (the square
root of the stored :math:`B_m^2`), i.e. a block-tridiagonal :math:`T_M` with
:math:`M = \mathrm{lld}-1` blocks. At each step the target blocks
:math:`S_m(j) = P_j^\dagger Q_m` are stored. Then

.. math::

   G_{ji}(z) = \sum_{m=1}^{M} S_m(j)\, X_m(z), \qquad
   X(z) = \text{first block column of } (z - T_M - \Sigma_M(z))^{-1},

solved by a backward Schur recursion (``source/block_tridiag.f90``).
:math:`\Sigma_M` is the same terminal self-energy as in the pairs route and
the on-site recursion (``sym_term`` selects the orbital-independent or
orbital-dependent semicircle; it acts on the last block only). For the reverse
block the advanced terminator :math:`\Sigma_M(z)^\dagger` is used. The
recursion runs on the CPU or the GPU; the projected solve, one per seed and
energy for all its targets, runs on the CPU.

**Accuracy.** The one-shot and pairs Chebyshev routes give the same results
up to rounding (4e-13 to 6e-13 relative in :math:`J_{ij}` on bcc Fe). The block
Haydock routes project onto different Krylov spaces at finite depth, from
:math:`P_i` and from :math:`(P_i + cP_j)/\sqrt{2}`, so off-site
:math:`G_{ij}`, and hence :math:`J_{ij}` and :math:`D_{ij}`, differ between them
and converge only with ``lld``. On the bcc Fe contour test case the first-neighbour
:math:`J_{ij}` differed by about 2 %; other pairs and the real-energy mesh
showed much larger differences, and the pairs route itself varied strongly with
``lld`` there. Self-pairs use the same Krylov space in both routes and agree
to rounding.

**Cost.** Per distinct first atom: one recursion of :math:`2\,\mathrm{lld}+1`
(Chebyshev) or :math:`\mathrm{lld}` (Haydock) Hamiltonian applications, shared
by all its targets, against four recursions per pair in the pairs route. On
bcc Fe with eight pairs sharing one first atom, the Chebyshev run took 15 s
against 124 s. Memory is currently the same as in the pairs route.

**Moment-space exchange integral** (``exchange_integration = 'moments'``, ``recur = 'chebyshev'``). Every exchange integrand is a sum of terms
:math:`\mathrm{Tr}[\Delta_i(E)\, X(E)\, \Delta_j(E)\, Y(E)]`, with :math:`X` built linearly from :math:`G_{ij}, G_{ji}` and
likewise :math:`Y`. With :math:`G(E)` the kernel-damped series of the stored moments :math:`\mu_n`, the energy integral is a
bilinear form in the moments,

.. math::

   \int f(E)\, E^p\, \mathrm{Tr}[A\, X(E)\, B\, Y(E)]\, dE = \sum_m \mathrm{Tr}[A\, X_m\, B\, \tilde Y^{(p)}_m], \qquad
   \tilde Y^{(p)}_m = \sum_n c^{(p)}_{mn} Y_n,

where :math:`c^{(p)}_{mn}` contains the kernel factors of :math:`m` and :math:`n` and the weight
:math:`-\int f(E) E^p e^{-i(m+n)\theta}/(a^2\sin^2\theta)\, dE`, :math:`E = a\cos\theta + b`. The exchange splitting
:math:`\Delta(E)` is linear in :math:`E`, :math:`\Delta = \alpha + \beta E`, so each pair needs the weights with
:math:`p = 0, 1, 2` (``use_onsite_delta``: :math:`p = 0` only). The moment sums are contracted as traces of all pairs of the 12 families of :math:`X` and :math:`Y`, from which every output is assembled. The limits are :math:`[E_{min}, E_F]` for :math:`T = 0` and :math:`f(E)` is the
Fermi function for :math:`T > 0`. Damping and inertia use :math:`G` and the analytic :math:`d^2G/dE^2` of the series at :math:`E_F`; the energy-resolved files are not written in this mode.

**Moment-space SCF integrals** (``scf_integration = 'moments'``, ``recur = 'chebyshev'``). The on-site integrals keep their
formulas :math:`\int f(E) E^p y(E)\,dE` with :math:`y` built from :math:`G`; only the points and the rule change. The
points are the Gauss-Legendre nodes in :math:`\theta` of ``cheb_energy_nodes`` on :math:`[E_{min}, E_F]` (Fermi function in the
weights for ``temperature > 0``), so :math:`\int f E^p y\,dE \approx \sum_q w_q e_q^p y(e_q)`. The Fermi level solves
:math:`N(E_F) = -\tfrac{1}{\pi}\mathrm{Im}\sum_m t_m \sum_q w_q g_m(e_q) = N_{val}`, with :math:`t_m` the diagonal on-site moments summed over atoms
and spins, :math:`g_m` the kernel-damped basis functions of the series.

Comparison: Direct Diagonalization vs. Recursion
=================================================

.. list-table::
   :header-rows: 1
   :align: left

   * - Aspect
     - Direct Diagonalization
     - Recursion Method
   * - Memory
     - :math:`O(N^2)`
     - :math:`O(N \times L)`
   * - Time (per energy)
     - :math:`O(N^3)`
     - :math:`O(N \times L)`
   * - Spectral resolution
     - Fixed by eigenvalues
     - Flexible
   * - Works for disorder
     - [OK] (harder)
     - [OK] (natural)
   * - Requires k-mesh
     - [OK]
     - [X]
   * - Code complexity
     - Low
     - High

where :math:`L` is the recursion depth ``lld``.

Provenance
==========

Main implementation locations:

- **Lanczos recursion:**
  
  - ``source/recursion.f90::hop()`` - Single step
  - ``source/recursion.f90::recur()`` - Full loop
  - ``source/recursion.f90::crecal()`` - Coefficient calculation

- **Chebyshev moments:**
  
  - ``source/recursion.f90::chebyshev_recur_ll()`` - One Chebyshev recursion step
  - ``source/density_of_states.f90::chebyshev_dos()`` - DOS from moments

- **Green's function:**
  
  - ``source/green.f90::sgreen()`` - Uses recursion coefficients
  - ``source/green.f90::block_green()`` - Block version

References
==========

- Pettifor & Weaire, "The Electronic Structure of Complex Systems" - Recursion method overview
- Lanczos, "Solution of Systems of Linear Equations by Minimized Iterations", J. Res. Natl. Bur. Stand. 49, 33 (1952)
- Weaire et al., "The use of continued fractions in the calculation of electronic structure", Comm. Phys. 1, 25 (1976)
- J. Inoue and Y. Ohta, J. Phys. C: Solid State Phys. 20, 1947 (1987) - Block recursion
- J. Inoue, A. Okada and Y. Ohta, J. Phys.: Condens. Matter 5, L465 (1993) - Block recursion with complex starting vectors (s-p and p-d elements of :math:`G_{ij}`)

See Also
========

- :ref:`theory/green_functions` - Green's function formalism
- :ref:`theory/scf_cycle` - Integration into SCF procedure
- :ref:`code_structure` - Module organization
