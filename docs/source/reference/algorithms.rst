.. _reference/algorithms:

=======================================
Key Algorithms
=======================================

SCF Loop Algorithm
==================

The self-consistent field cycle is the core of RS-LMTO.

**Pseudocode:**

.. code-block:: text

   for iteration = 1 to nstep
       ! 1. Build Hamiltonian with current potential
       build_bulkham (and build_lsham with SOC)

       ! 2. Recursion, selected by control%recur (default 'block')
       'lanczos'  : recur
       'block'    : recur_b
       'chebyshev': chebyshev_recur

       ! 3. Green function / DOS
       'lanczos'  : sgreen
       'block'    : zsqr, block_green
       'chebyshev': chebyshev_dos_dispatch

       ! 4. Fermi level, band moments (pl, ql)
       bands%calculate_fermi, bands%calculate_moments

       ! 5. Mix moments, then new potential by calctype
       mixpq
       'B': bulkpot   'S': surfpot   'I': imppot   'L': interfacepot

       ! 6. Convergence: delta = RMS difference of old and new packed
       !    (pl, ql) moments before mixing; converged if delta < conv_thr
       !    (and constraint / fixed-U residual checks when active)
       if converged: exit

   ! Loop ends after nstep iterations without a failure return.

**Implementation:**

- ``source/self.f90`` - Main loop and per-iteration dispatch
- ``source/hamiltonian_build.f90::build_bulkham()`` - Step 1
- ``source/recursion_haydock.f90::recur()``, ``recur_b()``;
  ``source/recursion_chebyshev.f90::chebyshev_recur()`` - Step 2
- ``source/green_lanczos.f90::sgreen()``, ``source/green_block.f90::block_green()`` - Step 3
- ``source/bands.f90::calculate_fermi()``, ``calculate_moments()`` - Step 4
- ``source/mix.f90::mixpq()``, ``source/charge.f90::bulkpot()`` - Step 5

Lanczos Recursion
=================

Generate tridiagonal Hamiltonian via Lanczos algorithm.

**Algorithm:**

.. code-block:: text

   a[n], b[n] = lanczos(H, v0, l_max)
       |v1⟩ = v0 / ||v0||
       for n = 1 to l_max
           a[n] ← ⟨vn | H | vn⟩
           |w⟩ ← (H - a[n]) |vn⟩ - b[n] |vn-1⟩
           b[n+1] ← ||w||
           if b[n+1] < ε (breakdown)
               return
           |vn+1⟩ ← |w⟩ / b[n+1]

**Tridiagonal matrix representation:**

.. math::

   H_{\text{tridiag}} = \begin{pmatrix}
       a_1 & b_1 & 0 & \cdots \\
       b_1 & a_2 & b_2 & \cdots \\
       0 & b_2 & a_3 & \cdots \\
       \vdots & \vdots & \vdots & \ddots
   \end{pmatrix}

**Implementation:**

- ``source/recursion_haydock.f90::recur()`` - Full loop
- ``source/recursion_haydock.f90::hop()`` - Single step (compute :math:`a_n`, :math:`b_n`)
- ``source/recursion_haydock.f90::crecal()`` - Coefficient recurrence

Chebyshev Recursion
===================

Compute spectral moments using Chebyshev polynomials.

**Algorithm:**

.. code-block:: text

   μ[n] = chebyshev_moment(H, n, E_min, E_max)
       ! Rescale energy to [-1, 1]
       a ← (E_max - E_min)/(2 - 0.3),  b ← (E_max + E_min)/2
       H_scaled ← (H - b)/a
       
       ! Chebyshev polynomial recursion: T_{n+1} = 2x T_n - T_{n-1}
       |ψ0⟩ ← unit orbital block on the site (SCF and DOS)
       |ψ1⟩ ← H_scaled |ψ0⟩
       
       μ[0] ← ⟨ψ0 | ψ0⟩
       μ[1] ← ⟨ψ0 | H_scaled | ψ0⟩
       
       for n = 2 to n_cheb
           |ψn⟩ ← 2 H_scaled |ψn-1⟩ - |ψn-2⟩
           μ[n] ← ⟨ψ0 | ψn⟩

**DOS from moments:**

.. math::

   \rho(E) = \frac{1}{\pi\sqrt{a^2-(E-b)^2}} \left[ g_0 \mu_0 + 2 \sum_{n \ge 1} g_n \mu_n T_n\!\left(\frac{E-b}{a}\right) \right]

with :math:`g_n` the Jackson kernel (``jackson_kernel``). Random-phase
starting vectors are used only for transport with
``cond_calctype = 'random_vec'``.

**Implementation:**

- ``source/recursion_chebyshev.f90::chebyshev_recur()`` - Moment calculation
- ``source/density_of_states.f90::chebyshev_dos()`` - DOS from moments

Green's Function Evaluation
============================

Compute :math:`G(E) = (E - H + i\eta)^{-1}` using continued fraction.

**From recursion coefficients:**

.. code-block:: text

   G(E) = 1 / (E - a[0] - b[1]² / (E - a[1] - b[2]² / (...)))

**Continued fraction evaluation:**

.. code-block:: text

   G = continued_fraction(E, a[], b[])
       ! Backward recursion (numerically stable)
       u[l_max+1] = 0
       u[l_max] = 1
       for n = l_max down to 1
           u[n-1] ← (E - a[n]) u[n] - b[n]² u[n+1]
       
       G ← 1 / u[0]

**Implementation:**

- ``source/green_lanczos.f90::sgreen()`` - On-site
- ``source/green_lanczos.f90::bgreen()`` - Block (multi-orbital)

Density Integration
===================

The Fermi level is found from the total DOS (``bands%fermi``):

.. code-block:: text

   integral ← 0
   for each pair of mesh intervals (Simpson's rule)
       integral ← integral + h (y[i-1] + 4 y[i] + y[i+1]) / 3
       if integral ≥ qqv (valence count): stop
   E_F ← linear interpolation inside the last interval

There is no temperature in this search. ``calculate_fermi`` runs it only
when ``fix_fermi`` is false and ``calctype`` is 'B' or 'L'.

**Implementation:**

- ``source/bands.f90::calculate_fermi()`` - Fermi level search

Density Mixing
==============

**Linear mixing:**

.. code-block:: text

   q_new ← α q_out + (1-α) q_in
   α ∈ [0, 1] : smaller α -> more conservative
   
   Advantage: Simple, always stable
   Disadvantage: Slow convergence

**Broyden mixing:**

.. code-block:: text

   function broyden_mix(q_in, q_out, history)
       ! Minimize residual over history
       r = q_out - q_in  (residual)
       
       if first iteration
           return α × r + q_in
       
       ! Construct Jacobian estimate from history
       J ← estimate from previous r vectors
       
       ! Predict optimal mix
       Δq_optimal ← -J⁻¹ r
       q_new ← q_in + Δq_optimal
       
       ! Store in history for next iteration
       return q_new

**Implementation:**

- ``source/mix.f90::mixpq()`` - Main routine

Exchange Interaction Calculation
=================================

For each pair :math:`(i, j)` in ``ijpair``, ``calculate_exchange``
integrates

.. math::

   J_{ij} = \frac{10^3}{4\pi} \int^{E_F} dE \,
   \text{Im}\,\text{Tr}\Big[ \Delta_i G_i^0 \Delta_j G_j^0
   - \sum_\alpha \Delta_i G_i^\alpha \Delta_j G_j^\alpha \Big]

with Simpson's rule and the Fermi weight (:math:`G^0`, :math:`G^\alpha`:
the ``Ginmag``/``Gjnmag`` and ``Gix``... ``Gjz`` arrays). See
:doc:`exchange_module` for D, A and the output files.

**Implementation:**

- ``source/exchange.f90::calculate_exchange()`` - Main routine
- Uses ``green.f90`` for inter-site Green's functions

Conductivity Calculation
========================

Chebyshev expansion of the Kubo-Bastin formula:

.. math::

   \sigma_{\alpha\beta}(E) = \frac{16}{\pi \Delta E^2} \sum_{n,m}
   \Gamma_{nm}(E) \mu_{nm}^{\alpha\beta}

followed by a zero-temperature Fermi-weighted cumulative Simpson integral
over energy. See :doc:`conductivity_module`.

**Implementation:**

- ``source/conductivity.f90::calculate_conductivity_tensor()`` - Main routine
- ``source/hamiltonian_build.f90::build_realspace_velocity_operators()`` - Velocity ops

Cluster Geometry Construction
=============================

``lattice%bravais`` builds the cluster from the lattice vectors ``a``
and the repetition counts ``n1``, ``n2``, ``n3``, keeping sites within the
cutoff.

**Implementation:**

- ``source/lattice_lifecycle.f90::bravais()`` - Cluster generation
- ``source/math.f90`` - Distance/geometry utilities

Hamiltonian Construction
========================

Build electronic Hamiltonian from TB parameters.

**Algorithm:**

.. code-block:: text

   H_ij = hamiltonian_element(i, j, lattice, potential, charge)
       
       if i == j (same atom):
           ! On-site: orbital energy + potential
           H_ii ← C_l + V_eff[charge_i]
       
       else (i ≠ j):
           ! Inter-site hopping
           R_ij ← position_j - position_i
           
           ! Distance-dependent hopping
           t_ij ← compute_hopping(R_ij, lmax, potential)
           
           ! Madelung corrections
           M_ij ← madelung_correction(R_ij, charge)
           
           H_ij ← t_ij + M_ij
       
       ! Spin-orbit coupling (if nsp = 2 or 4)
       if include_soc:
           H_ij += L·S term

**Implementation:**

- ``source/hamiltonian_build.f90::build_bulkham()`` - Bulk
- ``source/hamiltonian_build.f90::build_locham()`` - Local
- ``source/hamiltonian_build.f90::build_lsham()`` - Spin-orbit

Provenance
==========

All algorithms documented in:

- Theory sections of this documentation
- Source code comments (especially pseudocode in main routines)
- Published papers cited in theory pages

Key source files:

- ``source/self.f90`` - SCF loop driver
- ``source/recursion.f90`` - Lanczos/Chebyshev recursion
- ``source/green.f90`` - Green's function computation
- ``source/density_of_states.f90`` - DOS and integration
- ``source/mix.f90`` - Density mixing
- ``source/hamiltonian.f90`` - Hamiltonian construction

See Also
========

- :ref:`theory/scf_cycle` - SCF theory
- :ref:`theory/recursion_method` - Recursion theory
- :ref:`theory/green_functions` - Green's function theory
- :doc:`module_overview` - Module organization
