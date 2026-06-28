.. _reference/algorithms:

=======================================
Key Algorithms
=======================================

SCF Loop Algorithm
==================

The self-consistent field cycle is the core of RS-LMTO.

**Pseudocode:**

.. code-block:: text

   function scf_loop(control, lattice, hamiltonian, charge, mixing)
       for iteration = 1 to max_iterations
           ! 1. Build Hamiltonian with current potential
           H ← build_hamiltonian(lattice, charge, control)
           
           ! 2. Compute recursion coefficients
           a, b² ← recursion_lanczos(H, llsp, lld)
           
           ! 3. Calculate Green's functions
           G(E) ← green_from_recursion(a, b², E)
           
           ! 4. Integrate DOS below Fermi energy
           ρ(E) ← -Im[G(E)] / π
           n_new ← ∫ dE ρ(E) f(E - E_F)
           
           ! 5. Update charge and potential
           charge_new ← n_new
           V_eff ← V_xc[charge_new] + V_Hartree[charge_new]
           
           ! 6. Mix charge densities
           charge ← (1-α)·charge + α·charge_new   (linear)
                  or broyden_mix(charge, charge_new, history)
           
           ! 7. Check convergence
           ΔQ ← || charge - charge_old ||
           if ΔQ < dq_tol
               return success
           
       return failure (max iterations reached)

**Implementation:**

- ``source/self.f90::process()`` - Main loop
- ``source/self.f90::iterate()`` - One iteration
- ``source/hamiltonian.f90::build_bulkham()`` - Step 1
- ``source/recursion.f90::recur()`` - Step 2
- ``source/green.f90::bgreen()`` - Step 3
- ``source/density_of_states.f90::density()`` - Step 4
- ``source/charge.f90::calculate_potential()`` - Step 5
- ``source/mix.f90::mix_charge()`` - Step 6

Lanczos Recursion
=================

Generate tridiagonal Hamiltonian via Lanczos algorithm.

**Algorithm:**

.. code-block:: text

   a[n], b[n] = lanczos(H, v0, l_max)
       |v1⟩ = v0 / ||v0||
       for n = 1 to l_max
           a[n] ← ⟨vn | H | vn⟩
           |w⟩ ← (H - a[n]) |vn⟩ - b[n]² |vn-1⟩
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

- ``source/recursion.f90::recur()`` - Full loop
- ``source/recursion.f90::hop()`` - Single step (compute $a_n$, $b_n$)
- ``source/recursion.f90::crecal()`` - Coefficient recurrence

Chebyshev Recursion
===================

Compute spectral moments using Chebyshev polynomials.

**Algorithm:**

.. code-block:: text

   μ[n] = chebyshev_moment(H, n, E_min, E_max)
       ! Rescale energy to [-1, 1]
       H_scaled ← 2(H - E_min)/(E_max - E_min) - 1
       
       ! Chebyshev polynomial recursion: T_{n+1} = 2x T_n - T_{n-1}
       |ψ0⟩ ← random vector
       |ψ1⟩ ← H_scaled |ψ0⟩
       
       μ[0] ← ⟨ψ0 | ψ0⟩
       μ[1] ← ⟨ψ0 | H_scaled | ψ0⟩
       
       for n = 2 to n_cheb
           |ψn⟩ ← 2 H_scaled |ψn-1⟩ - |ψn-2⟩
           μ[n] ← ⟨ψ0 | ψn⟩

**DOS from moments:**

.. math::

   \rho(E) = \frac{1}{\pi\sqrt{1-(E/W)^2}} \sum_n T_n(E/W) \mu_n

**Implementation:**

- ``source/recursion.f90::chebyshev_recur_ll()`` - Moment calculation
- ``source/density_of_states.f90::chebyshev_dos()`` - DOS from moments

Green's Function Evaluation
============================

Compute $G(E) = (E - H + i\eta)^{-1}$ using continued fraction.

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

- ``source/green.f90::sgreen()`` - On-site
- ``source/green.f90::bgreen()`` - Block (multi-orbital)

Density Integration
===================

Convert DOS to charge density by integrating below Fermi level.

**Algorithm:**

.. code-block:: text

   n = integrate_dos(ρ(E), E_F, T)
       ! Fermi-Dirac distribution
       f(E) = 1 / (exp((E - E_F) / k_B T) + 1)
       
       n ← ∫_{-∞}^{∞} dE ρ(E) f(E)
       
       ! Numerical integration (trapezoidal rule on energy mesh)
       n = 0
       for i = 1 to n_energy
           dE = energy_mesh[i] - energy_mesh[i-1]
           weight = f(energy_mesh[i])
           n += ρ[i] × weight × dE

**Fermi level determination:**

.. code-block:: text

   E_F = find_fermi(ρ(E), N_target, T)
       ! Binary search or Newton-Raphson
       do while |n(E_F) - N_target| > tolerance
           dEdq ← dN/dE_F from DOS derivative
           E_F ← E_F - (n(E_F) - N_target) / dEdq

**Implementation:**

- ``source/density_of_states.f90::density()`` - Integration
- ``source/energy.f90::find_fermi()`` - Fermi level search

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

- ``source/mix.f90::mix_charge()`` - Main routine

Exchange Interaction Calculation
=================================

Compute Heisenberg exchange using :math:`J_{ij} = \langle i | H_\text{eff} | j \rangle`.

**Algorithm:**

.. code-block:: text

   J[i,j] = exchange_interaction(G_ij, ε_F)
       ! For each pair of atoms i, j
       
       ! 1. Compute inter-site Green's function
       G_ij(E) ← green_function(i, j)
       
       ! 2. Extract exchange coupling via:
       J_ij ← ∫ dE [Im(G_ij⁺⁻(E)) - Im(G_ij⁻⁺(E))] Θ(E_F - E)
       
       or via magnetic response:
       χ_ij ← derivative of Green function
       J_ij ← χ_ij²
   
   ! Output list of J vs. distance

**Implementation:**

- ``source/exchange.f90::calculate()`` - Main routine
- Uses ``green.f90`` for inter-site Green's functions

Conductivity Calculation
========================

Compute conductivity tensor using Kubo-Greenwood formula.

**Algorithm:**

.. code-block:: text

   σ_αβ(E) = (e²/2π) Im ⟨v_α G(E) v_β G(E)⟩
   
   where v_α = ∂H/∂k_α (velocity operator)
   
   For each energy point:
   1. Compute velocity operators from Hamiltonian derivatives
   2. Construct: v_α G(E) v_β G(E)
   3. Take imaginary part
   4. Multiply by (e²/2π)
   
   Transport coefficient at temperature T:
   σ(T) = ∫ dE σ(E) (-∂f/∂E)

**Implementation:**

- ``source/conductivity.f90::calculate()`` - Main routine
- ``source/hamiltonian.f90::build_velocity_operators()`` - Velocity ops

Cluster Geometry Construction
=============================

Build atomic cluster from lattice parameters.

**Algorithm:**

.. code-block:: text

   function build_cluster(alat, nx, ny, nz, r2)
       atoms = []
       
       ! Generate atoms in supercell
       for i = 1 to nbulk
           for ix = -nx/2 to nx/2
               for iy = -ny/2 to ny/2
                   for iz = -nz/2 to nz/2
                       ! Lattice position
                       R = (ix, iy, iz) × alat + R_basis[i]
                       
                       ! Check if within cutoff
                       if |R - R_center| < sqrt(r2)
                           atoms.append(R)
       
       return atoms

**Implementation:**

- ``source/lattice.f90::setup_cluster()`` - Cluster generation
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

- ``source/hamiltonian.f90::build_bulkham()`` - Bulk
- ``source/hamiltonian.f90::build_locham()`` - Local
- ``source/hamiltonian.f90::build_lsham()`` - Spin-orbit

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
