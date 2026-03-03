.. _theory/scf_cycle:

=======================================
Self-Consistent Field (SCF) Cycle
=======================================

Introduction
============

The **Self-Consistent Field (SCF)** cycle is the iterative procedure that solves the 
Kohn-Sham equations of DFT. RS-LMTO adjusts the effective one-electron potential 
until the charge density and potential reach self-consistency.

SCF Loop Overview
=================

The SCF cycle implements:

.. math::

   n_{\text{in}} \rightarrow V_{\text{eff}}(n_{\text{in}}) \rightarrow H(V_{\text{eff}}) \rightarrow n_{\text{out}} \rightarrow \text{check convergence}

**Iterative procedure:**

1. **Input charge density** $n^{(i)}$
2. **Construct effective potential** $V_{\text{eff}}(n^{(i)})$ from:
   
   - Exchange-correlation energy/potential
   - Hartree (Coulomb) potential
   - External potential (fixed nuclei)

3. **Build Hamiltonian** $H$ with potential $V_{\text{eff}}$
4. **Compute electronic structure** (Green's functions, DOS)
5. **Integrate DOS** to get new charge density $n^{(i+1)}$
6. **Mix densities** to improve convergence:
   
   .. math::

      n^{(i+1/2)} = \alpha n^{(i+1)} + (1-\alpha) n^{(i)}

7. **Check convergence** and repeat or stop

Key Equations
=============

**Charge Density:**

From Green's function (DOS integration below Fermi level):

.. math::

   n_i^{(i+1)} = \int_{-\infty}^{E_F} dE \, \rho_i(E) \, f(E - \mu)

where:

- $\rho_i(E)$ is local density of states
- $f(E - \mu)$ is Fermi-Dirac distribution (includes temperature)
- $E_F$ is Fermi energy determined from electron count

**Effective Potential:**

.. math::

   V_{\text{eff}} = V_{\text{XC}} + V_{\text{Hartree}} + V_{\text{ext}}

where:

- $V_{\text{XC}}$ - Exchange-correlation (from XC functional)
- $V_{\text{Hartree}}$ - Coulomb potential from charge distribution
- $V_{\text{ext}}$ - External potential (fixed in ground state)

**Energy Functional:**

Total energy (to be minimized):

.. math::

   E[\{n_i\}, V] = T[G(E)] + E_{\text{XC}}[n] + E_{\text{Hartree}}[n] + \mu N_e

where $T$ is kinetic energy derived from Green's function moments.

Implementation in RS-LMTO
=========================

**Main module:**

``source/self.f90`` contains ``type :: self`` with procedures:

.. code-block:: fortran

   type :: self
      ! ... pointers to other modules (lattice, charge, control, etc.)
   contains
      procedure :: iterate       ! Single SCF iteration
      procedure :: process       ! Full SCF loop
      procedure :: converge_scf  ! Check convergence
   end type self

**Typical SCF loop (pseudocode from ``self.f90::process``):**

.. code-block:: fortran

   do iteration = 1, max_iterations
      
      ! 1. Build Hamiltonian with current potential
      call hamiltonian%build_bulkham(lattice, charge, ...)
      
      ! 2. Compute recursion coefficients
      call recursion%recur(hamiltonian, ...)
      
      ! 3. Compute Green's functions
      call green%bgreen(recursion, energy, ...)
      
      ! 4. Integrate DOS for charge density
      call density_of_states%density(green, recursion, ...)
      
      ! 5. Update charge and potential
      call charge%calculate_density(dos, ...)
      call charge%calculate_potential(...)
      
      ! 6. Mix charge densities
      call mix%mix_charge(charge, ...)
      
      ! 7. Check convergence
      if (converge_scf()) exit
      
   end do

SCF Convergence Criteria
========================

**Typical Convergence Checks:**

1. **Charge density difference:**
   
   .. math::

      \Delta q = \sum_i |q_i^{(i+1)} - q_i^{(i)}| < \text{tolerance}

   Control: ``control%dq_tol`` (default: $10^{-5}$)

2. **Potential difference:**
   
   .. math::

      \Delta V = \sum_i |V_i^{(i+1)} - V_i^{(i)}| < \text{tolerance}

3. **Energy convergence:**
   
   .. math::

      \Delta E = |E^{(i+1)} - E^{(i)}| < \text{tolerance}

4. **Force/torque on atoms** (if applicable)

All criteria must be satisfied. Convergence typically requires 10-100 iterations.

Density Mixing Strategies
=========================

Pure **linear mixing** (often too slow):

.. math::

   n^{(i+1/2)} = \alpha n_{\text{new}}^{(i+1)} + (1 - \alpha) n_{\text{old}}^{(i)}

where mixing parameter $\alpha \in (0, 1]$ controls step size.

**Broyden Mixing** (more efficient):

Uses information from multiple previous iterations to predict better mix:

.. math::

   n^{(i+1/2)} = n^{(i)} + \alpha \Delta n^{(i)} + \sum_j c_j^{(i)} u_j^{(i)}

where:

- $\Delta n^{(i)} = n_{\text{new}} - n_{\text{old}}$
- $u_j$ are residual vectors from previous iterations
- $c_j$ are optimized coefficients (solve least-squares problem)

**Implementation in code:**

``source/mix.f90`` - contains mixing algorithms

.. code-block:: fortran

   type :: mix
      character(len=sl) :: mixing_type  ! 'linear' or 'broyden'
      real(rp) :: alpha                 ! Mixing parameter
      integer :: n_history              ! History length for Broyden
   contains
      procedure :: mix_charge           ! Apply mixing
   end type mix

**Control Parameters:**

See :ref:`keywords/scf_settings` for input keywords:

- ``mixing`` - Mixing type
- ``alpha`` - Mixing parameter
- ``broyden_history`` - Number of previous iterations to keep

Convergence Acceleration
========================

**Techniques used:**

1. **Variable mixing parameter** - Increase $\alpha$ as convergence improves
2. **Residual norm weighting** - Weight recent residuals more heavily
3. **Extrapolation** - Predict next density from trend

**When convergence stalls:**

- Reduce mixing parameter $\alpha$
- Increase Broyden history length
- Check for oscillations between states
- Ensure recursion cutoff is adequate

Fermi Level Adjustment
======================

During SCF, the **Fermi level** must be adjusted to maintain correct electron count:

**Target electron count:**

For each atom type and spin:

.. math::

   N_e = \sum_i n_i

**Fermi level search:**

.. math::

   \mu = E_F(\text{search to match } N_e)

**Methods:**

1. **Bisection** - Slow but reliable
2. **Newton-Raphson** - Faster, requires DOS derivative
3. **Fixed Fermi** (``fix_fermi = .true.``) - Skip adjustment (useful for impurities)

**Implementation:**

``energy.f90::find_fermi()`` performs energy mesh and searches for correct $E_F$.

Temperature Effects
===================

At **finite temperature**, use Fermi-Dirac distribution:

.. math::

   f(E) = \frac{1}{e^{(E - \mu)/k_B T} + 1}

**Effects on charge density:**

- Smears out Fermi surface
- Electrons excite above $E_F$
- Changes magnetic moments and stability

**Control:**

Keyword ``temperature`` in input (Kelvin or Ry). Default: T=0 (zero temperature).

Energy Calculation and Forces
==============================

The **total energy** at each SCF iteration includes:

.. math::

   E_{\text{tot}} = T[G] + E_{\text{XC}}[n] + E_{\text{Hartree}}[n] - \mu N_e + \text{const}

where $T[G]$ is kinetic energy from Green's function:

.. math::

   T = \int dE \, E \, \rho(E) \, f(E)

**Forces on atoms** (if doing geometry optimization):

.. math::

   F_I = -\frac{\partial E_{\text{tot}}}{\partial \mathbf{R}_I}

Computed via Hellmann-Feynman theorem (handled by ``force`` module if available).

Magnetic Convergence
====================

For **magnetic systems** (``nsp > 1``):

- Separate SCF loops for spin-up and spin-down
- Or: fully magnetized iterative solution
- Convergence checked on **magnetization** as well as charge

**Constraints:**

- **Fixed magnetization** - Constrain $m_z$ per atom (exchange-field approach)
- **Constrained moments** - For magnetic impurity problems

Implemented via exchange field: additional $\pm B_{\text{exch}} \sigma_z$ term in Hamiltonian.

Non-Collinear Magnetism
~~~~~~~~~~~~~~~~~~~~~~~

For **non-collinear** calculations (``nsp = 3`` or ``4``):

- Magnetization is a **3D vector** per atom
- Magnetic moment rotates during SCF
- More complex convergence behavior

**Acceleration:**

Use torque extrapolation or spin dynamics predictor-corrector.

When SCF Fails to Converge
===========================

**Common issues:**

1. **Too aggressive mixing** - Reduce ``alpha``
2. **Inadequate recursion depth** - Increase ``llsp``, ``lld``
3. **Energy mesh too coarse** - Increase ``channels_ldos``
4. **Unstable Hamiltonian** - Check element/potential parameters
5. **Bad initial guess** - Start from converged nearby configuration

**Debugging strategies:**

- Monitor energy evolution (should decrease monotonically)
- Plot charge density per iteration
- Check for charge oscillations
- Increase output verbosity (``verbose = .true.``)

Example SCF Convergence
=======================

Typical convergence for bulk Si (from ``example/bulk/Si/``):

.. code-block:: text

   Iteration   ΔQ (mRy)    ΔE (mRy)    Wall time (s)
   --------    --------    --------    --------
   1           0.523       -0.287      2.1
   2           0.156        0.089      1.9
   3           0.042        0.024      1.8
   4           0.009       -0.002      1.8
   5           0.002        0.000      1.8
   SCF converged in 5 iterations

Provenance
==========

Main code locations:

- **SCF loop:** ``source/self.f90::process()``
- **Single iteration:** ``source/self.f90::iterate()``
- **Convergence check:** ``source/self.f90::converge_scf()``
- **Charge integration:** ``source/density_of_states.f90::density()``
- **Density mixing:** ``source/mix.f90``
- **Potential update:** ``source/charge.f90::calculate_potential()``

Key control parameters: ``source/control.f90`` (max_iterations, convergence tolerances, etc.)

See Also
========

- :ref:`theory/green_functions` - Green's function calculation in SCF
- :ref:`theory/recursion_method` - Recursion method used in each iteration
- :ref:`keywords/scf_settings` - Input control parameters
- :ref:`user_guide/input_files` - SCF-related namelists
