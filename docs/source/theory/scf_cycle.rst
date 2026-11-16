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

1. **Input charge density** :math:`n^{(i)}`
2. **Construct effective potential** :math:`V_{\text{eff}}(n^{(i)})` from:
   
   - Exchange-correlation energy/potential
   - Hartree (Coulomb) potential
   - External potential (fixed nuclei)

3. **Build Hamiltonian** :math:`H` with potential :math:`V_{\text{eff}}`
4. **Compute electronic structure** (Green's functions, DOS)
5. **Integrate DOS** to get new charge density :math:`n^{(i+1)}`
6. **Mix densities** to improve convergence:
   
   .. math::

      n^{(i+1/2)} = \alpha n^{(i+1)} + (1-\alpha) n^{(i)}

7. **Check convergence** and repeat or stop

Key Equations
=============

**Charge Density:**

From Green's function (DOS integration below Fermi level):

.. math::

   n_i^{(i+1)} = \int_{-\infty}^{E_F} dE \, \rho_i(E)

where:

- :math:`\rho_i(E)` is local density of states
- :math:`E_F` is Fermi energy determined from electron count (zero temperature)

**Effective Potential:**

.. math::

   V_{\text{eff}} = V_{\text{XC}} + V_{\text{Hartree}} + V_{\text{ext}}

where:

- :math:`V_{\text{XC}}` - Exchange-correlation (from XC functional)
- :math:`V_{\text{Hartree}}` - Coulomb potential from charge distribution
- :math:`V_{\text{ext}}` - External potential (fixed in ground state)

**Energy Functional:**

Total energy (to be minimized):

.. math::

   E[\{n_i\}, V] = T[G(E)] + E_{\text{XC}}[n] + E_{\text{Hartree}}[n] + \mu N_e

where :math:`T` is kinetic energy derived from Green's function moments.

Implementation in RS-LMTO
=========================

**Main module:**

``source/self.f90`` contains ``type :: self`` with procedures:

.. code-block:: fortran

   type :: self
      ! ... pointers to other modules (lattice, charge, control, etc.)
   contains
      procedure :: run           ! SCF loop (the convergence test is inside it)
      procedure :: report        ! Final report
   end type self

**Typical SCF loop (simplified pseudocode of ``self.f90::run``):**

.. code-block:: fortran

   do iteration = 1, nstep
      
      ! 1. Build Hamiltonian with current potential
      call hamiltonian%build_bulkham(lattice, charge, ...)
      
      ! 2. Compute recursion coefficients
      call recursion%recur(hamiltonian, ...)
      
      ! 3. Compute Green's functions
      call green%bgreen(recursion, energy, ...)
      
      ! 4. Integrate DOS for charge density
      call density_of_states%density(green, recursion, ...)
      
      ! 5. Update potential
      call charge%bulkpot()
      
      ! 6. Mix charge densities
      call mix%mixpq(...)
      
      ! 7. Check convergence
      if (mix%delta < conv_thr) exit    ! 'Converged!' in the log
      
   end do

SCF Convergence Criteria
========================

The run is converged when ``mix%delta`` (the RMS difference between old and
new moments, printed as ``Diff``) is below ``conv_thr`` in ``&self`` (default
:math:`0.5\times 10^{-8}`). When active, the fixed-U residual test and the
constraint-angle test (``constraints_enable``) must also pass
(``self%is_converged`` in ``source/self.f90``). There is no potential, energy
or force test.

Density Mixing Strategies
=========================

Pure **linear mixing** (often too slow):

.. math::

   n^{(i+1/2)} = \alpha n_{\text{new}}^{(i+1)} + (1 - \alpha) n_{\text{old}}^{(i)}

where mixing parameter :math:`\alpha \in (0, 1]` controls step size.

**Broyden Mixing** (more efficient):

Uses information from multiple previous iterations to predict better mix:

.. math::

   n^{(i+1/2)} = n^{(i)} + \alpha \Delta n^{(i)} + \sum_j c_j^{(i)} u_j^{(i)}

where:

- :math:`\Delta n^{(i)} = n_{\text{new}} - n_{\text{old}}`
- :math:`u_j` are residual vectors from previous iterations
- :math:`c_j` are optimized coefficients (solve least-squares problem)

**Implementation in code:**

``source/mix.f90`` - contains mixing algorithms

.. code-block:: fortran

   type :: mix
      character(len=7) :: mixtype       ! 'linear' or 'broyden'
      real(rp) :: beta                  ! Mixing parameter
   contains
      procedure :: mixpq                ! Mix the Q and P moments
   end type mix

**Control Parameters:**

See :ref:`keywords/scf_settings` for input keywords:

- ``mixtype`` - Mixing type (``&mix``)
- ``beta`` - Mixing parameter (``&mix``)

Convergence Acceleration
========================

**When convergence stalls:**

- Reduce mixing parameter ``beta``
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

**Method:**

``bands%fermi`` integrates the DOS cumulatively with Simpson's rule up to the
valence count and interpolates linearly inside the last interval. It runs only
for ``calctype`` ``'B'`` or ``'L'`` with ``fix_fermi = .false.``; with
``fix_fermi = .true.`` the input Fermi level is kept (useful for impurities).

Temperature
===========

The real-space SCF and Fermi level are at zero temperature. ``&energy
temperature`` (K) only smears the energy integral of the exchange (:math:`J_{ij}`)
calculation; ``&reciprocal temperature`` (K) is used by the k-space route.

Energy Calculation
==================

The **total energy** at each SCF iteration includes:

.. math::

   E_{\text{tot}} = T[G] + E_{\text{XC}}[n] + E_{\text{Hartree}}[n] - \mu N_e + \text{const}

where :math:`T[G]` is kinetic energy from Green's function:

.. math::

   T = \int_{-\infty}^{E_F} dE \, E \, \rho(E)

No forces are computed.

Magnetic Convergence
====================

For **magnetic systems** (when the local density has unequal spin channels;
``control%nsp`` is a global mode, not a magnetization test):

- Separate SCF loops for spin-up and spin-down
- Or: fully magnetized iterative solution
- Convergence checked on **magnetization** as well as charge

The atomic/XC radial path uses two local spin-density channels even for the
scalar-relativistic collinear mode ``nsp=1``. In non-collinear modes, these
are local spin eigenchannels rather than globally collinear spin axes.

**Constraints:**

- **Fixed magnetization** - Constrain :math:`m_z` per atom (exchange-field approach)
- **Constrained moments** - For magnetic impurity problems

Implemented via exchange field: additional :math:`\pm B_{\text{exch}} \sigma_z` term in Hamiltonian.

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

1. **Too aggressive mixing** - Reduce ``beta``
2. **Inadequate recursion depth** - Increase ``lld``
3. **Energy mesh too coarse** - Increase ``channels_ldos``
4. **Unstable Hamiltonian** - Check element/potential parameters
5. **Bad initial guess** - Start from converged nearby configuration

**Debugging strategies:**

- Monitor energy evolution (should decrease monotonically)
- Plot charge density per iteration
- Check for charge oscillations
- Increase output verbosity (``verbose = .true.``)

SCF Log
=======

Each iteration logs ``Not converged! Diff= ...``; the last one logs
``Converged!`` followed by the final ``Diff``.

Provenance
==========

Main code locations:

- **SCF loop and convergence check:** ``source/self.f90::run()``
- **Charge integration:** ``source/density_of_states.f90::density()``
- **Density mixing:** ``source/mix.f90``
- **Potential update:** ``source/charge.f90::bulkpot()``

Key SCF parameters: ``source/self.f90`` (``nstep``, ``conv_thr``)

See Also
========

- :ref:`theory/green_functions` - Green's function calculation in SCF
- :ref:`theory/recursion_method` - Recursion method used in each iteration
- :ref:`keywords/scf_settings` - Input control parameters
- :ref:`user_guide/input_files` - SCF-related namelists
