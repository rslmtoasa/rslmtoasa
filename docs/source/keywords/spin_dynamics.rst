.. _keywords/spin_dynamics:

======================================
Spin Dynamics Parameters (&sd namelist)
======================================

Overview
========

The ``&sd`` namelist controls the **Atomistic Spin Dynamics (ASD)** simulation
that is launched after a self-consistent calculation when ``do_asd = .true.`` is
set in ``&control``.  The module integrates the Landau-Lifshitz-Gilbert (LLG)
equation for the local magnetic moment vectors :math:`\{\mathbf{m}_i(t)\}`:

.. math::

   \frac{d\mathbf{m}_i}{dt} =
     -\gamma \mathbf{m}_i \times \mathbf{B}_i^{\text{eff}}
     + \alpha \, \mathbf{m}_i \times \frac{d\mathbf{m}_i}{dt}
     + \boldsymbol{\xi}_i(t),

where :math:`\mathbf{B}_i^{\text{eff}}` is the effective field (from the exchange
Hamiltonian), :math:`\alpha` is the Gilbert damping, and
:math:`\boldsymbol{\xi}_i(t)` is a stochastic Langevin noise field.

For the underlying LLG theory see :ref:`theory/spin_dynamics`.

Parameters
==========

t_i
---

**Type:** Real

**Units:** ps (picoseconds)

**Default:** ``0.0``

**Purpose:** Initial simulation time.  Normally zero for a fresh start; can be
set to a non-zero value when restarting from a checkpoint.

**Example:**

.. code-block:: fortran

   t_i = 0.0

t_f
---

**Type:** Real

**Units:** ps

**Purpose:** Final simulation time.  The simulation runs from ``t_i`` to
``t_f``.

**Example:**

.. code-block:: fortran

   t_f = 10.0   ! 10 ps simulation

dt
--

**Type:** Real

**Units:** ps

**Purpose:** Integration time step.

**Typical range:** ``1.0e-4`` to ``1.0e-2`` ps

**Example:**

.. code-block:: fortran

   dt = 1.0e-3   ! 1 fs step

**Notes:**

- Smaller steps give higher accuracy but increase cost in proportion.
- The Depondt integrator (``integrator = 'depon'``) is symplectic and
  typically allows larger steps than a plain Euler scheme.
- A convergence check: halve ``dt`` and verify that observables are unchanged.

nt
--

**Type:** Integer

**Purpose:** Number of time steps to perform.  Alternatively the simulation
length can be specified via ``t_f`` and ``dt``; if all three are set,
``nt`` takes precedence.

**Example:**

.. code-block:: fortran

   nt = 10000

alpha
-----

**Type:** Real

**Default:** ``0.01``

**Purpose:** Gilbert damping parameter :math:`\alpha` in the LLG equation.
Controls the rate at which energy is dissipated from the spin system.

**Example:**

.. code-block:: fortran

   alpha = 0.05   ! Moderate damping

**Notes:**

- :math:`\alpha = 0`: undamped precession (energy-conserving).
- Small :math:`\alpha` (0.001–0.01): weakly damped, realistic for many metals.
- Large :math:`\alpha` (0.1–1.0): overdamped relaxation toward equilibrium.
- Experimentally measured values for common ferromagnets: Fe ≈ 0.002,
  Co ≈ 0.005, Ni ≈ 0.025.

sd_temp
-------

**Type:** Real

**Units:** Kelvin

**Default:** ``0.0``

**Purpose:** Simulation temperature for the stochastic Langevin noise term
:math:`\boldsymbol{\xi}_i(t)`.  At ``sd_temp = 0`` the noise is absent and the
LLG equation is deterministic.

**Example:**

.. code-block:: fortran

   sd_temp = 300.0   ! Room temperature

**Notes:**

- The noise amplitude is related to the temperature via the
  fluctuation-dissipation theorem:

  .. math::

     \langle \xi_i^\mu(t)\, \xi_j^\nu(t') \rangle =
       \frac{2 \alpha k_B T}{\gamma |\mathbf{m}_i|}\,
       \delta_{ij}\, \delta^{\mu\nu}\, \delta(t - t').

- Finite-temperature ASD is useful for simulating magnon spectra, Curie
  temperature estimates, and thermally activated switching.

b_stochastic
------------

**Type:** Real array, shape ``(3, natoms)``

**Purpose:** Optional fixed stochastic magnetic field added to each atom's
effective field, supplementing the thermal noise.  Primarily used for testing
or applying a spatially varying external perturbation.

**Notes:**

- For standard thermal ASD use ``sd_temp`` instead.

integrator
----------

**Type:** Character string (len=5)

**Default:** ``'depon'``

**Purpose:** Selects the numerical integration scheme for the LLG equation

**Allowed values:**

.. list-table::
   :widths: 15 70
   :header-rows: 1

   * - Value
     - Description
   * - ``'depon'``
     - **Depondt-Mertens** symplectic integrator [#depondt]_.
       Second-order accurate, nearly energy-conserving for small :math:`\alpha`.
       Recommended for most simulations.
   * - ``'heun'``
     - **Heun** (second-order Runge-Kutta) predictor-corrector.  More accurate
       per step but roughly twice the cost.  Useful for stochastic
       (finite-temperature) simulations where the Stratonovich convention is
       required.

**Example:**

.. code-block:: fortran

   integrator = 'depon'

.. [#depondt] P. Depondt and F. G. Mertens, *J. Phys.: Condens. Matter* **21**,
   336005 (2009).  DOI: `10.1088/0953-8984/21/33/336005 <https://doi.org/10.1088/0953-8984/21/33/336005>`_

asd_step
--------

**Type:** Integer

**Purpose:** Internal step counter; written to checkpoint files and read on
restart.  Normally managed automatically by the code — do not set manually
unless continuing a previously interrupted run.

Minimal Working Example
=======================

A 5 ps simulation of bcc Fe at room temperature with Depondt integrator:

.. code-block:: fortran

   &control
     do_asd  = .true.
     asd_jij = .true.   ! use computed J_ij as effective field
   /

   &sd
     t_i        = 0.0
     t_f        = 5.0
     dt         = 5.0e-4
     alpha      = 0.01
     sd_temp    = 300.0
     integrator = 'depon'
   /

Provenance
==========

- **Namelist definition:** ``source/include_codes/namelists/spin_dynamics.f90``
- **Module implementation:** ``source/spin_dynamics.f90``
- **Control flags:** ``source/control.f90`` — ``do_asd``, ``asd_jij``
- **Exchange coupling input:** ``source/exchange.f90``

See Also
========

- :ref:`theory/spin_dynamics` — LLG equation, effective field, integrator derivations
- :ref:`keywords/control_parameters` — ``do_asd``, ``asd_jij`` flags
- :ref:`reference/exchange_module` — Computing :math:`J_{ij}` for ASD input
