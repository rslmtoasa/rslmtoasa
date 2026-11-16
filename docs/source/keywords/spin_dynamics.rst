.. _keywords/spin_dynamics:

=======================================
Spin Dynamics Parameters (&sd namelist)
=======================================

Overview
========

The ``&sd`` namelist controls the **Atomistic Spin Dynamics (ASD)** simulation
that is launched after a self-consistent calculation when ``processing = 'sd'``
is set in ``&calculation``.  The module integrates the Landau-Lifshitz-Gilbert (LLG)
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

``t_i``, ``t_f``, ``nt``, ``integrator`` and ``b_stochastic`` are read from ``&sd`` but not used by
the current code (``source/spin_dynamics.f90``).

t_i
---

**Type:** Real

**Units:** ps (picoseconds)

**Default:** ``0.0``

Read but currently not used.

t_f
---

**Type:** Real

**Units:** ps

**Default:** ``0.0``

Read but currently not used.

dt
--

**Type:** Real

**Units:** ps

**Default:** ``1.0e-16``

**Purpose:** Integration time step.

**Typical range:** ``1.0e-4`` to ``1.0e-2`` ps

**Example:**

.. code-block:: fortran

   dt = 1.0e-3   ! 1 fs step

**Notes:**

- Smaller steps give higher accuracy but increase cost in proportion.
- The Depondt integrator (the only one in the code) is symplectic and
  typically allows larger steps than a plain Euler scheme.
- A convergence check: halve ``dt`` and verify that observables are unchanged.

nt
--

**Type:** Integer

Read but currently not used. The number of steps is ``asd_step``.

alpha
-----

**Type:** Real

**Default:** ``0.0``

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

**Type:** Real array, shape ``(3, nrec)``

Read but currently not used. It is allocated and set to zero only. For thermal ASD use
``sd_temp``.

integrator
----------

**Type:** Character string (len=5)

**Default:** ``'none'``

Read but currently not used. The spin-dynamics step always calls the Depondt-Mertens solver
[#depondt]_ (``source/spin_dynamics.f90``); no other integrator exists in the code.

.. [#depondt] P. Depondt and F. G. Mertens, *J. Phys.: Condens. Matter* **21**,
   336005 (2009).  DOI: `10.1088/0953-8984/21/33/336005 <https://doi.org/10.1088/0953-8984/21/33/336005>`_

asd_step
--------

**Type:** Integer

**Default:** ``0``

**Purpose:** Number of ASD steps: the loop of ``source/spin_dynamics.f90`` runs
``do n_step = 1, asd_step``. With the default, no step is performed.

Minimal Working Example
=======================

A 5 ps simulation of bcc Fe at room temperature with Depondt integrator:

.. code-block:: fortran

   &calculation
     processing = 'sd'
   /

   &sd
     dt         = 5.0e-4
     asd_step   = 10000   ! 10000 steps of 5.0e-4 ps = 5 ps
     alpha      = 0.01
     sd_temp    = 300.0
   /

Provenance
==========

- **Namelist definition:** ``source/include_codes/namelists/spin_dynamics.f90``
- **Module implementation:** ``source/spin_dynamics.f90``
- **Trigger:** ``processing = 'sd'`` in ``source/calculation.f90``
- **Exchange coupling input:** ``source/exchange.f90``

See Also
========

- :ref:`theory/spin_dynamics` — LLG equation, effective field, integrator derivations
- :ref:`keywords/control_parameters` — ``do_asd``, ``asd_jij`` (read but not used)
- :ref:`reference/exchange_module` — Computing :math:`J_{ij}` for ASD input
