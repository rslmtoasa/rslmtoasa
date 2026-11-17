.. _theory/spin_dynamics:

=======================================
Atomistic Spin Dynamics (ASD)
=======================================

Overview
========

The spin dynamics module implements **atomistic spin dynamics** simulations for studying
time-dependent magnetic phenomena in materials. This approach treats each atomic magnetic
moment as a classical spin vector that evolves according to the Landau-Lifshitz-Gilbert (LLG)
equation or related stochastic variants.

**Key applications:**

- Magnetic response to external fields
- Spin relaxation and damping
- Temperature-dependent magnetization dynamics
- Magnon excitations
- Magnetic domain wall motion

Landau-Lifshitz-Gilbert Equation
=================================

The LLG equation governs the time evolution of magnetic moments:

.. math::

   \frac{d\mathbf{m}_i}{dt} = -\gamma \mathbf{m}_i \times \mathbf{B}_i^{\text{eff}} 
   + \alpha \mathbf{m}_i \times \frac{d\mathbf{m}_i}{dt}

where:

- :math:`\mathbf{m}_i` is the unit magnetic moment vector at site :math:`i`
- :math:`\gamma` is the gyromagnetic ratio
- :math:`\mathbf{B}_i^{\text{eff}}` is the effective field (derivative of Hamiltonian)
- :math:`\alpha` is the Gilbert damping parameter

**Effective field:**

The field is not built from a :math:`J_{ij}` spin Hamiltonian. At each
predictor and corrector step the code runs the electronic SCF
(``self%run()``), computes the torque field ``bands%mag_for`` with
``bands%calculate_magnetic_torques()`` and converts it to Tesla
(``sd_run`` in ``source/spin_dynamics.f90``). The thermal field is generated
inside ``depondt.f90``.

Stochastic LLG (Langevin Dynamics)
===================================

For finite temperature simulations, thermal fluctuations are included via a stochastic field:

.. math::

   \frac{d\mathbf{m}_i}{dt} = -\gamma \mathbf{m}_i \times \mathbf{B}_i^{\text{eff}} 
   + \alpha \mathbf{m}_i \times \frac{d\mathbf{m}_i}{dt}
   + \gamma \mathbf{m}_i \times \mathbf{B}_i^{\text{th}}(t)

**Thermal field properties:**

.. math::

   \langle \mathbf{B}_i^{\text{th}}(t) \rangle = 0

   \langle B_i^{\text{th},\mu}(t) B_j^{\text{th},\nu}(t') \rangle 
   = 2D \delta_{ij} \delta_{\mu\nu} \delta(t-t')

with noise strength:

.. math::

   D = \frac{\alpha k_B T}{\gamma \mu_i}

Numerical Integration Methods
==============================

**Depondt-Mertens Integrator**

The only integrator is the Depondt-Mertens solver
(``source/include_codes/abspinlib/depondt.f90``): a Heun-type
predictor-corrector that rotates each moment about the effective field instead
of adding increments, so the moment magnitude is kept
(``depondt_evolve_first`` / ``depondt_evolve_second``).

Implementation in RS-LMTO
==========================

**Module structure (spin_dynamics.f90):**

Type: spin_dynamics
-------------------

**Key members** (read from ``&sd``):

- ``asd_step`` - Number of spin-dynamics steps
- ``dt`` - Time step (the code passes ``dt*1.0d-12`` to the integrator)
- ``alpha`` - Gilbert damping
- ``sd_temp`` - Temperature
- ``integrator``, ``nt``, ``t_i``, ``t_f`` - Read but not used

**Procedures:**

- ``sd_run()`` - Spin-dynamics loop
- ``asd_pred()``, ``asd_corr()`` - Predictor and corrector steps
- ``writelammpstraj()`` - Trajectory output
- ``build_from_file()``, ``restore_to_default()`` - Input and defaults

Typical Workflow
================

``sd_run`` repeats for ``asd_step`` steps:

1. **Field** - Electronic SCF, then torque field ``bands%mag_for``
2. **Predictor** - Depondt step with that field
3. **Field** - Electronic SCF at the predicted moments, new torque field
4. **Corrector** - Depondt step, then electronic SCF at the corrected moments

Input Parameters
================

In ``&sd`` (``source/include_codes/namelists/spin_dynamics.f90``):

- ``asd_step`` - Number of steps (default: 0)
- ``dt`` - Time step (default: ``1.0d-16``; see the note on units below)
- ``sd_temp`` - Temperature (default: 0)
- ``alpha`` - Gilbert damping :math:`\alpha` (default: 0)

``integrator``, ``nt``, ``t_i`` and ``t_f`` are read but ignored. The run is
selected with ``processing = 'sd'`` in ``&calculation``. The unit of ``dt``
is not settled: the source comment says fs, while the code multiplies ``dt``
by ``1.0d-12`` (``ps_to_s``) before calling the integrator.

Example Input
=============

.. code-block:: fortran

   &calculation
      processing = 'sd'
   /

   &sd
      asd_step = 100
      sd_temp = 300.0
      alpha = 0.1
   /

Output Files
============

**output.lammpstrj**

LAMMPS dump written after each step: timestep, number of atoms, and per atom
the atomic number, position and unit spin vector (``ITEM: ATOMS type x y z vx
vy vz``). The usual SCF report is also written.

Provenance
==========

**Implementation:**

- ``source/spin_dynamics.f90`` - Main ASD module
- ``source/exchange.f90`` - Exchange parameter calculation
- ``source/include_codes/abspinlib/depondt.f90`` - Depondt-Mertens integrator
- ``source/include_codes/abspinlib/randomnumbers.f90`` - Thermal noise generation

**Theory references:**

[1] Landau, L. D., & Lifshitz, E. (1935). Theory of magnetic domains. *Phys. Z. Sowjetunion*, 8, 153.

[2] Gilbert, T. L. (2004). A phenomenological theory of damping in ferromagnetic materials. *IEEE Trans. Magn.*, 40(6), 3443-3449.

[3] Mentink, J. H., et al. (2010). Stable and fast semi-implicit integration of the stochastic Landau-Lifshitz equation. *J. Phys.: Condens. Matter*, 22(17), 176001.

[4] Evans, R. F., et al. (2014). Atomistic spin model simulations of magnetic nanomaterials. *J. Phys.: Condens. Matter*, 26(10), 103202.

See Also
========

- :doc:`../user_guide/examples` - Spin dynamics examples
- :ref:`keywords/control_parameters` - Control parameters
- :doc:`green_functions` - Green's function formalism
- ``source/exchange.f90`` - Exchange calculation details
