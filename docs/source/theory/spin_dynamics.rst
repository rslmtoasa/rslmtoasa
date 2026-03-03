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

.. math::

   \mathbf{B}_i^{\text{eff}} = -\frac{1}{\mu_i} \frac{\partial H}{\partial \mathbf{m}_i}

Includes contributions from:

- Exchange interactions: :math:`J_{ij} \mathbf{m}_i \cdot \mathbf{m}_j`
- External field: :math:`\mathbf{B}_{\text{ext}}`
- Anisotropy: :math:`K_i (\mathbf{m}_i \cdot \hat{\mathbf{e}})^2`
- Thermal fluctuations: :math:`\mathbf{B}_i^{\text{th}}`

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

A predictor-corrector method preserving :math:`|\mathbf{m}_i| = 1`:

1. **Predictor step:**

   .. math::

      \mathbf{m}_i^{n+1/2} = \mathbf{m}_i^n + \frac{\Delta t}{2} \frac{d\mathbf{m}_i^n}{dt}

2. **Normalize:**

   .. math::

      \tilde{\mathbf{m}}_i^{n+1/2} = \frac{\mathbf{m}_i^{n+1/2}}{|\mathbf{m}_i^{n+1/2}|}

3. **Corrector step:**

   .. math::

      \mathbf{m}_i^{n+1} = \mathbf{m}_i^n + \Delta t \frac{d\tilde{\mathbf{m}}_i^{n+1/2}}{dt}

4. **Final normalization**

**Heun Method**

Second-order Runge-Kutta for better accuracy:

.. math::

   \mathbf{k}_1 = f(\mathbf{m}_i^n)

   \mathbf{k}_2 = f(\mathbf{m}_i^n + \Delta t \mathbf{k}_1)

   \mathbf{m}_i^{n+1} = \mathbf{m}_i^n + \frac{\Delta t}{2}(\mathbf{k}_1 + \mathbf{k}_2)

followed by normalization.

Implementation in RS-LMTO
==========================

**Module structure (spin_dynamics.f90):**

Type: spin_dynamics
-------------------

**Key members:**

- ``integrator`` - Integration scheme ('depon' or 'heun')
- ``nt`` - Number of time steps
- ``dt`` - Time step (femtoseconds)
- ``t_i, t_f`` - Initial and final times
- ``asd_step`` - Current ASD step number
- Temperature - Simulation temperature

**Main procedures:**

- ``process()`` - Run full spin dynamics simulation
- ``iterate()`` - Single time step integration
- ``calculate_effective_field()`` - Compute :math:`\mathbf{B}_i^{\text{eff}}`
- ``langevin_field()`` - Generate thermal noise
- ``update_moments()`` - Integrate LLG equation

Exchange Coupling from DFT
===========================

The exchange parameters :math:`J_{ij}` are calculated from DFT using the magnetic
force theorem and Green's function approach:

.. math::

   J_{ij} = -\frac{1}{4\pi} \int_{-\infty}^{E_F} dE \, 
   \text{Im} \text{Tr}[t_i G_{ij}^{\uparrow}(E) t_j G_{ji}^{\downarrow}(E)]

where:

- :math:`G_{ij}^{\sigma}(E)` is the inter-site Green's function for spin :math:`\sigma`
- :math:`t_i` is the scattering matrix at site :math:`i`

See ``exchange.f90`` for implementation.

Observables
===========

**Magnetization:**

.. math::

   \mathbf{M}(t) = \frac{1}{N} \sum_i \mu_i \mathbf{m}_i(t)

**Energy:**

.. math::

   E(t) = -\sum_{ij} J_{ij} \mathbf{m}_i(t) \cdot \mathbf{m}_j(t)

**Autocorrelation function:**

.. math::

   C(t) = \langle \mathbf{M}(0) \cdot \mathbf{M}(t) \rangle

**Magnetic susceptibility:**

Fourier transform of correlation function relates to susceptibility:

.. math::

   \chi(\omega) = \int_0^\infty dt \, e^{i\omega t} C(t)

Typical Workflow
================

1. **SCF calculation** - Obtain ground state electronic structure
2. **Exchange calculation** - Compute :math:`J_{ij}` parameters
3. **Initial configuration** - Set initial magnetic moments
4. **Time evolution** - Integrate LLG for many time steps
5. **Analysis** - Extract magnetization dynamics, correlations

Input Parameters
================

In ``&control`` namelist:

- ``asd_nsteps`` - Number of ASD time steps (default: 1000)
- ``asd_dt`` - Time step in femtoseconds (default: 0.1 fs)
- ``asd_temp`` - Temperature in Kelvin (default: 0)
- ``asd_damping`` - Gilbert damping :math:`\alpha` (default: 0.1)
- ``asd_integrator`` - Integration method ('depondt' or 'heun')

Example Input
=============

.. code-block:: fortran

   &control
      nsp = 2                  ! Spin-polarized
      processing = 'sd'        ! Spin dynamics
      asd_nsteps = 10000
      asd_dt = 0.1             ! 0.1 fs time step
      asd_temp = 300.0         ! 300 K
      asd_damping = 0.1        ! Gilbert damping
      asd_integrator = 'depondt'
   /

Output Files
============

**asd_moments.dat**

Magnetic moments at each time step:

.. code-block:: text

   # Time(fs)  mx_1   my_1   mz_1   mx_2   my_2   mz_2  ...
   0.0         0.0    0.0    1.0    0.0    0.0    -1.0
   0.1         0.01   0.0    0.999  -0.01  0.0    -0.999
   ...

**asd_magnetization.dat**

Total magnetization vs time:

.. code-block:: text

   # Time(fs)  Mx     My     Mz     |M|
   0.0         0.0    0.0    1.0    1.0
   0.1         0.0    0.0    0.998  0.998
   ...

**asd_energy.dat**

System energy vs time:

.. code-block:: text

   # Time(fs)  E_exchange  E_anisotropy  E_total
   0.0         -10.5       0.0           -10.5
   0.1         -10.4       0.01          -10.39
   ...

Performance Considerations
==========================

**Time step selection:**

- Too large: Numerical instability, moments not conserved
- Too small: Computationally expensive
- Typical: 0.01-1 fs depending on exchange strength

**Scaling:**

- Memory: :math:`O(N)` for :math:`N` atoms
- Compute: :math:`O(N^2)` for full exchange matrix
- Can use cutoff radius for :math:`O(N)` scaling

**Parallelization:**

- OpenMP parallelization over atomic sites
- MPI for large systems (domain decomposition)

Provenance
==========

**Implementation:**

- ``source/spin_dynamics.f90`` - Main ASD module
- ``source/exchange.f90`` - Exchange parameter calculation
- ``source/Depondt.f90`` - Depondt-Mertens integrator (if present)
- ``source/RandomNumbers.f90`` - Thermal noise generation (if present)

**Theory references:**

[1] Landau, L. D., & Lifshitz, E. (1935). Theory of magnetic domains. *Phys. Z. Sowjetunion*, 8, 153.

[2] Gilbert, T. L. (2004). A phenomenological theory of damping in ferromagnetic materials. *IEEE Trans. Magn.*, 40(6), 3443-3449.

[3] Mentink, J. H., et al. (2010). Stable and fast semi-implicit integration of the stochastic Landau-Lifshitz equation. *J. Phys.: Condens. Matter*, 22(17), 176001.

[4] Evans, R. F., et al. (2014). Atomistic spin model simulations of magnetic nanomaterials. *J. Phys.: Condens. Matter*, 26(10), 103202.

See Also
========

- :doc:`../user_guide/examples` - Spin dynamics examples
- :ref:`keywords/control_parameters` - ASD input parameters
- :doc:`green_functions` - Green's function formalism
- ``source/exchange.f90`` - Exchange calculation details
