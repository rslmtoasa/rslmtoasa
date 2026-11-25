.. _examples/conductivity_fcccu:

===================================================
Post-Processing: Conductivity of Cu
===================================================

**Location:** ``example/conductivity/fccCu/``

**System:** Electrical conductivity tensor for fcc Cu

**Physics:** Calculate conductivity using Kubo-Bastin formula with Chebyshev expansion (KPM).

Overview
========

This example calculates the electrical conductivity tensor σ_αβ(E) as a function of energy using:

- Kubo-Bastin formula for disordered systems
- Chebyshev polynomial expansion (Kernel Polynomial Method)
- Velocity operator matrix elements

Files
=====

- ``input.nml`` - Conductivity calculation setup
- ``Cu.nml`` - Cu tight-binding parameters

Input File: input.nml
=====================

.. code-block:: fortran

   &calculation
   pre_processing = 'bravais'
   verbose = T
   /
   
   &lattice
   rc = 60
   alat = 3.61411                 ! FCC Cu lattice constant
   crystal_sym = 'fcc'
   wav = 1.41237
   pbc = .true.
   b1 = .true.                    ! ← Periodic cluster
   b2 = .true.
   b3 = .true.
   n1 = 20
   n2 = 20
   n3 = 20
   ntype = 1
   ct(1) = 4.0
   r2 = 16.0
   /
   
   &atoms
   database = './'
   label(1) = 'Cu'
   /
   
   &self
   ws_all = .true.
   nstep = 100
   /
   
   &hamiltonian
   v_alpha = 0, 1, 0              ! ← Current direction (y)
   v_beta = 1, 0, 0               ! ← Field direction (x)
   pol_alpha = 'z'                 ! ← Spin current (z)
   pol_beta = 'z'
   /
   
   &energy
   fermi = -0.088977
   energy_min = -1.0
   energy_max = 1.2
   channels_ldos = 2500
   /
   
   &control
   calctype = 'B'
   nsp = 2                        ! Collinear, with spin-orbit coupling
   cond_ll = 500                  ! ← Chebyshev expansion order
   recur = 'chebyshev'            ! ← Use Chebyshev recursion
   lld = 150
   linear_in = 'charge'
   linear_out = 'spin'             ! Spin contribution
   cond_calctype = 'per_type'     ! ← Per-atom-type resolution
   /
   
   &mix
   mixtype = 'broyden'
   beta = 0.4
   /

Critical Conductivity Parameters
---------------------------------

- ``cond_ll = 500``: Chebyshev expansion order (higher -> better resolution, but slower)
- ``linear_out = 'spin'`` or ``'orbital'``: Calculate spin or orbital contribution (``linear_in = 'charge'``)
- ``cond_calctype = 'per_type'``: Atom-type resolved conductivity
- ``v_alpha, v_beta``: Velocity operator directions for tensor component σ_αβ
- ``b1 = b2 = b3 = .true.``: Periodic cluster; the linear-response calculation does not need open boundaries, and periodic boundaries give better numerics
- ``recur = 'chebyshev'``: Use Chebyshev recursion (KPM method)

Velocity Operator: Selecting Tensor Component
==============================================

The ``&hamiltonian`` namelist specifies which conductivity tensor element to calculate:

.. code-block:: text

   v_alpha = 0, 1, 0    ! y-direction
   v_beta = 1, 0, 0     ! x-direction
   => Calculates sigma_yx

**Common configurations:**

- **Longitudinal:** ``v_alpha = (1,0,0), v_beta = (1,0,0)`` -> σ_xx
- **Hall:** ``v_alpha = (0,1,0), v_beta = (0,0,1)`` -> σ_yz
- **Cross terms:** Any combination

For cubic Cu, by symmetry: σ_xx = σ_yy = σ_zz and off-diagonal terms are zero.

Running the Calculation
=======================

The calculation has two steps. With ``post_processing = 'none'`` (the default) the file above runs
the SCF only:

.. code-block:: bash

   mkdir -p runs
   cp -r example/conductivity/fccCu runs/fccCu
   cd runs/fccCu
   ../../build/bin/rslmto.x

The conductivity step does not repeat the SCF. The folders ``sigma_xy_z-spin`` (``linear_out = 'spin'``)
and ``sigma_xy_z-orb`` (``linear_out = 'orbital'``) hold an ``input.nml`` that differs from the file
above only in ``post_processing = 'conductivity'`` (in place of ``pre_processing = 'bravais'``) and, for
the orbital case, ``linear_out``. Copy the converged ``Cu_out.nml`` (read before ``Cu.nml``) into the
subfolder and run there:

.. code-block:: bash

   cp Cu_out.nml sigma_xy_z-spin/
   cd sigma_xy_z-spin
   ../../../build/bin/rslmto.x

Expected Output
===============

**Output files:**

- ``cond_total.out`` - Total response vs energy (raw kernel, no header)
- ``Cu_cond.out`` - Cu contribution (if per_type)
- ``cond_total_orb_real.out``, ``cond_total_orb_im.out``, ``Cu_cond_orb_real.out``, ``Cu_cond_orb_im.out`` - Orbital-resolved parts
- ``response_total.out``, ``Cu_response.out``, ``response_units.txt`` - Response in physical units and the units used

**Expected results:**

The shipped deck (``v_alpha = 0,1,0``, ``v_beta = 1,0,0``, ``pol = 'z'``, ``linear_out = 'spin'``)
computes the :math:`\sigma_{yx}` spin-current response, not the charge conductivity. The value
below refers to the charge conductivity :math:`\sigma_{xx}`, which needs ``linear_out = 'charge'`` and
``v_alpha = v_beta = 1, 0, 0``. For fcc Cu at Fermi level:

.. code-block:: text

   σ_xx(E_F) ≈ 6.0 × 10⁷ (Ω·m)⁻¹
   
   Experimental: 5.96 × 10⁷ (Ω·m)⁻¹
   Agreement: ~1%!

**View results:**

.. code-block:: bash

   head cond_total.out
   
   # Plot the real part vs energy
   gnuplot -e "plot 'cond_total.out' u 1:2 w l; pause -1"

**Output format:**

.. code-block:: text

   # cond_total.out (no header; energies relative to E_F)
   E-E_F (Ry)      Re              Im

   # response_total.out
   # E-Ef [eV], Re, Im; units: ...

**Expected runtime:** 10-20 minutes.

Physical Interpretation
========================

**Conductivity peak at Fermi level:**

- High DOS at E_F (d-band contribution)
- Long scattering time (weak disorder)
- Cu is excellent conductor

**Energy dependence:**

- Below E_F: Decreasing σ (fewer available states)
- Above E_F: Also decreasing (fewer occupied states)
- Maximum at E_F for metals

Troubleshooting
===============

**Poor energy resolution:**

- Increase ``cond_ll`` to 1000 or higher
- Increase ``lld`` to 200-300

**Wrong conductivity value:**

- Check Fermi energy is correct
- Verify velocity operators are normalized
- Ensure ``pbc`` and ``b1 = b2 = b3 = .true.`` (periodic cluster)

See Also
========

- :doc:`conductivity_fccpt` - Spin Hall effect in Pt
- :doc:`conductivity_bccfe` - Anomalous Hall effect in Fe
- :doc:`../../reference/conductivity_module` - Conductivity calculation details
- :doc:`../../keywords/control_parameters` - cond_ll, linear_in, linear_out parameters
