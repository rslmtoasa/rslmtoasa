.. _examples/conductivity_fccpt:

===================================================
Post-Processing: Spin Hall Effect in Pt
===================================================

**Location:** ``example/conductivity/fccPt/``

**System:** Spin Hall conductivity of fcc Pt

**Physics:** Large spin Hall effect due to strong spin-orbit coupling in 5d element.

Overview
========

Platinum exhibits a large spin Hall effect - converting charge current to spin current via spin-orbit coupling. This example calculates the spin Hall conductivity using a collinear calculation with spin-orbit coupling (``nsp = 2``).

Files
=====

- ``input.nml`` - Conductivity with nsp=2
- ``Pt.nml`` - Pt parameters

Input File: input.nml
=====================

The SCF input in ``example/conductivity/fccPt/input.nml``:

.. code-block:: fortran

   &calculation
   pre_processing = 'bravais'
   verbose = T
   /
   &lattice
   rc = 060
   alat = 3.9264951876
   crystal_sym = 'fcc'
   wav = 1.5344598728
   pbc = .true.
   b1 = .false.
   b2 = .false.
   b3 = .false.
   n1 = 20
   n2 = 20
   n3 = 20
   ntype = 1
   ct(1) = 4.0d0
   r2 = 16.0d0
   /
   &atoms
   database = './'
   label(1) = 'Pt'
   /
   &self
   ws_all = .true.
   nstep = 100
   /
   &hamiltonian
   v_alpha = 0, 1, 0              ! Charge current (y)
   v_beta =  1, 0, 0              ! Electric field (x)
   pol_alpha = 'z'                ! Spin current (z)
   pol_beta = 'z'
   /
   &energy
   fermi = -0.089509
   energy_min = -1.0
   energy_max = 1.2
   channels_ldos = 2500
   /
   &control
   calctype = 'B'
   nsp = 2                        ! Collinear, with spin-orbit coupling
   cond_ll = 500
   lld = 150
   recur = 'chebyshev'
   linear_in = 'charge'
   linear_out = 'spin'            ! Spin current response
   cond_calctype = 'per_type'
   /
   &mix
   mixtype = 'broyden'
   beta = 0.4
   /

Spin-Orbit Coupling: nsp = 2
----------------------------

``nsp = 2`` is the collinear calculation with spin-orbit coupling. ``nsp = 4`` (non-collinear with
spin-orbit coupling) is not used in this example.

Spin Hall Effect
================

The spin Hall conductivity relates charge current to transverse spin current:

.. math::

   \mathbf{j}_s^z = \sigma_{SH} \frac{\hbar}{2e} \mathbf{E} \times \hat{z}

where the spin Hall angle is:

.. math::

   \theta_{SH} = \frac{\sigma_{SH}}{\sigma_{charge}}

Running the Calculation
=======================

The calculation has two steps. With ``post_processing = 'none'`` (the default) the file above runs
the SCF only:

.. code-block:: bash

   mkdir -p runs
   cp -r example/conductivity/fccPt runs/fccPt
   cd runs/fccPt
   ../../build/bin/rslmto.x

The conductivity step does not repeat the SCF. The folders ``sigma_xy_z-spin`` (``linear_out = 'spin'``)
and ``sigma_xy_z-orb`` (``linear_out = 'orbital'``) hold an ``input.nml`` that differs from the file
above only in ``post_processing = 'conductivity'`` (in place of ``pre_processing = 'bravais'``) and, for
the orbital case, ``linear_out``. Copy the converged ``Pt_out.nml`` (read before ``Pt.nml``) into the
subfolder and run there:

.. code-block:: bash

   cp Pt_out.nml sigma_xy_z-spin/
   cd sigma_xy_z-spin
   ../../../build/bin/rslmto.x

Expected Output
===============

**Spin Hall angle:**

For Pt at room temperature:

.. code-block:: text

   θ_SH ≈ 0.05-0.10 (5-10%)
   
   Much larger than light metals:
   - Cu: ~0.001 (0.1%)
   - Al: ~0.0001 (0.01%)

**Expected runtime:** 20-30 minutes.

Physical Interpretation
========================

**Why Pt has large spin Hall effect:**

1. Strong spin-orbit coupling in 5d orbitals (scales as Z⁴)
2. Large DOS at Fermi level
3. Favorable band structure near E_F

**Applications:**

- **Spin current source** in spintronics
- **Spin-orbit torque** devices
- **Magnetic memory and logic**

Pt is the most widely used spin Hall material in experimental devices.

See Also
========

- :doc:`bulk_mn3sn` - Non-collinear (nsp=4) example
- :doc:`conductivity_bccfe` - Anomalous Hall effect
- :doc:`conductivity_fcccu` - Conductivity of Cu
- :doc:`../../keywords/control_parameters` - nsp parameter
