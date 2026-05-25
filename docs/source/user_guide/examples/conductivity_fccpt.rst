.. _examples/conductivity_fccpt:

===================================================
Post-Processing: Spin Hall Effect in Pt
===================================================

**Location:** ``example/conductivity/fccPt/``

**System:** Spin Hall conductivity of fcc Pt

**Physics:** Large spin Hall effect due to strong spin-orbit coupling in 5d element.

Overview
========

Platinum exhibits a large spin Hall effect - converting charge current to spin current via spin-orbit coupling. This example calculates the spin Hall conductivity using **fully relativistic** treatment (``nsp = 4``).

Files
=====

- ``input.nml`` - Conductivity with nsp=4
- ``Pt_scf.nml`` - Pt parameters with spin-orbit coupling

Input File: input.nml
=====================

.. code-block:: fortran

   &calculation
   pre_processing = 'bravais'
   verbose = T
   /
   
   &lattice
   rc = 60
   alat = 3.923                   ! FCC Pt lattice constant
   crystal_sym = 'fcc'
   wav = 1.4844
   ntype = 1
   ct(1) = 4.0
   r2 = 16.0
   /
   
   &atoms
   database = './'
   label(1) = 'Pt'
   /
   
   &self
   nstep = 100
   /
   
   &hamiltonian
   v_alpha = 0, 1, 0              ! Charge current (y)
   v_beta = 1, 0, 0               ! Electric field (x)
   js_alpha = 'z'                 ! Spin current (z)
   /
   
   &energy
   fermi = -0.095
   energy_min = -1.0
   energy_max = 1.5
   channels_ldos = 2500
   /
   
   &control
   calctype = 'B'
   nsp = 4                        ! ← Fully relativistic (strong SOC!)
   cond_ll = 500
   recur = 'chebyshev'
   lld = 150
   cond_type = 'orb'              ! ← Orbital contribution
   /
   
   &mix
   beta = 0.3
   mixtype = 'broyden'
   /

Key Difference: nsp = 4
-----------------------

``nsp = 4`` enables **fully relativistic** calculation:

- 4-component spinors
- Full spin-orbit coupling
- Essential for heavy elements (Z > 70)
- Required for spin Hall effect

Spin Hall Effect
================

The spin Hall conductivity relates charge current to transverse spin current:

.. math::

   \\mathbf{j}_s^z = \\sigma_{SH} \\frac{\\hbar}{2e} \\mathbf{E} \\times \\hat{z}

where the spin Hall angle is:

.. math::

   \\theta_{SH} = \\frac{\\sigma_{SH}}{\\sigma_{charge}}

Running the Calculation
=======================

.. code-block:: bash

   cd example/conductivity/fccPt
   ../../../../build/bin/rslmto.x

Expected Output
===============

**Spin Hall angle:**

For Pt at room temperature:

.. code-block:: text

   θ_SH ≈ 0.05-0.10 (5-10%)
   
   Much larger than light metals:
   - Cu: ~0.001 (0.1%)
   - Al: ~0.0001 (0.01%)

**Expected runtime:** 20-30 minutes (nsp=4 is computationally expensive).

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

- :doc:`bulk_mn3sn` - Another nsp=4 example
- :doc:`conductivity_bccfe` - Anomalous Hall effect
- :doc:`conductivity_fcccu` - Non-relativistic conductivity
- :doc:`../../keywords/control_parameters` - nsp parameter
