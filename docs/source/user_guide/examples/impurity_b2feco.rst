.. _examples/impurity_b2feco:

================================================
Impurity: Fe in B2 FeCo
================================================

**Location:** ``example/impurity/B2FeCo/``

**System:** Single Fe impurity replacing Co in ordered B2 FeCo

**Physics:** Demonstrates embedded cluster impurity calculation with local magnetic perturbation.

Overview
========

This example calculates a single Fe impurity in B2-ordered FeCo alloy using:

- ``pre_processing = 'newclubulk'`` - Build cluster with impurity
- ``calctype = 'I'`` - Impurity calculation
- ``nclu = 1`` - One impurity site
- ``inclu`` - Impurity position specification

Files
=====

- ``input.nml`` - Main parameters
- ``Fe_scf.nml`` - Bulk Fe site parameters
- ``Co_scf.nml`` - Bulk Co site parameters
- ``Fe-imp_scf.nml`` - Fe impurity parameters (different from bulk Fe!)

Input File: input.nml
=====================

.. code-block:: fortran

   &calculation
   pre_processing = 'newclubulk'  ! ← Build cluster with impurity
   verbose = T
   /
   
   &lattice
   rc = 60
   alat = 2.859                   ! B2 lattice constant
   crystal_sym = 'b2'             ! B2 (CsCl-type) structure
   wav = 1.406216
   ntype = 3                      ! Fe, Co, Fe-imp
   nclu = 1                       ! ← One impurity
   ct(:) = 3*3.0d0
   r2 = 9.0d0
   inclu(1, :) = 0.0, 0.0, 0.0    ! ← Impurity at origin
   /
   
   &atoms
   database = './'
   label(1) = 'Fe'                ! Bulk Fe atoms
   label(2) = 'Co'                ! Bulk Co atoms
   label(3) = 'Fe-imp'            ! Fe impurity (replaces Co)
   /
   
   &self
   nstep = 100
   /
   
   &energy
   fermi = -0.041385
   energy_min = -2.5              ! Wider window for impurity states
   energy_max = 1.5
   channels_ldos = 2500
   /
   
   &control
   calctype = 'I'                 ! ← Impurity calculation
   nsp = 2
   lld = 21
   recur = 'block'
   /
   
   &mix
   beta = 0.10
   mixtype = 'broyden'
   /

Critical Impurity Parameters
-----------------------------

- ``pre_processing = 'newclubulk'``: Builds finite cluster embedding impurity in infinite bulk
- ``calctype = 'I'``: **Impurity** calculation (embedded cluster method)
- ``nclu = 1``: Number of impurity sites (can have multiple impurities)
- ``inclu(1, :) = 0.0, 0.0, 0.0``: Position of impurity in fractional coordinates
- Wider ``energy_min/max``: Impurity may induce resonances outside bulk band range

B2 Structure: Auto-Generated
=============================

``crystal_sym = 'b2'`` generates CsCl-type structure:

.. code-block:: text

   Simple cubic with 2-atom basis:
   - Fe at (0, 0, 0)
   - Co at (0.5, 0.5, 0.5)
   
   Impurity setup:
   - Replace Co at origin with Fe-imp
   - Surrounded by bulk Fe and Co atoms

No lattice.nml needed.

Atom Parameters: Fe-imp_scf.nml
================================

The impurity atom has **modified** parameters vs bulk Fe:

.. code-block:: fortran

   &element
   symbol = 'Fe-imp'              ! Distinct from bulk Fe
   atomic_number = 26
   core = 18
   valence = 8
   /
   
   &par
   lmax = 2
   
   ! Modified parameters due to:
   ! - Co nearest neighbors (instead of Fe+Co mix)
   ! - Local lattice distortion
   ! - Different magnetic environment
   
   center_band(:, 1) = -0.313, 0.338, -0.244
   center_band(:, 2) = -0.274, 0.389, -0.041
   
   width_band(:, 1) = 0.401, 0.261, 0.115
   width_band(:, 2) = 0.400, 0.265, 0.139
   
   ql(1, :, 1) = 0.335, 0.372, 4.180  ! Modified from bulk Fe
   ql(1, :, 2) = 0.360, 0.448, 2.285
   
   ws_r = 2.657329
   /

**Key point:** Even though chemically Fe, the impurity has different tight-binding parameters than bulk Fe due to local environment.

Running the Calculation
=======================

.. code-block:: bash

   cd example/impurity/B2FeCo
   ../../../../build/bin/rslmto.x

Expected Output
===============

**Convergence:** 30-60 iterations.

**Impurity magnetic moment:**

.. code-block:: bash

   grep "magnetic moment" Fe-imp_out.nml
   grep "magnetic moment" Fe_out.nml      # Compare with bulk Fe

Expected: Impurity moment ~2.5-2.8 μ_B (enhanced vs bulk Fe ~2.2 μ_B).

**Local charge redistribution:**

.. code-block:: bash

   grep "charge" Fe-imp_out.nml
   grep "charge" Fe_out.nml

**Expected runtime:** 2-4 minutes.

Physical Interpretation
========================

Fe impurity in FeCo shows:

- **Enhanced magnetic moment** due to Co neighbors (more electrons available for polarization)
- **Local DOS modification** (impurity-induced resonances)
- **Scattering effects** (contributes to residual resistivity)

This is relevant for:

- Understanding disorder in alloys
- Magnetic damping mechanisms
- Residual resistivity calculations

See Also
========

- :doc:`bulk_bccfe` - Bulk magnetic calculation
- :doc:`surface_fcccu001` - Another type of inhomogeneous system
- :doc:`../../keywords/lattice_geometry` - nclu, inclu parameters
- :doc:`../../theory/green_functions` - Embedded cluster method
