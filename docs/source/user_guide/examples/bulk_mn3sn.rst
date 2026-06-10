.. _examples/bulk_mn3sn:

==========================================
Bulk Mn₃Sn (Non-Collinear Magnetism)
==========================================

**Location:** ``example/bulk/Mn3Sn/120/``

**System:** Mn₃Sn with hexagonal structure and 120° non-collinear magnetic order

**Physics:** Demonstrates fully relativistic calculation with non-collinear magnetism (``nsp = 4``).

Overview
========

Mn₃Sn forms a kagome lattice with antiferromagnetic 120° spin structure. This is a topological material showing anomalous Hall effect despite near-zero net magnetization. The calculation requires:

- ``nsp = 4`` - Fully relativistic (spinor formalism)
- Complex lattice structure (8 atoms per unit cell)
- Non-collinear magnetic ordering

Files
=====

- ``input.nml`` - Main parameters
- ``lattice.nml`` - Hexagonal structure with 8 atoms
- ``Mn1.nml`` through ``Mn6.nml`` - Six inequivalent Mn sites
- ``Sn1.nml``, ``Sn2.nml`` - Two Sn sites

Input File: input.nml
=====================

.. code-block:: fortran

   &calculation
   pre_processing = 'bravais'
   verbose = T
   /
   
   &lattice
   rc = 60
   alat = 5.677000               ! Hexagonal lattice constant
   crystal_sym = 'file'          ! Complex structure from lattice.nml
   wav = 1.5261136945
   ntype = 8                     ! 6 Mn + 2 Sn
   ct(1:8) = 4.0d0
   r2 = 16.0d0
   /
   
   &atoms
   database = './'
   label(1) = 'Sn1'
   label(2) = 'Sn2'
   label(3) = 'Mn1'              ! Six Mn sites with different
   label(4) = 'Mn2'              ! moment directions
   label(5) = 'Mn3'
   label(6) = 'Mn4'
   label(7) = 'Mn5'
   label(8) = 'Mn6'
   /
   
   &self
   nstep = 100
   /
   
   &energy
   fermi = -0.091448
   energy_min = -1.2
   energy_max = 1.0
   channels_ldos = 2500
   /
   
   &control
   calctype = 'B'
   nsp = 4                       ! ← Fully relativistic (4 spinor components)
   lld = 21
   recur = 'block'
   sym_term = .true.
   /
   
   &mix
   beta = 0.10
   mixtype = 'broyden'
   magbeta(1:8) = 8*0.0d0        ! No separate magnetic mixing
   /

Key Parameter: nsp = 4
----------------------

``nsp = 4`` enables fully relativistic treatment with 4-component spinors. This is **required** for non-collinear magnetism where magnetic moments point in different directions (not just ↑ or ↓).

Lattice File: lattice.nml
==========================

The hexagonal structure with kagome planes:

.. code-block:: fortran

   &lattice
   nbulk_bulk = 8
   ntot = 8
   nbas = 8
   nrec = 8
   
   ! Hexagonal lattice vectors
   a(:, 1) = 0.866025403784, -0.500000000000,  0.000000000000
   a(:, 2) = 0.000000000000,  1.000000000000,  0.000000000000
   a(:, 3) = 0.000000000000,  0.000000000000,  0.798661264753
   
   ! Atomic positions (8 atoms: 6 Mn + 2 Sn)
   crd(:, 1) = 0.577350269190,  0.000000000000,  0.199665316188  ! Sn1
   crd(:, 2) = 0.288675134595,  0.500000000000,  0.598995948564  ! Sn2
   crd(:, 3) = 0.144337567297,  0.750000000000,  0.199665316188  ! Mn1
   crd(:, 4) = 0.144337567297,  0.250000000000,  0.199665316188  ! Mn2
   crd(:, 5) = 0.577350269190,  0.500000000000,  0.199665316188  ! Mn3
   crd(:, 6) = 0.721687836487, -0.250000000000,  0.598995948564  ! Mn4
   crd(:, 7) = 0.721687836487,  0.250000000000,  0.598995948564  ! Mn5
   crd(:, 8) = 0.288675134595,  0.000000000000,  0.598995948564  ! Mn6
   
   izp(1:8) = 1, 2, 3, 4, 5, 6, 7, 8
   ! ... (additional indexing arrays)
   /

**Structure:** 

- Two kagome planes of Mn at different z-coordinates
- Sn atoms between planes
- Each Mn forms triangular network in xy-plane

Atom Parameters: Example Mn1.nml
=================================

.. code-block:: fortran

   &element
   symbol = 'Mn1'
   atomic_number = 25
   core = 18
   valence = 7                   ! 3d⁵4s²
   /
   
   &par
   lmax = 2
   
   ! Spin-up and spin-down parameters
   pl(:, 1) = 4.5899147462, 4.3260351472, 3.6404205929
   pl(:, 2) = 4.6974968083, 4.3466020354, 3.8509487007
   
   center_band(:, 1) = -0.164125105, 0.360414798, 0.144109265
   center_band(:, 2) = -0.191082914, 0.350365981, -0.013903302
   
   width_band(:, 1) = 0.392011514, 0.237280570, 0.153006356
   width_band(:, 2) = 0.424721639, 0.244849087, 0.136863500
   
   ql(1, :, 1) = 0.2874611428, 0.2851651002, 4.2983834574
   ql(1, :, 2) = 0.2999985272, 0.2737097583, 1.0855169725
   
   ws_r = 2.656529
   /

**Magnetic moment:** Each Mn site has local moment ~3 μ_B, but directions vary following 120° structure.

Running the Calculation
=======================

.. code-block:: bash

   cd example/bulk/Mn3Sn/120
   ../../../../build/bin/rslmto.x

Expected Output
===============

**Convergence:** 40-80 iterations (non-collinear is computationally demanding).

**Magnetic structure:**

.. code-block:: bash

   grep "moment" Mn*_out.nml

Each Mn should show moment components in x, y, z (not just along one axis).

**Expected runtime:** 5-10 minutes.

Physical Interpretation
========================

Mn₃Sn exhibits:

- **Large anomalous Hall effect** despite small net magnetization
- **Weyl points** in electronic structure
- **Berry curvature** from non-collinear spins and spin-orbit coupling

This is a prototypical topological antiferromagnet.

See Also
========

- :doc:`bulk_bccfe` - Collinear magnetism (nsp=2)
- :doc:`conductivity_bccfe` - Anomalous Hall effect
- :doc:`../../keywords/control_parameters` - nsp parameter
- :doc:`../../theory/spin_dynamics` - Non-collinear magnetism theory
