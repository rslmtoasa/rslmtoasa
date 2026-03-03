.. _examples/bulk_bccfe:

===================================
Bulk bcc Iron (Ferromagnetic)
===================================

**Location:** ``example/bulk/bccFe/``

**System:** Body-centered cubic Fe with collinear ferromagnetism

**Physics:** Demonstrates spin-polarized SCF calculation for a ferromagnetic transition metal.

Overview
========

Iron crystallizes in the BCC structure at room temperature and is ferromagnetic with a magnetic moment of approximately 2.2 Bohr magnetons per atom. This example demonstrates:

- Spin-polarized calculation with ``nsp = 2``
- Magnetic SCF convergence
- BCC lattice generation (no separate lattice.nml needed)
- Careful mixing for magnetic stability

Files
=====

- ``input.nml`` - Main calculation parameters
- ``Fe.nml`` - Spin-dependent tight-binding parameters

Input File: input.nml
=====================

.. code-block:: fortran

   &calculation
   pre_processing = 'bravais'
   verbose = T
   /
   
   &lattice
   rc = 80                       ! Larger cluster for magnetic system
   alat = 2.86120                ! BCC lattice constant (Å)
   crystal_sym = 'bcc'           ! Body-centered cubic (auto-generated)
   wav = 1.40880                 ! Wigner-Seitz radius
   ntype = 1                     ! Single atom type
   ct(1) = 3.0d0                 ! Charge transfer parameter
   r2 = 9.00d0                   ! Screening radius
   /
   
   &atoms
   database = './'
   label(1) = 'Fe'
   /
   
   &self
   nstep = 100                   ! More iterations for magnetic convergence
   /
   
   &energy
   fermi = -0.069282
   energy_min = -1.0
   energy_max = 0.5
   channels_ldos = 2500
   /
   
   &control
   calctype = 'B'
   nsp = 2                       ! Spin-polarized (up and down channels)
   lld = 21
   recur = 'block'
   /
   
   &mix
   beta = 0.01                   ! Small mixing for magnetic stability
   mixtype = 'broyden'
   /

Key Parameters for Magnetic Systems
------------------------------------

**nsp = 2:** Enables spin-polarization with separate up and down electron channels. Essential for ferromagnets.

**beta = 0.01:** Very small mixing parameter. Magnetic systems are sensitive to mixing - too large beta causes oscillations and non-convergence.

**nstep = 100:** Magnetic systems typically need more iterations than non-magnetic ones.

**crystal_sym = 'bcc':** BCC structure is generated automatically - no lattice.nml file needed.

Lattice: Auto-Generated BCC
============================

With ``crystal_sym = 'bcc'``, the code generates:

- Simple cubic cell with 1-atom basis
- Atom at origin: (0,0,0)
- 8 nearest neighbors at distance a*sqrt(3)/2 = 2.48 Å

Atom Parameters: Fe.nml
========================

Fe.nml contains **spin-dependent** tight-binding parameters showing up and down spin channels with different d-electron occupations.

Running the Calculation
=======================

.. code-block:: bash

   cd example/bulk/bccFe
   ../../../build/bin/rslmto.x

Expected Output
===============

**Convergence:** 30-50 iterations (slower than non-magnetic Si).

**Magnetic moment:**

.. code-block:: bash

   grep "magnetic moment" Fe_out.nml

Expected: 2.15-2.25 Bohr magnetons per Fe atom (experimental: 2.22).

**Fermi energy:**

.. code-block:: bash

   grep "fermi" input_out.nml

Expected: fermi = -0.069 ± 0.005 Ry.

**Spin-split DOS:** Plot DOS to see exchange splitting (up and down bands shifted).

Troubleshooting
===============

**Oscillating magnetic moment:**

- Reduce beta to 0.005 or even 0.001
- Increase nstep to 200

**Wrong magnetic moment:**

- Check Fe.nml has correct spin-dependent parameters
- Verify fermi energy is reasonable
- Try starting from different initial fermi value

See Also
========

- :doc:`bulk_si` - Non-magnetic bulk example
- :doc:`bulk_mn3sn` - Non-collinear magnetism
- :doc:`exchange_bccfe` - Calculate exchange interactions for this system
- :doc:`../../keywords/control_parameters` - nsp parameter details
