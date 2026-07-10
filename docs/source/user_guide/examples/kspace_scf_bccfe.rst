.. _examples/kspace_scf_bccfe:

==========================================
K-Space-Enabled SCF for bcc Iron
==========================================

**Based on:** ``example/bulk/bccFe/`` (see :doc:`bulk_bccfe`)

**System:** Body-centered cubic Fe, collinear ferromagnet

**Physics:** Demonstrates the optional reciprocal-space SCF branch, in which
the DOS, band moments, Fermi level, and any LDA+U density-matrix updates are
taken from a k-space diagonalization workflow instead of the real-space
recursion. Useful as a bulk validation cross-check against the recursion path.

Overview
========

By default RS-LMTO computes the SCF quantities with the real-space recursion
method. For a periodic bulk system you can instead route the per-iteration DOS
and moments through the reciprocal-space workflow by setting
``use_kspace = .true.`` in ``&self``. The recursion path remains the production
default; this branch is intended for bulk validation.

See :ref:`theory/kspace_modes` for the full k-space option inventory.

Input File: input.nml
=====================

This is the :doc:`bulk_bccfe` input with an added ``use_kspace`` switch and a
``&reciprocal`` namelist controlling the k-mesh:

.. code-block:: fortran

   &calculation
   pre_processing = 'bravais'
   verbose = T
   /

   &lattice
   rc = 80
   alat = 2.86120
   crystal_sym = 'bcc'
   wav = 1.40880
   ntype = 1
   ct(1) = 3.0d0
   r2 = 9.00d0
   /

   &atoms
   database = './'
   label(1) = 'Fe'
   /

   &self
   nstep = 100
   use_kspace = .true.           ! Take DOS/moments/Fermi from the k-space route
   /

   &reciprocal
   kspace_ham_order = 'second'   ! Second-order ASA H(k)
   nk1 = 16, nk2 = 16, nk3 = 16  ! Monkhorst-Pack mesh
   use_symmetry_reduction = .true.
   use_time_reversal = .true.
   dos_method = 'tetrahedron'    ! or 'blochl' / 'gaussian'
   /

   &energy
   fermi = -0.069282
   energy_min = -1.0
   energy_max = 0.5
   channels_ldos = 2500
   /

   &control
   calctype = 'B'
   nsp = 2
   lld = 21
   recur = 'block'
   /

   &mix
   beta = 0.01
   mixtype = 'broyden'
   /

Key Parameters
==============

**use_kspace = .true.** (``&self``)
   Enables the optional k-space DOS/moment branch. With it set, the per-iteration
   DOS, band moments, Fermi level and LDA+U LDM updates come from the
   reciprocal-space workflow rather than the recursion.

**kspace_ham_order = 'second'** (``&reciprocal``)
   Builds the second-order ASA Hamiltonian
   :math:`H(k) = h(k) - [hoh](k) + E_\nu + H_{LS}`, adding the Bloch-summed
   ``Hcc_2c(k)`` combined correction when ``ccor_2c = .true.``. Use ``'first'``
   for the lower-order ``h(k) + H_LS`` form.

**nk1/nk2/nk3**
   Monkhorst-Pack mesh density. Denser meshes improve DOS/Fermi-level accuracy
   at the cost of more diagonalizations.

**dos_method**
   ``'tetrahedron'``, ``'blochl'`` (Blöchl-improved tetrahedron), or
   ``'gaussian'`` broadening.

Running the Calculation
=======================

.. code-block:: bash

   cd example/bulk/bccFe
   # edit input.nml to add use_kspace / &reciprocal as above
   ../../../build/bin/rslmto.x

Expected Output
===============

The converged magnetic moment and Fermi energy should agree with the recursion
run in :doc:`bulk_bccfe` (about 2.2 μ_B/atom, ``fermi ≈ -0.069 Ry``) to within
the k-mesh discretization error — that agreement is the point of the
cross-check. The k-space DOS is written to ``dos_kspace.dat``.

Notes
=====

- This branch is best suited to periodic bulk systems; surface/impurity
  geometries stay on the recursion route.
- Requires ``ENABLE_SPGLIB=ON`` for symmetry-reduced meshes and k-path
  generation (``use_symmetry_reduction = .true.``).

See Also
========

- :ref:`theory/kspace_modes` - K-space Hamiltonian options and inventory
- :doc:`bulk_bccfe` - The real-space recursion version of this system
- :doc:`../../keywords/hamiltonian_parameters` - ``ccor_2c`` combined correction
