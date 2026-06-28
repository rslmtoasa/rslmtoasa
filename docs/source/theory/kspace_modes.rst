.. _theory/kspace_modes:

================================
K-Space Hamiltonian Options
================================

Overview
========

The reciprocal-space workflow is controlled from ``&reciprocal`` with
``kspace_ham_order``. All active paths build a standard Hermitian
Hamiltonian and solve with identity overlap.

Available options:

- ``kspace_ham_order='second'``: second-order ASA Hamiltonian.
- ``kspace_ham_order='first'``: first-order ``h(k) + H_LS`` Hamiltonian.
- ``kspace_ham_order='auto'``: use the default production path.
- ``kspace_ham_order='proper'``: deprecated alias for ``second``.

Current Functionality Inventory
===============================

Implemented now
---------------

- Reciprocal lattice generation and Monkhorst-Pack meshes.
- Optional symmetry-reduced k-point mesh generation when symmetry support is available.
- High-symmetry k-path generation and band-structure output.
- Real-space to reciprocal-space Fourier transform for multi-site ``H(k)`` assembly.
- Dense k-point diagonalization with LAPACK ``ZHEEV``.
- DOS workflows:
  
  - tetrahedron
  - Blochl-improved tetrahedron
  - Gaussian broadening

- Orbital/site projected DOS from eigenvectors.
- Fermi-level search from DOS and band-moment integration.
- Optional diagnostics:
  
  - k-space mode summary
  - ``H(Gamma)`` bounds diagnostics
  - experimental finite ``HALL`` diagonalization check

Wired entry points
------------------

- ``post_processing='band_structure'`` routes through reciprocal-space workflow.
- ``post_processing='density_of_states'`` routes through reciprocal-space workflow.
- ``&self use_kspace=.true.`` enables an optional k-space SCF DOS/moment branch.

Hamiltonian Semantics
=====================

``second``
----------

- Builds the second-order base Hamiltonian:

  .. math::

     H(k) = h(k) - [hoh](k) + E_\nu + H_{LS}

- Adds the Bloch-summed ``Hcc_2c(k)`` when ``ccor_2c=.true.``.
- Solves a standard Hermitian eigenproblem with identity overlap.

``first``
---------

- Builds the first-order Hamiltonian ``h(k) + H_LS``.
- Adds the Bloch-summed ``Hcc_2c(k)`` when ``ccor_2c=.true.``.
- Useful for comparisons against the lower-order reciprocal-space path.

``proper``
----------

- Deprecated compatibility spelling for ``second``.

Minimal Example
===============

.. code-block:: fortran

   &calculation
      post_processing = 'band_structure'
   /

   &reciprocal
      kspace_ham_order = 'second'
      nk1 = 16, nk2 = 16, nk3 = 16
      use_symmetry_reduction = .true.
      use_time_reversal = .true.
      kspace_diagnostics = .true.
      suppress_internal_logs = .true.
   /

Optional SCF K-Space Branch
===========================

The SCF module has an optional k-space DOS/moment branch through ``&self``:

.. code-block:: fortran

   &self
      use_kspace = .true.
   /

Notes:

- Recursion remains the default production SCF path.
- With ``use_kspace=.true.``, DOS, band moments, Fermi level, and LDA+U LDM
  updates are taken from the reciprocal-space workflow.
- This branch is best suited for bulk validation runs.

For Chebyshev bounds used by recursion/KPM, see :ref:`theory/recursion_method`.
