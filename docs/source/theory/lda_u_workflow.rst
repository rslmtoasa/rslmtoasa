.. _theory_lda_u_workflow:

=======================================
LDA+U and LDA+U+V Workflow
=======================================

Purpose
=======

This page documents the current LDA+U and LDA+U+V implementation in RS-LMTO-ASA:

1. Reference equations and conventions used in development.
2. Developer-oriented mapping from equations to Fortran routines and variables.
3. Practical ``input.nml`` examples for the supported modes.

Primary authorship credit
=========================

The LDA+U/LDA+U+V development documented here is based on work carried out by
**Viktor Frilen** and **Emil Beiersdorf** (project documentation, June 2024).
Please keep this attribution when reusing or extending this workflow description.

Credits and references
======================

This implementation and documentation build on the internal project document by
Viktor Frilen and Emil Beiersdorf (June 2024), and on the standard literature:

- Liechtenstein, Anisimov, Zaanen (rotationally invariant DFT+U).
- Dudarev et al. (``U_eff = U - J`` form).
- Campo and Cococcioni (extended DFT+U+V form).
- Agapito, Curtarolo, Nardelli (ACBN0 framework).

Reference equations
===================

Rotationally invariant on-site interaction (Liechtenstein form):

.. math::

   E_U = \frac{1}{2}\sum_{m m' m'' m''' \sigma}
   \left[
   \langle m m''|V_{ee}|m' m'''\rangle n^{\sigma}_{m m'} n^{-\sigma}_{m'' m'''}
   +
   \left(
   \langle m m''|V_{ee}|m' m'''\rangle
   -
   \langle m m''|V_{ee}|m''' m'\rangle
   \right)
   n^{\sigma}_{m m'} n^{\sigma}_{m'' m'''}
   \right].

Double counting (FLL):

.. math::

   E_{dc} =
   \frac{U}{2}n(n-1) -
   \frac{J}{2}\sum_{\sigma} n^{\sigma}(n^{\sigma}-1).

Resulting potential matrix (on-site channel):

.. math::

   V^{\sigma}_{m m'} =
   \frac{\partial (E_U - E_{dc})}{\partial n^{\sigma}_{m m'}}.

Simplified Dudarev form:

.. math::

   E_U - E_{dc} = \frac{U_{\mathrm{eff}}}{2}\sum_{\sigma}
   \mathrm{Tr}\left[n^{\sigma}(1-n^{\sigma})\right], \quad
   U_{\mathrm{eff}} = U - J.

Extended +V correction (nearest-neighbor inter-site contribution):

.. math::

   E_{U+V} - E_{dc} =
   \sum_{I,\sigma}\frac{U^I}{2}\mathrm{Tr}[n^{I\sigma}(1-n^{I\sigma})]
   -
   \sum_{I,J,\sigma}^{*}\frac{V^{IJ}}{2}\mathrm{Tr}[n^{IJ\sigma}n^{JI\sigma}].

Implementation status
=====================

- ``liechtenstein`` potential path:
  fully wired and default.
- ``acbn0`` potential path:
  optional path is available, but should currently be treated as an
  implementation-in-progress path for validation studies.
- ``+V``:
  implemented as a nearest-neighbor orbital-channel correction using local
  occupations and a practical diagonal-in-m approximation.

Code mapping (developer view)
=============================

Main entry points
-----------------

- ``source/hamiltonian.f90``: ``calculate_hubbard_u_potential_general`` builds
  the on-site Hubbard potential matrix.
- ``source/hamiltonian.f90``: ``calculate_hubbard_v_potential`` builds the
  inter-site nearest-neighbor +V potential contribution.
- ``source/hamiltonian.f90``: ``build_bulkham`` and ``build_locham`` inject
  ``hubbard_u_pot`` and ``hubbard_v_pot`` into Hamiltonian blocks.
- ``source/bands.f90``: ``calculate_hubbard_u_sc`` performs the self-consistent
  ``U_eff`` update loop for channels selected by mask.
- ``source/math.f90``: ``Coulomb_mat`` and ``a_k`` provide Coulomb matrix
  elements in real spherical-harmonic form.
- ``source/math.f90``: ``tabulated_slater_integrals`` provides tabulated radial
  integral values used by SC-U flow.

Core arrays and variables
-------------------------

- ``ldm(na,l,ispin,m,m')``: local density matrix per atom type and channel.
- ``hubbard_u(l)``, ``hubbard_j(l)``: on-site U/J per ``l`` channel.
- ``hubbard_u_pot(i,j,na)``: assembled on-site correction in basis space.
- ``hubbard_v(i,j,li,lj)`` and ``hubbard_v_pot``: inter-site correction data.
- ``hubbard_u_sc(itype,l)``: mask enabling SC-U update for a channel.
- ``hubbard_u_potential_form``:
  ``'liechtenstein'`` (default) or ``'acbn0'``.

Namelist inputs
===============

The relevant keys are in ``&hamiltonian``.

1. Fixed on-site U/J (recommended baseline)
-------------------------------------------

.. code-block:: fortran

   &hamiltonian
     hubbard_u_general(1,:) = 0.0, 0.0, 2.5, 0.0
     hubbard_j_general(1,:) = 0.0, 0.0, 0.9, 0.0
     hubbard_u_potential_form = 'liechtenstein'
   /

2. Self-consistent U update (SC-U)
----------------------------------

.. code-block:: fortran

   &hamiltonian
     ! Enable SC-U on d channel for atom type 1.
     hubbard_u_sc(1,3) = 1
     ! Do not combine with explicit hubbard_u_general/hubbard_j_general.
     hubbard_u_potential_form = 'liechtenstein'
   /

3. Optional ACBN0-like potential path
-------------------------------------

.. code-block:: fortran

   &hamiltonian
     hubbard_u_sc(1,3) = 1
     hubbard_u_potential_form = 'acbn0'
   /

4. Add nearest-neighbor +V
--------------------------

.. code-block:: fortran

   &hamiltonian
     hubbard_u_general(1,:) = 0.0, 0.0, 2.5, 0.0
     hubbard_j_general(1,:) = 0.0, 0.0, 0.9, 0.0
     ! Example: V between type 1 d and type 1 d (l-index: s=1, p=2, d=3, f=4)
     hubbard_v(1,1,3,3) = 0.6
   /

Notes and constraints
=====================

- Internal energy unit is Ry; user inputs in namelists are interpreted in eV and converted.
- ``hubbard_u_sc`` and explicit ``hubbard_u_general/hubbard_j_general`` are mutually exclusive in the current flow.
- For ``lmax=3`` (spdf), orbital blocks are sized from the active basis dimensions (no hardwired spd-only matrix extents in the current implementation path).
