.. _keywords/hamiltonian_parameters:

===================================================
Hamiltonian Parameters (&hamiltonian Namelist)
===================================================

Overview
========

The ``&hamiltonian`` namelist controls how the TB-LMTO Hamiltonian is
assembled: the operator form (two-centre vs. H-O-H), the combined-correction
(CCOR) terms, Hubbard :math:`+U`/:math:`+V` corrections, spin-spiral / local-axis
options, and Hamiltonian export. It is read by
``hamiltonian%build_from_file`` (``source/hamiltonian_build.f90``).

.. note::

   The Hubbard :math:`+U`/:math:`+V` keywords (``hubbard_u``, ``hubbard_j``,
   ``hubbard_v``, ``hubbard_u_potential_form``, and the impurity/general
   variants) are documented alongside the workflow that uses them — see
   :ref:`theory_lda_u_workflow`. This page covers the Hamiltonian *form* and
   correction switches. Hubbard inputs in this namelist are given in **eV**;
   internal storage is in Ry.

Hamiltonian form
================

hoh
---

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Use the "H-O-H" (overlap-corrected, three-centre) form of the
Hamiltonian operator instead of the plain two-centre form.

**Notes:**

- Changes the operator applied in the recursion — the fast Chebyshev and
  block-Lanczos kernels implement a dedicated two-sweep (bond + onsite + hop)
  ``hoh`` path.
- The k-space Fourier assembly also branches on ``hoh``.

**Related code:** ``source/hamiltonian_build.f90``, ``source/chebyshev_fast.f90``

**See also:** :ref:`theory/recursion_method`

local_axis
----------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Build the Hamiltonian in a per-site local spin-quantization axis
(for non-collinear magnetism).

orb_pol
-------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Add orbital-polarization corrections to the Hamiltonian to enhance
orbital magnetism.

Combined correction (CCOR)
==========================

The combined correction restores accuracy lost by the ASA at higher energies
by adding second-order correction blocks to the Hamiltonian
(``source/hamiltonian_ccor.f90``). On the k-space route the Bloch-summed
``Hcc_2c(k)`` term is added when ``ccor_2c = .true.`` (see
:ref:`theory/kspace_modes`).

ccor_2c
-------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Enable the second-order combined-correction term.

**Example:**

.. code-block:: fortran

   &hamiltonian
      ccor_2c = .true.
   /

ccor_elin
---------

**Type:** Real (Ry)

**Default:** ``0.0``

**Purpose:** Linearization energy about which the combined correction is
expanded.

ccor_vmt_mode
-------------

**Type:** Character string

**Default:** ``'surface_scalar'``

**Allowed values:** ``'surface_scalar'``, ``'vmad_scalar'``, ``'pair_surface'``

**Purpose:** Selects how the muffin-tin zero / potential shift entering the
combined correction is evaluated. An invalid value is fatal at namelist read.

ccor_debug / ccor_strict
------------------------

**Type:** Logical

**Default:** ``.false.`` / ``.false.``

**Purpose:** ``ccor_debug`` emits additional CCOR diagnostics; ``ccor_strict``
enforces stricter consistency checks on the CCOR inputs.

Spin-spiral / non-collinear
===========================

- ``q_ss`` (real(3), default ``[0,0,0]``) — spin-spiral wavevector, in
  Cartesian units of :math:`2\pi/\text{alat}` (``q_ss=0.5`` along a cubic
  axis is the zone boundary at :math:`\pi/\text{alat}`).
- ``theta_ss`` (real, default ``0``) — spin-spiral cone angle, degrees in
  the namelist (converted to radians internally).
- ``magnetic_representation`` (character, default ``'periodic_nc'``) — selects
  ``'periodic_nc'``, ``'gbt_single_q'``, or ``'explicit_texture'`` independently
  of whether recursion or reciprocal space solves the resulting operator.

  .. warning::

     ``q_ss`` is only read by ``'gbt_single_q'`` (as the bond gauge phase) and
     ``'explicit_texture'`` (as the site-dependent moment direction). Under the
     **default** ``'periodic_nc'`` it is stored and never used, so a spiral
     sweep would return the same collinear answer at every ``q``. Setting a
     nonzero ``q_ss`` without selecting one of the two spiral representations
     is therefore **fatal** at ``build_bulkham``; the two live interpretations
     are different physics (the gauge trick versus the real-space texture it is
     validated against), so the code names all three exits rather than guessing.
     ``theta_ss`` alone does not trigger it: a cone angle at ``q = 0`` is a
     global spin rotation and is energy-invariant without SOC.
     ``post_processing='frozen_magnon'`` selects ``'gbt_single_q'`` itself and
     is unaffected.
- ``v_alpha`` / ``v_beta`` (real(3)) — quantization-axis vectors.
- ``js_alpha`` / ``jl_alpha`` (character, default ``'z'``) — spin / orbital
  projection axes.
- ``gbt_kspace`` (logical, default ``.false.``) — deprecated compatibility
  input. ``.true.`` emits a warning and is ignored; select
  ``magnetic_representation='gbt_single_q'`` explicitly. Both solvers consume the same
  linked real-space operator; reciprocal space always uses the ordinary
  Fourier transform.

Hamiltonian export
==================

export
------

**Type:** Character string

**Default:** ``'none'``

**Allowed values:** ``'none'``, ``'python'``, ``'rs2pao'``

**Purpose:** After building the real-space Hamiltonian, optionally write it to
disk in a PAOFLOW-compatible format (the same "legacy7" record layout that
``post_processing = 'paoflow2rs'`` reads back).

**Notes:**

- ``'python'`` writes ``<base>_paoham.dat`` via ``export_rs_tb_all``.
- ``'rs2pao'`` writes via the ``rs2pao`` path.
- ``'none'`` disables export.

**Related code:** ``source/hamiltonian_paoflow_io.f90``,
``source/calculation.f90::pre_processing_bravais()``

**See also:** :ref:`code_structure`

Provenance
==========

- **Definition:** ``source/hamiltonian.f90::type hamiltonian``
- **Namelist:** ``source/include_codes/namelists/hamiltonian.f90``
- **Reading:** ``source/hamiltonian_build.f90::build_from_file()``

See Also
========

- :ref:`theory/kspace_modes` - CCOR on the k-space route
- :ref:`theory_lda_u_workflow` - Hubbard :math:`+U`/:math:`+V` inputs
- :ref:`keywords/control_parameters` - recursion / GPU switches
