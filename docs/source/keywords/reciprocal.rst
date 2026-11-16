.. _keywords/reciprocal:

=====================================================
K-Space Parameters (&reciprocal, &kpath)
=====================================================

Overview
========

``&reciprocal`` configures every k-space route: the k-space SCF
(``use_kspace = .true.`` in ``&self``), the post-processing modes
``'band_structure'``, ``'density_of_states'``, ``'fermi_surface'`` and
``'bsf'``, and the k-space Green functions used by ``gf_route = 'lehmann'`` or
``'dyson'`` (see :ref:`keywords/output_options`). ``&kpath`` sets the band-structure
path. Both namelists are optional; a read error in either stops the run.

The Hamiltonian forms selected by ``kspace_ham_order`` and ``reciprocal_mode``
are described in :ref:`theory/kspace_modes`. All energies are in Ry.

.. code-block:: fortran

   &reciprocal
      nk1 = 16, nk2 = 16, nk3 = 16
      use_symmetry_reduction = .true.
      dos_method = 'blochl'
      dos_energy_min = -0.8
      dos_energy_max = 0.9
      n_energy_points = 2000
   /

   &kpath
      auto_kpath = .true.
      nk_per_segment = 40
   /

A finite spin-spiral vector (``q_ss`` in ``&hamiltonian``) switches
``use_symmetry_reduction`` and ``use_time_reversal`` off for the run.

K-point mesh and symmetry
=========================

.. list-table::
   :header-rows: 1
   :widths: 24 12 18 46

   * - Key
     - Type
     - Default
     - Meaning
   * - ``nk1``, ``nk2``, ``nk3``
     - integer
     - 8, 8, 8
     - Monkhorst-Pack mesh along the three reciprocal vectors.
   * - ``k_offset_x``, ``k_offset_y``, ``k_offset_z``
     - real
     - 0.0
     - All zero gives the Gamma-centred mesh; any nonzero value selects the shifted
       Monkhorst-Pack mesh. Only whether an offset is nonzero is used; the values are
       written to the Fermi-surface metadata.
   * - ``use_symmetry_reduction``
     - logical
     - ``.true.``
     - Reduce the mesh to the irreducible wedge with the crystal symmetry (spglib).
   * - ``use_time_reversal``
     - logical
     - ``.true.``
     - Also use time reversal in the reduction.
   * - ``strict_symmetry_checks``
     - logical
     - ``.false.``
     - Stop, instead of warning, when the full-to-irreducible k map is inconsistent.
   * - ``dump_symmetry_kmap``
     - logical
     - ``.false.``
     - Write the full-to-irreducible k map for inspection.
   * - ``tetra_symmetry_mode``
     - character
     - ``'irreducible_native'``
     - Tetrahedron integration on the irreducible mesh (``'irreducible_native'``) or on the
       expanded full mesh (``'full_expand_ref'``). Other values stop the run.

Density of states
=================

.. list-table::
   :header-rows: 1
   :widths: 24 12 18 46

   * - Key
     - Type
     - Default
     - Meaning
   * - ``dos_method``
     - character
     - ``'tetrahedron'``
     - ``'tetrahedron'``, ``'blochl'`` (corrected tetrahedron) or ``'gaussian'``.
   * - ``n_energy_points``
     - integer
     - 1000
     - Points of the DOS energy grid.
   * - ``dos_energy_min``, ``dos_energy_max``
     - real (Ry)
     - -1.0, 1.0
     - Energy window of the DOS.
   * - ``gaussian_sigma``
     - real (Ry)
     - 0.01
     - Gaussian width, used only with ``dos_method = 'gaussian'``.
   * - ``temperature``
     - real (K)
     - 300.0
     - Electronic temperature of the k-space occupations (Fermi function).
   * - ``auto_find_fermi``
     - logical
     - ``.true.``
     - Find the Fermi level from the electron count.
   * - ``total_electrons``
     - real
     - 0.0
     - Read but not used: the electron count comes from the valence of the atoms
       (see ``tests/KNOWN_ISSUES.md``).

Fermi surface export
====================

.. list-table::
   :header-rows: 1
   :widths: 24 12 18 46

   * - Key
     - Type
     - Default
     - Meaning
   * - ``fs_nk1``, ``fs_nk2``, ``fs_nk3``
     - integer
     - 0
     - Dense mesh of ``post_processing = 'fermi_surface'``; values ``<= 0`` use ``nk1..nk3``.
   * - ``fs_use_symmetry_reduction``
     - logical
     - ``.false.``
     - Read but ignored: the Fermi-surface export always covers the full Brillouin zone.
   * - ``eigenpair_output_file``
     - character
     - ``'fermi_surface'``
     - Base name of the export; the writer appends ``.bin`` and ``.meta``.
   * - ``write_eigenpair_projections``
     - logical
     - ``.false.``
     - Also write site, orbital and local-spin weights of every eigenpair.
   * - ``write_spin_resolved_eigenpairs``
     - logical
     - ``.false.``
     - Also write the up and down sectors, each solved separately.

Hamiltonian and solver
======================

.. list-table::
   :header-rows: 1
   :widths: 24 12 18 46

   * - Key
     - Type
     - Default
     - Meaning
   * - ``kspace_ham_order``
     - character
     - ``'auto'``
     - ``'auto'``, ``'first'``, ``'second'`` or ``'orthogonal'``; ``'proper'`` is a deprecated
       alias of ``'second'``. See :ref:`theory/kspace_modes`.
   * - ``reciprocal_mode``
     - character
     - ``'ham_only'``
     - ``'ham_only'`` or ``'generalized_overlap_proxy'``; ``'generalized_overlap'`` is a
       deprecated alias of the latter, and ``'generalized_overlap_kanpur'`` stops the run.
   * - ``reciprocal_backend``
     - character
     - ``'lapack'``
     - Eigensolver: ``'lapack'`` (``'cpu'`` is accepted as an alias) or ``'cuda'``, which needs
       a build with ``ENABLE_CUDA_RECIPROCAL``.
   * - ``green_eta``
     - real (Ry)
     - 0.02
     - Broadening of the k-space Green function, :math:`z = E + i\eta`.
   * - ``green_backend``
     - character
     - ``'lehmann'``
     - ``'lehmann'`` or ``'dyson'``; selects the k-space Green-function construction.

Diagnostics
===========

.. list-table::
   :header-rows: 1
   :widths: 24 12 18 46

   * - Key
     - Type
     - Default
     - Meaning
   * - ``suppress_internal_logs``
     - logical
     - ``.true.``
     - Suppress the detailed symmetry-analysis log.
   * - ``kanpur_diagnostics``
     - logical
     - ``.true.``
     - Print the mapping of the Hamiltonian blocks to the Kanpur notation.
   * - ``gamma_bounds_diagnostics``
     - logical
     - ``.false.``
     - Print spectral bounds of :math:`H(\Gamma)`.
   * - ``hall_diag_experimental``
     - logical
     - ``.false.``
     - Diagonalise the finite real-space Hamiltonian as a check (experimental).

K-path (&kpath)
===============

.. list-table::
   :header-rows: 1
   :widths: 24 12 18 46

   * - Key
     - Type
     - Default
     - Meaning
   * - ``auto_kpath``
     - logical
     - ``.true.``
     - Build the band-structure path from the space group (spglib).
   * - ``nk_per_segment``
     - integer
     - 40
     - k points per path segment.

**Related code:** ``source/reciprocal_lifecycle.f90`` (``build_from_file``,
``restore_to_default``), ``source/include_codes/namelists/reciprocal.f90``.
