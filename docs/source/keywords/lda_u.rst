.. _keywords/lda_u:

=====================================
LDA+U Parameters (&hamiltonian, &par)
=====================================

Overview
========

The Hubbard keys live in ``&hamiltonian``; the mixing of the occupation matrix
is ``ldm_beta`` in ``&mix``. A converged U run writes its occupation matrices
and the U/J it used into the potential files (``&par``) and, optionally, into
a separate density-matrix file (``&density_matrix``). The workflow, the
energy functional and the radial-moment routes are described in
:ref:`theory_lda_u_workflow`.

All U, J and V values in ``input.nml`` are in **eV**; the code stores them in Ry.
The orbital index ``l`` of every array runs 1..4 for s, p, d, f.

.. code-block:: fortran

   &hamiltonian
      hubbard_u_general(1,:) = 0.0, 0.0, 2.0, 0.0
      hubbard_j_general(1,:) = 0.0, 0.0, 0.9, 0.0
      hubbard_u_potential_form = 'liechtenstein'
   /

   &mix
      ldm_beta = 0.3
   /

&hamiltonian keys
=================

.. list-table::
   :header-rows: 1
   :widths: 34 10 18 38

   * - Key
     - Shape
     - Default
     - Meaning
   * - ``hubbard_u_general``, ``hubbard_j_general``
     - (ntype, 4), eV
     - unset
     - U and J per atom type and l. This is the main fixed-U input.
   * - ``hubbard_u_impurity``, ``hubbard_j_impurity``
     - (nrec, 4), eV
     - unset
     - U and J per recursion site of an impurity or surface cluster (``calctype`` I or S).
   * - ``hubbard_u``, ``hubbard_j``
     - (nrec, 4), eV
     - unset
     - Older input form, used together with ``uj_orb``.
   * - ``uj_orb``
     - (nrec), character
     - ``''``
     - Orbital letters (``s``, ``p``, ``d``, ``f``) that ``hubbard_u``/``hubbard_j`` apply to,
       per recursion site.
   * - ``hubbard_u_sc``
     - (ntype, 4), integer
     - 0
     - 1 switches on the self-consistent U update for that type and l. Cannot be combined
       with ``hubbard_u_general``/``hubbard_j_general``.
   * - ``hubbard_v``
     - (ntype, ntype, 4, 4), eV
     - unset
     - Inter-site V between type pairs and l channels. Experimental: needs
       ``hubbard_experimental``.
   * - ``hubbard_u_potential_form``
     - character
     - ``'liechtenstein'``
     - ``'liechtenstein'`` or ``'acbn0'``; ``'acbn0'`` needs ``hubbard_experimental``.
   * - ``hubbard_radial_moments``
     - character
     - ``'auto'``
     - Route of the d-shell radial moments: ``'local'``, ``'spectral'`` or ``'auto'`` (``local``
       where implemented for the input, otherwise ``spectral``). See
       :ref:`theory_lda_u_workflow`.
   * - ``hubbard_experimental``
     - logical
     - ``.false.``
     - Allow the experimental ACBN0 and +V inputs.
   * - ``hubbard_soc_experimental``
     - logical
     - ``.false.``
     - Accepted and ignored: spinor U at ``nsp >= 2`` is on whenever fixed U is given.
   * - ``hubbard_spinor_dm_diagnostic``
     - logical
     - ``.false.``
     - Extract the full spinor d-shell occupation matrix at ``nsp >= 2`` for the log; not used
       by the Hamiltonian.
   * - ``hubbard_log_level``
     - character
     - ``'summary'``
     - ``'summary'``, ``'detailed'`` or ``'debug'``.

Occupation-matrix mixing (&mix)
===============================

.. list-table::
   :header-rows: 1
   :widths: 34 10 18 38

   * - Key
     - Type
     - Default
     - Meaning
   * - ``ldm_beta``
     - real
     - 1.0
     - Linear mixing of the occupation matrix that builds the U potential; clipped to
       [0, 1]. The same value mixes the spinor occupation at ``nsp >= 2``.

Hubbard data in potential files (&par)
======================================

Written by the code into ``<label>_out.nml`` and read back on restart; not normally
edited by hand.

.. list-table::
   :header-rows: 1
   :widths: 32 14 54

   * - Key
     - Type
     - Meaning
   * - ``hubbard_u``, ``hubbard_j``
     - real (lmax+1)
     - U and J of this type per l (Ry).
   * - ``hubbard_radial_moments``, ``hubbard_backend``
     - character
     - Radial-moment route and backend of the run that wrote the file.
   * - ``hubbard_nsp``
     - integer
     - ``nsp`` of the run that wrote the file; 0 = not recorded.
   * - ``hubbard_axis``
     - real (3)
     - Spin axis of the radial channels of that run; zero = not recorded.
   * - ``ldm_flatten``, ``ldm_flatten_imag``
     - real
     - Real and imaginary parts of the occupation matrix in flattened form; the imaginary
       part is absent in older files.
   * - ``dm_file``
     - character
     - Name of a density-matrix file (``&density_matrix``, below) next to the potential file.
       A missing file gives a warning and the run continues from ``ldm_flatten``.
   * - ``dm_sha256``
     - character
     - Written by older versions; read and ignored.

Density-matrix file (&density_matrix)
=====================================

The file named by ``dm_file`` holds one ``&density_matrix`` group. It is written by
the code. A missing file named in ``dm_file`` gives a warning and ``dm_file`` is
cleared; a file that exists but cannot be opened or read as a ``&density_matrix``
group stops the run.

.. list-table::
   :header-rows: 1
   :widths: 32 14 54

   * - Key
     - Type
     - Meaning
   * - ``label``
     - character
     - Atom label the matrix belongs to.
   * - ``radial_mode``, ``backend``
     - character
     - Radial-moment route and backend of the writing run.
   * - ``hubbard_u_ev``, ``hubbard_j_ev``
     - real (4)
     - U and J per l of the writing run (eV).
   * - ``run_id``, ``iteration``
     - character, integer
     - Identifier and SCF iteration of the writing run.
   * - ``dm_re``, ``dm_im``
     - real (l, spin, spin, m, m)
     - Real and imaginary parts of the occupation matrix.
   * - ``spinor_d``
     - integer
     - 1 when the d block is stored as a full spinor matrix, 0 otherwise.

**Related code:** ``source/hamiltonian_build.f90`` (``build_from_file``),
``source/potential.f90`` (``read_density_matrix_file``), ``source/mix.f90``.
