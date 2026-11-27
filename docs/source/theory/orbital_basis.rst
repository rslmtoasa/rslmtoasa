Orbital basis conventions
=========================

This page lists the orbital orderings used in the code. Two bases occur: the
real cubic order of the two-centre construction and the spherical basis of
everything else.

Spherical basis (Hamiltonian, Green function, density matrix)
-------------------------------------------------------------

The matrices ``ee``, ``hall``, ``lsham``, ``obarm``, ``enim``, ``g0``, ``ldm``,
the k-space Hamiltonian, the canonical export (``rs_tb.tb``, ``rs_tb.kspace``)
and ``nonlocal_dm*.out`` use complex harmonics.

* Basis functions: :math:`(-1)^m Y_{lm}` with Condon-Shortley :math:`Y_{lm}`,
  equivalently :math:`\overline{Y_{lm}}` with :math:`m` ascending.
* Order within a shell: :math:`m = +l, \ldots, -l` (descending).
* Spin layout: per site the orbitals ``1..(lmax+1)**2`` are spin up, followed by
  the same list for spin down (``spin_off = (lmax+1)**2``).
* Shells follow each other in order of :math:`l`; shell :math:`l` occupies the
  indices ``l*l+1 .. (l+1)**2``.

.. list-table::
   :header-rows: 1

   * - index
     - l
     - m
     - cubic orbitals it mixes
   * - 1
     - 0
     - 0
     - s
   * - 2
     - 1
     - +1
     - x, y
   * - 3
     - 1
     - 0
     - z
   * - 4
     - 1
     - -1
     - x, y
   * - 5
     - 2
     - +2
     - xy, x2-y2
   * - 6
     - 2
     - +1
     - yz, zx
   * - 7
     - 2
     - 0
     - 3z2-r2
   * - 8
     - 2
     - -1
     - yz, zx
   * - 9
     - 2
     - -2
     - xy, x2-y2

For f the indices 10..16 run over :math:`m = +3, \ldots, -3`; the m = 0 member
(index 13) is ``fz3``.

Cubic order (two-centre blocks, legacy export)
----------------------------------------------

Real orbitals: ``s | x y z | xy yz zx x2-y2 3z2-r2 | fz3 fxz2 fyz2 fz(x2-y2) fxyz fx3 fy3``
(``canso`` in ``lattice_strux.f90``). It is used by the structure constants,
the hopping blocks ``hhh``/``hhmag``, the pair blocks of the GBT builders and
the real matrices ``L_x, L_y, L_z`` before conversion.

Conversion points:

* ``hcpx(..., 'cart2sph')`` (``math.f90``) is the single conversion; it acts on
  one spin block of ``norb x norb`` at a time. Hopping blocks are converted in
  ``chbar_nc`` (``hamiltonian_build.f90``); on-site operators (``obarm``,
  ``enim``, L matrices, ``lsham``) and the k-space GBT blocks at their build.
* The legacy seven-column export (``rs_tb_paoham.dat``, ``rs2paoham.dat``) is
  written in cubic order: ``hamiltonian_export.f90`` applies
  ``hcpx(..., 'sph2cart')`` before writing.
* The point-group symmetrization of the site density matrix
  (``reciprocal_symmetry.f90``) goes to cubic order, averages with real
  rotation matrices, and converts back.

Real m order
------------

Real harmonics in the order :math:`m = -l, \ldots, l` (l = 1: y z x) are used by
``strux_lib`` (mapped to cubic order by ``build_orbital_map``), by ``realgaunt``
and ``Coulomb_mat`` (``math.f90``), and by ``disp_matrix`` for :math:`l \ge 3`
(for :math:`l \le 2` it returns cubic order).
``Coulomb_mat_spherical`` is the complex-harmonic version.

Importers
---------

.. list-table::
   :header-rows: 1

   * - importer
     - file basis expected
     - transform on read
   * - ``paoflow2rs`` (``build_from_paoflow_opt``)
     - spherical
     - none
   * - ``effective_import = 'legacy'``
     - cubic
     - ``hcpx('cart2sph')``
   * - ``effective_import = 'canonical'``
     - spherical
     - none

The native legacy export is in cubic order, so it is read back correctly by the
``legacy`` importer and not by ``paoflow2rs``. The basis of external PAOFLOW
files is not documented here.

Known mismatches
----------------

Details and measurements are in ``tests/KNOWN_ISSUES.md``, entry "Orbital basis
audit".

* ``ROTMAT`` (``local_axis``): the orbital block is the conjugate of the
  rotation in the spherical basis for azimuth different from zero.
* ACBN0 and SC-U contract the real-harmonic ``Coulomb_mat`` with the
  spherical-basis density matrix.
* ``paoflow2rs`` takes its file in the spherical basis without a transform.
* ``disp_matrix`` (``symbolic_atom.f90``) is in cubic or real m order and is
  multiplied with spherical-basis Green functions in ``calculate_jijk``; its
  header comment lists the cubic order only.
* The metadata strings ``..._l_then_m_increasing`` in the export files name the
  opposite in-shell order; the files are in the order given above.
