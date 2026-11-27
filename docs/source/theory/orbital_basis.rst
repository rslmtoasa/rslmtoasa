.. _theory/orbital_basis:

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
(``math.f90``), and by ``disp_matrix`` for :math:`l \ge 3`
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

On-site multipoles (``multipole_route = 'general'``)
-----------------------------------------------------

The on-site density matrix :math:`n_{ab}` is taken in the spherical basis above, with :math:`b_a = (-1)^m Y_{lm}`,
and the sphere density is :math:`\rho(\mathbf r) = \sum_{ab} \psi_a\, n_{ab}\, \overline{\psi_b}`,
:math:`\psi_a = (u_{l_a}(r)/r)\, b_a(\hat r)`. The multipoles are
:math:`Q_{LM} = \int \rho\, r^L R_{LM}\, d^3r` (electrons bohr\ :sup:`L`), with real harmonics
:math:`R_{L0} = Y_{L0}`, :math:`R_{LM} = \sqrt2 (-1)^M \mathrm{Re}\, Y_{LM}` (:math:`M>0`),
:math:`R_{LM} = \sqrt2 (-1)^M \mathrm{Im}\, Y_{L|M|}` (:math:`M<0`), so that

.. math::

   Q_{LM} = \sum_{ab} n_{ab}\,\mathcal R_{l_a l_b L}\int b_a \overline{b_b} R_{LM}\, d\Omega .

The radial functions are expanded to second order in :math:`E - E_l` about the energy :math:`E_l` of the partial wave
of shell :math:`l` (``EVAL`` of ``NEWRHO``, shifted by ``vmad`` to the Hamiltonian energy scale),

.. math::

   \psi_l(E) = \phi + (E - E_l)\dot\phi + \tfrac12 (E - E_l)^2 \ddot\phi .

With the raw energy moments :math:`M_p = \int f(E) E^p n\, dE` (:math:`p = 0,1,2`, about :math:`E = 0`) and
:math:`n^{(ij)}_{ab} = \int (E - E_a)^i (E - E_b)^j n_{ab}\, dE` (:math:`E_a` of the shell of :math:`a`):
:math:`n^{(00)} = M_0`, :math:`n^{(10)} = M_1 - E_a M_0`, :math:`n^{(01)} = M_1 - E_b M_0`,
:math:`n^{(11)} = M_2 - (E_a + E_b) M_1 + E_a E_b M_0`, :math:`n^{(20)} = M_2 - 2E_a M_1 + E_a^2 M_0`,
:math:`n^{(02)} = M_2 - 2E_b M_1 + E_b^2 M_0`. The radial pair function is

.. math::

   n^{(00)} \phi\phi' + n^{(10)} \dot\phi\phi' + n^{(01)} \phi\dot\phi' + n^{(11)} \dot\phi\dot\phi'
   + \tfrac12\left(n^{(20)} \ddot\phi\phi' + n^{(02)} \phi\ddot\phi'\right),

where each product of two functions is the large component weighted with :math:`\sqrt{\mathrm{GFAC}_a \mathrm{GFAC}_b}` plus the
small component. For :math:`a, b` in the same shell it is :math:`Q_0\phi^2 + 2Q_1\phi\dot\phi + Q_2(\dot\phi^2 + \phi\ddot\phi)`,
the spherical density of ``NEWRHO``. Approximations: second order in :math:`E - E_l` (third and higher orders neglected);
the radial functions are spin-resolved (up and down blocks each with their own :math:`\phi`) and the spin off-diagonal blocks
of :math:`n` do not contribute to the charge density; the cross-shell (:math:`l_a \ne l_b`) weight
:math:`\sqrt{\mathrm{GFAC}_a \mathrm{GFAC}_b}` of the large component is a choice, not derived.

Known mismatches
----------------

Details and measurements are in ``tests/KNOWN_ISSUES.md``, entry "Orbital basis
audit".

* ``paoflow2rs`` takes its file in the spherical basis without a transform.
* ``disp_matrix`` (``symbolic_atom.f90``) is in cubic order for l <= 2 and in
  real m order for f; ``calculate_jijk`` converts it for l <= 2 and stops with
  a fatal for an f shell.
* Export files written before the rename carry ``..._l_then_m_increasing`` in
  the ``orbital_order`` metadata for the layout given above (m descending); the
  effective importer accepts both strings, the writers emit
  ``..._l_then_m_decreasing``.
