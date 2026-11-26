.. _keywords/output_options:

=======================================
Output Control Parameters
=======================================

Overview
========

Several parameters control what output files are generated and their verbosity level.

Parameters
==========

verbose
-------

**Type:** Logical

**Purpose:** Enable verbose output during calculation

**Default:** .false.

**Example:**

.. code-block:: fortran

   verbose = .true.

**Effect:**

- Prints detailed information after each setup phase
- Calls ``print_state()`` routines for structure, control, etc.
- Useful for debugging or understanding parameter values
- Increases screen output significantly

**Where to find:**

- In the ``&calculation`` namelist only (an invalid key in ``&control``)

**Related code:** ``source/calculation.f90::pre_processing_*()``

idos
----

**Type:** Integer

**Default:** 0

Read but currently not used (``&control`` key). The real-space SCF step writes
``totaldos.out`` (two columns: energy minus Fermi energy in Ry, total DOS) and
``<symbol>_dos.out`` per recursion atom (``source/bands.f90``) independently of
this key.

post_processing
---------------

**Type:** Character string

**Purpose:** Type of post-processing/property calculation

**Allowed values:**

.. list-table::
   :header-rows: 1

   * - post_processing
     - Action
   * - ``'none'``
     - No post-processing (default)
   * - ``'hamiltonian_export'``
     - Export the frozen native Hamiltonian
   * - ``'dos_frozen'``
     - Calculate frozen native real-space DOS
   * - ``'dos_import'``
     - Calculate DOS from a real-space Hamiltonian export
   * - ``'conductivity_import'``
     - Run the Kubo linear response with an imported real-space Hamiltonian
   * - ``'exchange_import'``
     - Calculate pair exchange, DMI and tensors using an imported Hamiltonian and the native pair list
   * - ``'density_of_states'``
     - Calculate native reciprocal DOS
   * - ``'band_structure'``
     - Calculate band structure
   * - ``'paoflow2rs'``
     - PAOFLOW-to-RS conversion workflow
   * - ``'exchange'``
     - Calculate exchange interactions
   * - ``'conductivity'``
     - Calculate transport properties
   * - ``'exchange_p2rs'``
     - Exchange in p2rs workflow
   * - ``'conductivity_p2rs'``
     - Conductivity in p2rs workflow
   * - ``'orbital_modern'``
     - Orbital-moment post-processing
   * - ``'bsf'``
     - Bloch spectral function A(k,E) (B3)
   * - ``'kspace_green'``
     - k-space Green's-function driver (B2)
   * - ``'frozen_magnon'``
     - Frozen-magnon E(q)/omega(q) sweep (B1)
   * - ``'fermi_surface'``
     - Dense k-space eigenpair export for Fermi-surface analysis
   * - ``'spin_response'``
     - Spin-response driver (see :ref:`keywords/spin_response`)

**Default:** 'none'

**Example:**

For standalone fixed Hamiltonian workflows, set ``pre_processing='none'``
and ``processing='none'``.
``post_processing='hamiltonian_export'`` builds the frozen native operator
and applies the ``&hamiltonian`` export settings.
``'dos_frozen'`` evaluates its real-space DOS with Chebyshev or block recursion.
The three import routes read the exported real-space Hamiltonian and its
metadata. Conductivity requires Chebyshev recursion. DOS and pair exchange
also support block recursion. ``'dos_effective'`` and
``'conductivity_effective'`` remain accepted aliases.
``'conductivity_p2rs'`` reads a standalone external ``paoham.dat``.
See the :download:`Hamiltonian export and import manual <../../HAMILTONIAN_EXPORT.md>`
for workflow selection and examples.

.. code-block:: fortran

   post_processing = 'exchange'

**Where to find:**

- In ``&calculation`` namelist

**Related code:** ``source/calculation.f90::process()`` - selects post-processing routine

**Output files generated:**

- 'density_of_states': ``dos_kspace.dat`` (k-space DOS). The real-space SCF
  step writes ``totaldos.out`` and ``<symbol>_dos.out`` (``source/bands.f90``).
- 'band_structure': ``band_structure.dat``
- 'exchange': ``jijso.out``, ``jijfo.out``, ``dijso.out``, ``dijfo.out``,
  ``aijso.out``, ``aijfo.out``, ``jtensso.out``, ``jtensfo.out`` and the
  ``jijparts.out``, ``dijparts.out``, ``aijparts.out`` files; ``jij.out``,
  ``dij.out``, ``aij.out``, ``jtens.out`` in the other exchange variants
  (``source/exchange.f90``)
- 'conductivity': ``cond_total.out``, ``cond_total_orb_real.out``,
  ``cond_total_orb_im.out``, ``response_total.out``,
  ``response_orbitals_real.out``, ``response_orbitals_imag.out``
  (``source/conductivity.f90``)
- 'bsf': ``bsf.dat`` and ``bsf_bands.dat`` (see
  ``source/reciprocal_bsf.f90::calculate_bsf()``)
- 'frozen_magnon': the file named by ``&frozen_magnon output_file``
  (default ``frozen_magnon.dat``), plus ``frozen_magnon_branches.dat``/
  ``frozen_magnon_modes.dat`` when ``branch_mode='auto'`` — see
  :ref:`keywords/frozen_magnon`.

gf_route
--------

**Type:** Character string

**Purpose:** Selects which engine fills the ``green%gij``/``gji`` (and
``gij_eta`` torque-ladder) arrays consumed by exchange, damping, and
conductivity post-processing (B2/B5) — a route-agnostic dispatch key so
those consumers do not need to know which engine produced the Green's
function.

**Allowed values:** ``'recursion'`` (real-space Haydock/Chebyshev, the
original route), ``'lehmann'`` (k-space, strict Lehmann representation,
:math:`\Sigma=0`, B2 backend E), ``'dyson'`` (k-space, direct Dyson
inversion :math:`[zS-H(k)-\Sigma]^{-1}`, B2 backend D)

**Default:** ``'recursion'``

**Example:**

.. code-block:: fortran

   gf_route = 'lehmann'

**Notes:**

- In ``&calculation`` namelist.
- At ``gf_route='recursion'`` this is a pure no-op relative to pre-B2
  behavior — bit-identical output.
- At ``'lehmann'``/``'dyson'`` the exchange/damping/conductivity pipeline
  is unchanged; only the Green's-function producer differs. See
  ``green_backend``/``green_eta`` (:ref:`keywords/reciprocal`)
  for the k-space engine's own settings, and
  ``docs/DECISIONS.md`` for the accepted convergence
  envelopes (Lehmann-route J_ij/damping are broadening-defined at the
  documented ``green_eta`` default, not identical to the recursion route).

**Related code:** ``source/calculation.f90::check_gf_route()``

exchange_delta
~~~~~~~~~~~~~~

Selects the exchange splitting for ``post_processing='exchange_import'``.
``'lmto'`` (default) uses the energy-dependent LMTO potential matrix.
``'onsite'`` uses the full real part of the onsite up-minus-down
Hamiltonian block for a collinear ``spd`` model, following the PAOFLOW
pair-exchange convention.

do_damping
----------

**Type:** Logical

**Purpose:** Enable Gilbert-damping (:math:`\alpha`) evaluation in the
exchange post-processing pipeline (B5.3), via the SOC-derivative
torque-correlation route already present in the ``gij_eta`` ladder.

**Default:** ``.false.``

**Example:**

.. code-block:: fortran

   do_damping = .true.

**Notes:**

- In ``&calculation`` namelist.
- Runs on whichever Green's function ``gf_route`` filled — recursion,
  Lehmann, or Dyson — with no separate wiring per route.

**Related code:** ``source/calculation.f90`` (search ``do_damping``)

do_chi0
-------

**Type:** Logical

**Purpose:** In the exchange post-processing, also evaluate the static transverse
bare susceptibility :math:`\chi^0_{ij}(\omega = 0)` of the d channel for every pair of
the pair list and write it to ``chi0_ij.dat``:

.. math::

   \chi^0_{ij} = -\frac{X - \bar Y}{2\pi i}, \quad
   X = \int f(E)\, \mathrm{Tr}_d[G^{\uparrow\uparrow}_{ji} G^{\downarrow\downarrow}_{ij}]\, dE, \quad
   Y = \int f(E)\, \mathrm{Tr}_d[G^{\downarrow\downarrow}_{ji} G^{\uparrow\uparrow}_{ij}]\, dE .

Same sign and units (1/Ry) as ``spin_response`` with ``method = 'mills'`` and
``channels = 'd'``. Without SOC :math:`X = Y` and :math:`\chi^0_{ij}` is real.

**Default:** ``.false.``

**Notes:**

- In ``&calculation`` namelist. Needs ``nsp = 1``, or ``nsp = 2`` with every
  moment along :math:`\pm z`; other runs stop.
- Runs on whichever ``gf_route`` filled ``G_ij``, and with both values of
  ``exchange_integration``. The on-site term :math:`\chi^0_{ii}` needs the pair
  :math:`(i, i)` in the pair list.
- With SOC the transverse channel alone is not closed: the coupling to the
  longitudinal and charge channels is not included.
- Also writes ``chi0_types.dat``, one row per atom type that has a self pair
  :math:`(i, i)`; the first such :math:`i` of the type is its seed. Columns: the
  d moment :math:`m_d = -\tfrac{1}{\pi}\mathrm{Im}\int f\,\mathrm{Tr}_d[G^{\uparrow\uparrow}_{ii} - G^{\downarrow\downarrow}_{ii}]\,dE`
  of the seed, :math:`U` from the Goldstone condition
  :math:`\sum_b \chi^0_{ab}(q=0)\, U_b m_b = -m_a` (holds by construction),
  :math:`U_\mathrm{Mills} = (C^\downarrow_d - C^\uparrow_d)/m_d` from ``center_band``, the
  Mills residual :math:`1 + \sum_b \chi^0_{ab} U^\mathrm{Mills}_b m_b / m_a`, and
  :math:`\chi^0_{ab}(q=0)`, the sum of :math:`\mathrm{Re}\,\chi^0_{ij}` over the listed pairs of
  the seed with :math:`j` of type :math:`b`. The sum runs over the pair list as given:
  list every atom of the shells to include, not one per shell. The three :math:`U`
  columns are ``n/a`` unless every type has a seed and a moment of at least
  1e-6 :math:`\mu_B`. The definitions follow ``spin_response`` with ``method = 'mills'``.
- In moment space (``exchange_integration = 'moments'``) :math:`m_d` uses the linear
  Chebyshev weights on the same interval, window and kernel as the bilinear ones.
- Jülich vertex (as ``spin_response`` with ``method = 'juelich'``): the d amplitudes are
  weighted by :math:`F(e) = \alpha + \beta e = (1 + (e - E_\nu)(E_F - E_\nu)p)/\sqrt{1 + (E_F - E_\nu)^2 p}`,
  :math:`p = 1/\mathrm{ppar}^2`, at the energy of each state, i.e. :math:`G` is replaced by
  :math:`\tilde G_{ab} = \langle a|F_a(H)\,G\,F_b(H)|b\rangle = F_a(z)F_b(z)G_{ab} - \beta_a\beta_b H_{ab}
  - (\alpha_a\beta_b + \beta_a\alpha_b + \beta_a\beta_b z)\delta_{ab}`. Multiplying :math:`G` by
  :math:`F(z)` alone is not equivalent. :math:`H_{ab}` comes from the first Chebyshev moment, so the
  Jülich columns of ``chi0_ij.dat`` and ``chi0_types.dat`` are filled on the real-space Chebyshev
  route only (``recur = 'chebyshev'``, ``gf_route = 'recursion'``) and are ``n/a`` otherwise. In moment
  space :math:`\tilde G` is expanded in its own moments :math:`\langle a|F_a T_m F_b|b\rangle` (two fewer
  than computed); on the energy mesh the identity above is applied to the damped :math:`G`, so the two
  integration modes differ at the level of the expansion error.
- Also writes ``chi0_adiabatic.dat``, the adiabatic spectrum at the q points of the
  ``&spin_response`` namelist (``n_q``, ``q_list``, ``q_file``; q = 0 when the namelist is absent):
  :math:`\chi^0_{ab}(\mathbf q) = \sum \chi^0_{ij}\, e^{2\pi i\, \mathbf q\cdot \mathbf r_{ij}}` over the listed
  pairs of the seed of type :math:`a` with :math:`j` of type :math:`b` (:math:`\mathbf r_{ij}` in direct
  coordinates; with every pair included this is the :math:`\chi^0_{ab}(\mathbf q)` of ``spin_response``), and
  :math:`\Omega(\mathbf q) = -N(\chi^0(\mathbf q)^{-1} + U)` with :math:`N = \mathrm{diag}(m)` and the Goldstone
  :math:`U`, so :math:`\Omega(0) = 0` by construction. Columns as ``<prefix>_adiabatic.dat`` of
  ``spin_response`` (frequencies, largest imaginary part, force-theorem frequencies), first with the
  coefficient vertex, then with the Jülich vertex (``n/a`` where its :math:`U` is not defined). The
  lattice sum is truncated to the pair list: list every atom within the radius to include.
- Also writes ``chi0_cluster.dat``, the adiabatic spectrum of a cluster of magnetic atoms (Lounis-type
  cluster RPA). The sites :math:`S` are the atoms :math:`i` with a self pair :math:`(i, i)` in the pair list
  (atoms, not types), in pair-list order; :math:`\chi^0_S(a, b)` is :math:`\chi^0` of the listed pair
  :math:`(S_a, S_b)`, or the complex conjugate of :math:`(S_b, S_a)` if only that one is listed (the static
  :math:`\chi^{0,+-}` matrix is Hermitian, :math:`\chi^0_{ji} = \overline{\chi^0_{ij}}`), or 0 if neither is
  listed (truncation to the pair list). :math:`m_a` is the d moment of the self pair of :math:`S_a`, :math:`U_a` the
  Goldstone :math:`U` of :math:`\mathrm{Re}\,\chi^0_S` and :math:`m` over the cluster, and
  :math:`\Omega = -N(\chi^{0\,-1}_S + U)` with :math:`N = \mathrm{diag}(m)`, :math:`|S|` frequencies, the
  largest :math:`|\mathrm{Im}|` and the force-theorem frequencies, for the coefficient and the Jülich vertex.
  The file has three blocks, each after its own header line: one row per site (m, :math:`U`), one row per
  frequency index, one row with the largest imaginary parts. The values of a vertex are ``n/a`` when there
  is no site, a site has :math:`|m| < 10^{-6}\,\mu_B`, or (Jülich) the vertex is absent. The result is meaningful
  when every magnetic atom that couples is in :math:`S` (adatoms and clusters on non-magnetic hosts); for
  magnetic hosts the cluster Goldstone :math:`U` and the truncation make it approximate.
- With ``do_exchange`` also writes ``j0_sum_rule.dat``, a local check of the exchange integral. For a
  first-order Hamiltonian :math:`H^s = C^s + W^s S W^s` with spin-independent :math:`S`
  (``nsp = 1``, ``hoh = .false.``), :math:`\sum_j \mathrm{Tr}[\delta_i G^\uparrow_{ij}\delta_j G^\downarrow_{ji}]
  = \sum_k \delta_k (r_k G^\downarrow_{kk} - G^\uparrow_{kk}/r_k)`, with :math:`\delta_k(E)` the exchange
  vertex of ``d_matrix`` and :math:`r_k = W_\downarrow/W_\uparrow`, :math:`j` over the whole cluster and
  :math:`k` over the orbitals :math:`l \le 2` of atom :math:`i`; only the self pair :math:`G_{ii}` enters. The
  file lists this local :math:`J_0` next to the sum of the listed :math:`J_{ij}` of atom :math:`i`. The identity
  is only approximate with ``hoh``, with spin-orbit coupling, for a pair list shorter than the cluster and for
  a kernel-damped :math:`G`. Not written with ``exchange_delta = 'onsite'``.

- Dynamic :math:`\chi^0_{ij}(\omega)` on the energy-mesh route (coefficient vertex), from the Lehmann form with
  :math:`A = (G - G^\dagger)/(-2\pi i)` and :math:`G^A_{ji} = G_{ij}^\dagger`:

  .. math::

     \chi^0_{ij}(\omega) = \int f(E)\,\mathrm{Tr}_d[A^{\uparrow\uparrow}_{ji}(E)\, G^{\downarrow\downarrow}_{ij}(E+\omega)]\,dE
     + \int f(E)\,\mathrm{Tr}_d[G^{\uparrow\uparrow,A}_{ji}(E-\omega)\, A^{\downarrow\downarrow}_{ij}(E)]\,dE .

  At :math:`\omega = 0` this is algebraically the static expression above. Poles sit at
  :math:`\omega = e^\downarrow - e^\uparrow` (positive for a majority-up moment). The energy mesh must be uniform
  (otherwise a warning and no dynamic output); each :math:`\omega` of the ``&spin_response`` grid (``omega_min``,
  ``omega_max``, ``n_omega``; defaults 0, 0.05 Ry, 201, negative values allowed) is rounded to the nearest multiple
  :math:`k\,\Delta E` of the mesh step and the rounded value is what is written, so the resolution is :math:`\Delta E`.
  A term whose shifted energy :math:`E \pm \omega` lies outside the mesh is dropped. ``eta`` of ``&spin_response`` is not
  used: the broadening is that of the recursion (its terminator or Chebyshev kernel). Not available in moment space
  (``exchange_integration = 'moments'``: one info line, no dynamic files). Cost: one pass over the mesh per pair and
  :math:`\omega`. Writes ``chi0_dynamic.dat`` and ``chi_rs_q.dat``: for each q of ``&spin_response`` the
  :math:`\chi^0_{ab}(\mathbf q, \omega)` lattice sum at every :math:`\omega`, :math:`\chi = (I + \chi^0 K)^{-1}\chi^0` with
  :math:`K = \mathrm{diag}(U)` and the Goldstone :math:`U` of ``chi0_types.dat``, and
  :math:`\mathrm{tr}\,L = -\mathrm{Im}\,\mathrm{tr}\,\chi/\pi`. Coefficient vertex only (no Jülich). At q = 0 and
  :math:`\omega = 0` the matrix :math:`I + \chi^0 K` is singular (Goldstone) and :math:`\chi` is ``n/a``. For a self pair
  :math:`\int d\omega\,(-1/\pi)\,\mathrm{Im}\,\chi^0_{ii}(\omega)` over all :math:`\omega` equals :math:`m_d`.

**Related code:** ``source/spin_response_rs.f90``, ``source/exchange_dynamics.f90``
(``chi0_pair_mesh``, ``chi0_pair_moments``, ``j0_local_mesh``, ``j0_local_moments``, ``chi0_pair_dynamic``)

do_sigma
--------

**Type:** Logical

**Purpose:** After a converged SCF (``pre_processing = 'bravais'``, real-space Chebyshev), evaluate once the local Green
function :math:`G_{ii}(z) = \langle i|(z - H - \Sigma_t(z))^{-1}|i\rangle` of every recursion atom for a site-local,
non-Hermitian self-energy :math:`\Sigma_t` per atom type :math:`t`, read from ``sigma.dat``, and write it to ``gloc_sigma.dat``.

``sigma.dat``: a free header line; ``nz ntype nb``; then for each :math:`z` a line ``x y`` (:math:`z = E_F + x + iy`, Ry)
followed by ``ntype`` blocks of ``nb`` lines, each line holding the ``nb`` entries of one row of :math:`\Sigma_t` as
``re im`` pairs (Ry, same basis as the Hamiltonian). ``gloc_sigma.dat``: header line, ``nz nrec nb``, then per :math:`z` the
line ``x y`` and one :math:`G_{ii}` block (1/Ry) per recursion atom in the same row layout; ``n/a`` lines mark a block that did not converge.

**Default:** ``.false.``

**Notes:**

- In ``&calculation`` namelist; ``sigma_tol`` (default ``1.0e-6``) is the relative tail bound that marks a block converged.
- The expansion is in :math:`(H + \Sigma - b)/a` with a complex centre :math:`b` and uses at most ``lld`` moments (order
  :math:`n` costs :math:`n` operator applications). It converges when :math:`z` lies outside the ellipse enclosing the
  numerical range of :math:`H + \Sigma`: points too close to the real axis relative to :math:`\max|\mathrm{Im}\,\Sigma|`
  are reported ``n/a``.

**Related code:** ``source/self_energy.f90``, ``source/cheb_sigma_kernel.f90``

sigma_tol
---------

**Type:** Real

**Purpose:** Relative tail bound below which a :math:`G_{ii}(z)` block of ``do_sigma`` is written; see ``do_sigma``.

**Default:** ``1.0e-6``

green_backend
-------------

**Type:** Character string

**Purpose:** Selects the k-space engine ``gf_route='lehmann'``/``'dyson'``
dispatches to when filling ``green%gij``/``gji`` (B2, ``reciprocal_green.f90``).

**Allowed values:** ``'lehmann'`` (strict Lehmann representation,
:math:`\Sigma=0`; one diagonalization per k-point, all energies and pairs
amortized), ``'dyson'`` (direct Dyson inversion
:math:`[zS-H(k)-\Sigma(z)]^{-1}` per (k, z); required once :math:`\Sigma \neq 0`,
e.g. a future CPA/DMFT provider)

**Default:** ``'lehmann'``

**Example:**

.. code-block:: fortran

   green_backend = 'dyson'

**Notes:**

- In ``&reciprocal`` namelist (not ``&calculation`` — ``gf_route`` is set
  in ``&calculation`` and internally assigns ``reciprocal%green_backend``
  to match; see ``source/calculation.f90`` lines 815 and 899).
- With :math:`\Sigma=0`, ``'dyson'`` reproduces ``'lehmann'`` to solver
  tolerance — a permanent consistency check (gate G-B2-1, see
  ``docs/DECISIONS.md``).

**Related code:** ``source/reciprocal_green.f90``, ``source/lehmann_kernel.f90``,
``source/dyson_kernel.f90``

green_eta
---------

**Type:** Real (Ry)

**Purpose:** Retarded broadening :math:`\eta` for the k-space Green's
function :math:`z = E + i\eta` (B2). Larger :math:`\eta` gives faster
k-mesh convergence of intersite quantities (J_ij, damping) at the cost of
a broadening-defined (not recursion-matching) estimator.

**Default:** ``0.02`` Ry

**Example:**

.. code-block:: fortran

   green_eta = 0.02

**Notes:**

- In ``&reciprocal`` namelist.
- ``0.02`` Ry is the accepted working point from gate **G-B2-2**
  (signed 2026-07-16): documented as giving ~1% k-mesh convergence of
  J_ij at a 16³ mesh. Supersedes an earlier ``0.01`` Ry placeholder from
  gate G-B2-1. See ``docs/validation/reciprocal_green_convergence.md`` for the convergence study
  behind this number (the gate record G-B2-2 is at tag ``archive/fable_v4b-2026-11-05``).
- The :math:`\eta \to 0`, :math:`N_k \to \infty` limit is deliberately not
  the shipped default — it is expensive and currently gated by the
  backend-E full-unreduced-BZ requirement (no symmetry reduction / k-mesh
  parallelism yet; attaches with B4, not started).

**Related code:** ``source/reciprocal_lifecycle.f90`` (``restore_to_default``)

pre_processing
---------------

**Type:** Character string

**Purpose:** Pre-processing/geometry generation

**Allowed values:**

.. table::
   :align: left

   +--------------------+---------------------------------------+
   | pre_processing     | Action                                |
   +====================+=======================================+
   | 'none'             | No pre-processing (default)           |
   +--------------------+---------------------------------------+
   | 'bravais'          | Generate bulk cluster                 |
   +--------------------+---------------------------------------+
   | 'buildsurf'        | Generate surface from bulk            |
   +--------------------+---------------------------------------+
   | 'newclusurf'       | Insert impurity into surface          |
   +--------------------+---------------------------------------+
   | 'newclubulk'       | Insert impurity into bulk             |
   +--------------------+---------------------------------------+
   | 'buildinterface'   | Generate interface (layered) cluster  |
   +--------------------+---------------------------------------+

**Default:** 'none'

**Example:**

.. code-block:: fortran

   pre_processing = 'buildsurf'

**Where to find:**

- In ``&calculation`` namelist

**Related code:** ``source/calculation.f90::pre_processing_*()``

processing
-----------

**Type:** Character string

**Purpose:** Main calculation type

**Allowed values:**

.. table::
   :align: left

   +--------------------+---------------------------------------+
   | processing         | Type                                  |
   +====================+=======================================+
   | 'none'             | SCF only (default)                    |
   +--------------------+---------------------------------------+
   | 'sd'               | Spin dynamics                         |
   +--------------------+---------------------------------------+

**Default:** 'none'

**Where to find:**

- In ``&calculation`` namelist

**Related code:** ``source/calculation.f90::processing_*()``

Typical Output File Workflow
=============================

**Example 1: Simple SCF + DOS**

.. code-block:: fortran

   &calculation
      post_processing = 'density_of_states'
      verbose = .true.
   /

Creates:

- ``dos_kspace.dat`` - k-space density of states

**Example 2: Bulk + Surface + Impurity**

.. code-block:: fortran

   &calculation
      pre_processing = 'buildsurf'
   /
   ! (then run once to build surface)

   &calculation
      pre_processing = 'newclusurf'
      post_processing = 'exchange'
   /
   ! (run again on impurity system)

**Example 3: Complete magnetic study**

.. code-block:: fortran

   &control
      nsp = 2  ! Collinear FR with SOC
   /
   &calculation
      post_processing = 'exchange'
      verbose = .true.
   /

Creates:

- ``jijso.out``, ``jijfo.out``, ... - Exchange interactions

Controlling Output Size
=======================

**Large output:**

- Long ``channels_ldos`` -> large DOS files
- Long recursion cutoff -> slower but not larger files
- Large cluster -> more exchange pairs

**Disk space:**

Typical example calculations:

- Si (5×5×5 cluster, dos): ~100 KB - 1 MB
- Fe with exchange (7×7×7): ~1-5 MB
- Multiple energy windows: scale with number of runs

**Managing output:**

.. code-block:: bash

   # Compress after calculation
   gzip totaldos.out
   gzip jijso.out
   
   # View compressed files
   gunzip -c totaldos.out.gz | head -20

Parsing Output Files
====================

**Python example:**

.. code-block:: python

   import numpy as np
   
   # Read DOS (columns: E - E_F in Ry, total DOS)
   dos = np.loadtxt('totaldos.out')
   E = dos[:, 0]
   rho = dos[:, 1]

Provenance
==========

Output control parameters found in:

- **verbose:** ``source/calculation.f90::type calculation``
- **idos:** ``source/control.f90::type control`` (not used)
- **post_processing:** ``source/calculation.f90::type calculation``
- **pre_processing:** ``source/calculation.f90::type calculation``
- **Output writers:**
  
  - ``source/bands.f90``, ``source/reciprocal_dos.f90`` (DOS files)
  - ``source/exchange.f90`` (exchange files)
  - ``source/conductivity.f90`` (conductivity files)
  - ``source/self.f90`` (``report.out``, ``minfo.out``, ``linfo.out``)

Magnetic inertia (&calculation)
===============================

.. list-table::
   :header-rows: 1
   :widths: 34 10 18 38

   * - Key
     - Type
     - Default
     - Meaning
   * - ``do_inertia``
     - logical
     - ``.false.``
     - With ``post_processing = 'exchange'`` and ``nsp = 2``, also write the raw
       torque-correlation inertia kernel per pair (``allinertias.out``, ``inertia-energy.out``,
       ``example-real.out``, ``example-imag.out``). Needs at least three energy points.

See Also
========

- :ref:`keywords/control_parameters` - Related control options
- :doc:`../user_guide/output_files` - Output file formats and interpretation
- :doc:`../user_guide/examples` - Example output configurations
