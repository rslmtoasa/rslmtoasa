.. _keywords/control_parameters:

=======================================
Control Parameters (&control Namelist)
=======================================

Overview
========

The ``&control`` namelist contains the main calculation settings, including:

- Type of calculation (relativistic treatment, collinearity, etc.)
- SCF loop limits and tolerances
- Recursion cutoff depths
- Pre/post-processing options

Parameters
==========

nsp
---

**Type:** Integer

**Purpose:** Specifies calculation type (collinearity and relativistic treatment)

This is the global electronic-structure mode, not the number of spin-density
channels passed to an XC evaluator. The atomic radial/XC path keeps two local
spin-density channels even when ``nsp=1``.

**Allowed values:**

.. list-table::
   :widths: 10 30 60
   :header-rows: 1

   * - nsp
     - Type
     - Description
   * - 1
     - Scalar relativistic (SR) collinear
     - One Hamiltonian, fixed z-axis
   * - 2
     - Fully relativistic (FR) collinear
     - Spin-orbit coupling, z-axis
   * - 3
     - SR non-collinear
     - Spins arbitrary direction
   * - 4
     - FR non-collinear
     - Spin-orbit + arbitrary spin direction

**Default:** none. ``nsp`` is not initialised in ``control%restore_to_default``; give it in every input.

**Example:**

.. code-block:: fortran

   nsp = 1  ! Scalar relativistic, collinear global mode

**Related code:** ``source/control.f90``, ``source/self.f90::run()``

**See also:** :ref:`theory/lmto_asa_overview` for theory

llsp
----

**Type:** Integer

**Purpose:** Enters only the allocation size of the recursion coefficient
arrays, ``max(llsp, lld)``. The recursion depth is set by ``lld``.

**Default:** 16

**Notes:**

- No other part of the code reads ``llsp``; it does not change the recursion
  depth or the energy resolution.

**Example:**

.. code-block:: fortran

   llsp = 16

**Related code:** ``source/recursion_core.f90``

**See also:** :ref:`keywords/lattice_geometry`

lld
---

**Type:** Integer

**Purpose:** Recursion cutoff: number of recursion levels used by the
recursion and Chebyshev routines (``llsp`` only sets an array size)

**Typical range:** 50-150

**Default:** 16

**Notes:**

- Higher values resolve more spectral detail at higher cost
- For d-rich systems (Fe, Co, Ni), may need higher values

**Example:**

.. code-block:: fortran

   lld = 100

**Related code:** ``source/recursion.f90``

Iteration limit and convergence
-------------------------------

The SCF iteration limit and convergence threshold are ``nstep`` and ``conv_thr`` in the ``&self`` namelist (not in
``&control``); see :ref:`keywords/scf_settings`.

.. note::

   ``verbose`` is a key of the ``&calculation`` namelist, not of ``&control``
   (a ``verbose`` in ``&control`` stops the run); see
   :ref:`keywords/output_options`.

random_vec_num
--------------

**Type:** Integer

**Purpose:** Number of random vectors of the stochastic trace in the
conductivity (transport) moments

**Typical range:** 1-100

**Default:** 1

**Meaning:**

Used only with ``cond_calctype = 'random_vec'``. The trace of the transport
moments is estimated with random vectors:

.. math::

   \mu_n \approx \frac{1}{N_{\text{random}}} \sum_i \langle r_i | H^n | r_i \rangle

- Higher value = more accurate but slower
- The DOS moments do not use this key.

**Example:**

.. code-block:: fortran

   random_vec_num = 20  ! Accurate

**Related code:** ``source/recursion_transport.f90``, ``source/conductivity.f90``

**See also:** :ref:`theory/recursion_method`

mext
----

**Type:** Integer

**Default:** 0

Read but currently not used. The magnetic moments are mixed with
``mix_magnetic_moments`` and ``magbeta`` of ``&mix``.

lrot
----

**Type:** Logical

**Default:** .false.

Read but currently not used.

incorb
------

**Type:** Logical

**Default:** .false.

Read but currently not used.

idos
----

**Type:** Integer

**Default:** 0

Read but currently not used. The real-space DOS route writes ``totaldos.out``
and ``<symbol>_dos.out`` (``source/bands.f90``) independently of this key.

nlim
----

**Type:** Integer

**Default:** none (not set in ``control%restore_to_default``)

Read but currently not used.

npold
-----

**Type:** Integer

**Purpose:** Orbital dimension of the structure-constant arrays; the code uses
``max(norb, npold)``. Rarely set by hand.

**Default:** ``norb`` (orbitals per spin channel of the basis: 4, 9 or 16)

.. note::

   ``orb_pol`` is a key of the ``&hamiltonian`` namelist, not of ``&control``
   (an ``orb_pol`` in ``&control`` stops the run); see
   :ref:`keywords/hamiltonian_parameters`.

recur
-----

**Type:** Character string

**Purpose:** Selects the recursion algorithm used to compute the Green's function

**Allowed values:** ``'lanczos'``, ``'chebyshev'``, ``'block'``

**Default:** ``'block'``

**Example:**

.. code-block:: fortran

   recur = 'chebyshev'

**Notes:**

- ``'lanczos'``: Scalar Lanczos recursion — tridiagonalises the Hamiltonian
  up to the cutoff ``lld``.  Preferred for small clusters and
  impurity calculations.
- ``'block'``: Block-Lanczos recursion, processing several starting vectors
  together (fast kernels in ``haydock_fast.f90``).  This is the default
  production path.
- ``'chebyshev'``: Expands the Green's function in Chebyshev polynomials with
  stochastic trace evaluation (controlled by ``random_vec_num``).  More
  efficient for large disordered systems; fast kernels in
  ``chebyshev_fast.f90``, optionally offloaded to the GPU (see ``gpu_plugin``).

**Related code:** ``source/recursion.f90`` (+ ``recursion_core`` /
``recursion_haydock`` / ``recursion_chebyshev`` submodules)

**See also:** :ref:`theory/recursion_method`, :ref:`theory/gpu_acceleration`

gpu_plugin
----------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Offload the supported recursion kernels to an NVIDIA GPU via the
CUDA plugin.

**Example:**

.. code-block:: fortran

   gpu_plugin = .true.

**Notes:**

- Requires a build configured with ``ENABLE_CUDA_PLUGIN=ON`` and an available
  CUDA device.  When the plugin is not compiled in or no device is present,
  the run falls back to the CPU path (this keyword is then ignored).
- Covers Chebyshev, block-Lanczos, scalar Lanczos (``nsp=1`` only), stochastic,
  intersite and orbital-moment kernels — bypassing the CPU fast backends.

**Related code:** ``source/rsrec_cuda_plugin.f90``, ``source/recursion_core.f90``

**See also:** :ref:`theory/gpu_acceleration`

gpu_backend
-----------

**Type:** Character string

**Default:** ``'csr'``

**Purpose:** Selects the on-device sparse storage / matvec strategy used by the
CUDA plugin (only relevant when ``gpu_plugin = .true.``).

**Allowed values:** ``'csr'``, ``'bsr'``, ``'fft'``, ``'conv'``

**Notes:**

- ``'csr'`` / ``'bsr'``: compressed / block sparse row matvec (``'bsr'`` is
  usually faster for the dense per-atom orbital blocks).
- ``'fft'`` / ``'conv'``: periodic operator apply; **require a fully periodic
  lattice** (``pbc=.true.`` with ``b1``/``b2``/``b3``), otherwise the run
  warns and falls back to the CPU path.
- An unrecognized value is fatal (``control%check_all``); the same holds for ``gpu_precision``.

**Related code:** ``source/rsrec_cuda_plugin.f90::decode_gpu_backend()``

**See also:** :ref:`theory/gpu_acceleration`

.. note::

   The combined-correction (``ccor_2c``) and H-O-H (``hoh``) Hamiltonian
   switches live in the ``&hamiltonian`` namelist, not ``&control`` — see
   :ref:`keywords/hamiltonian_parameters`.

calctype
--------

**Type:** Character string

**Purpose:** Selects the cluster/embedding topology the run is built on —
**not** which post-SCF property is computed (that is ``post_processing``,
see :ref:`keywords/output_options`).

**Allowed values:** ``'B'`` (bulk), ``'S'`` (surface), ``'I'`` (impurity),
``'L'`` (layered/interface — ``buildinterface``, B7)

**Default:** none (``''``); the run stops unless one of the four values is given.

**Example:**

.. code-block:: fortran

   calctype = 'L'

**Notes:**

- ``'B'``: ordinary bulk/supercell cluster (``bravais`` pre-processing).
- ``'S'``: one-sided surface/slab cluster (``buildsurf``).
- ``'I'``: an impurity cluster embedded in a bulk or surface host
  (``newclubulk``/``newclusurf``).
- ``'L'``: interface/layered cluster with two frozen reference regions
  (``buildinterface``, B7 — see :doc:`../user_guide/examples/interface_fcccu111`
  for ``region_b_kind='metal'``/``'vacuum'`` and the alignment machinery).
- Consumed throughout ``lattice_cluster.f90``, ``calculation.f90``,
  ``bands.f90``, and ``energy.f90``; checked in
  ``control.f90::check_all`` (fatal on any other value).
- This entry previously listed post-processing values (``'bulk'``,
  ``'surface'``, ``'bands'``, ``'exchange'``, ``'conductivity'``, ``'dos'``)
  which do not correspond to ``calctype`` at all — those are
  ``post_processing`` values; see :ref:`keywords/output_options`.

**Related code:** ``source/control.f90::check_all()``, ``source/lattice_cluster.f90``

do_asd
------

**Type:** Logical

**Default:** ``.false.``

Read but currently not used. The ASD module is started by
``processing = 'sd'`` in ``&calculation``; its parameters are set in the
``&sd`` namelist (see :ref:`keywords/spin_dynamics`).

**Notes:**

- ``source/calculation.f90`` selects the ASD run from ``processing``, not
  from ``do_asd``.

**See also:** :ref:`theory/spin_dynamics`, :ref:`keywords/spin_dynamics`

asd_jij
-------

**Type:** Logical

**Default:** ``.false.``

Read but currently not used.

linear_in
---------

**Type:** Character string

**Purpose:** Operator of the linear-response input (the perturbation)

**Allowed values:** the same operator names as ``linear_out`` (below).
For ``linear_in /= 'charge'`` the conductivity output is written in raw
units (no S/cm conversion).

**Default:** ``'charge'``

**Related code:** ``source/recursion_transport.f90``, ``source/conductivity.f90``

linear_out
----------

**Type:** Character string

**Purpose:** Operator of the linear-response output (the measured response)

**Allowed values:** ``'charge'``, ``'spin'``, ``'orbital'``, ``'spin_accumulation'``,
``'orbital_accumulation'``, ``'spin_torque'``, ``'orbital_torque'``,
``'spin_soc_torque'`` (alias ``'soc_spin_torque'``), ``'quadrupole'``, ``'oap'``,
``'quadrupole_accumulation'``, ``'oap_accumulation'``, ``'octupole'``,
``'octupole_current'``, ``'octupole_accumulation'``; any other value is fatal

**Default:** ``'charge'``

**Related code:** ``source/conductivity.f90``

**See also:** :ref:`reference/conductivity_module`

intersite_mode
--------------

**Type:** Character string

**Purpose:** How the inter-site Green functions ``G_ij`` and ``G_ji`` of the pair list ``ijpair`` are obtained from the recursion (``recur = 'chebyshev'`` or ``'block'``).

**Allowed values:** ``'pairs'``, ``'one_shot'``

**Default:** ``'pairs'``

**Notes:**

- ``'pairs'``: four recursions per pair, started from :math:`(|i\rangle + c|j\rangle)/\sqrt{2}` with :math:`c = \pm 1, \pm i`; self-pairs use one recursion.
- ``'one_shot'``: one recursion per distinct first atom ``i`` of the pair list, started from site ``i`` alone. The blocks at every requested target ``j`` are read off during that recursion, and ``G_ij`` follows from ``G_ji`` by Hermitian conjugation. Pairs that share the first atom share one recursion, so list pairs grouped by ``i``. See :ref:`theory/recursion_method` (section "Inter-site Green functions: pairs and one-shot").
- With ``recur = 'chebyshev'`` both routes compute the same moments; measured differences in ``J_ij`` on bcc Fe were 4e-13 to 6e-13 relative. Works with ``cheb_lp`` and ``gpu_plugin`` (fp32 and fp64).
- With ``recur = 'block'`` the two routes build different Krylov spaces at finite ``lld``, so off-site ``G_ij``, ``J_ij`` and ``D_ij`` differ from ``'pairs'``; the difference was about 2 % for the first-neighbour ``J_ij`` on the bcc Fe contour test case and can be much larger for other pairs and depths. Self-pairs agree with ``'pairs'`` to rounding. The terminator (``sym_term``) is the same as for ``'pairs'``. Works with ``post_processing = 'exchange'`` and ``'exchange_p2rs'``, and with ``gpu_plugin`` (fp32 and fp64; the small projected solve runs on the CPU).

**Related code:** ``source/intersite_request.f90``, ``source/recursion_transport.f90`` (``chebyshev_recur_ij_oneshot``, ``recur_b_ij_oneshot``), ``source/green_block.f90`` (``block_green_ij_oneshot``), ``source/block_tridiag.f90``, ``source/calculation_reciprocal.f90``

exchange_integration
--------------------

**Type:** Character string

**Purpose:** How the energy integral of the exchange observables (``J_ij``, ``D_ij``, ``A_ij`` and the two-index ``so``/``fo`` parts) is evaluated on the real-space recursion route (``post_processing = 'exchange'``, ``gf_route = 'recursion'``).

**Allowed values:** ``'dos'``, ``'moments'``

**Default:** ``'dos'``

**Notes:**

- ``'dos'``: ``G_ij`` is reconstructed on the energy mesh and integrated with the Simpson rule.
- ``'moments'``: the integral over :math:`[E_{min}, E_F]` (Fermi function for ``temperature > 0``) is evaluated from the Chebyshev moments of ``G_ij`` and ``G_ji``. Requires ``recur = 'chebyshev'`` and ``gf_route = 'recursion'``. ``do_damping`` and ``do_inertia`` use the series at :math:`E_F` and do not write the energy-resolved files. See :ref:`theory/recursion_method` (section "Inter-site Green functions: pairs and one-shot").

**Related code:** ``source/green_chebyshev.f90``, ``source/exchange_dynamics.f90``

scf_integration
---------------

**Type:** Character string

**Purpose:** How the on-site energy integrals of the real-space SCF (charges, spin and orbital moments, band energy, density matrices) are evaluated.

**Allowed values:** ``'dos'``, ``'moments'``

**Default:** ``'dos'``

**Notes:**

- ``'dos'``: ``G`` is reconstructed on the energy mesh and integrated with the Simpson rule.
- ``'moments'``: ``E_F`` solves :math:`N(E_F)` = valence from the on-site Chebyshev moments (bisection, 1e-12 Ry), ``G`` is reconstructed at Gauss-Legendre nodes of :math:`[E_{min}, E_F]` and the integrals are weighted sums over the nodes (Fermi function at the ``&energy`` ``temperature``; 0 gives the step). Requires ``recur = 'chebyshev'``. The energy-resolved DOS and cumulative-moment files are not written. See :ref:`theory/recursion_method`.

**Related code:** ``source/bands.f90`` (``calculate_fermi``, ``integrate``), ``source/green_chebyshev.f90`` (``chebyshev_scf_nodes``)

multipole_route
---------------

**Type:** Character string

**Purpose:** How the on-site multipole moments of the sphere charge are obtained.

**Allowed values:** ``'q10'``, ``'general'``

**Default:** ``'q10'``

**Notes:**

- ``'q10'``: today's route; the l = 1 dipole :math:`Q_{10}` for ``dipole_electrostatics`` from the on-site Green function (``compute_dipole_moments``). Only :math:`n^{(00)}` and the spin-averaged large-component :math:`\phi` of each site enter, with the elements s-p(m = 0), p(m = 0)-d(m = 0) and p(m = +-1)-d(m = +-1); the radial functions are those of ``rad_fun`` of the site. Bit-identical to a run without the keyword.
- ``'general'``: all :math:`Q_{LM}`, :math:`L = 0..2 l_{max}`, from the raw energy moments :math:`M_p = \int f(E) E^p n\,dE`, :math:`p = 0..2`, of the on-site density matrix and the radial functions of ``NEWRHO`` (second order in :math:`E - E_\nu`). Requires ``scf_integration = 'moments'`` (fatal otherwise). Runs every SCF step from the second one (the radial functions are those of the previous atomic step), also without ``dipole_electrostatics``.
- Each SCF step writes ``multipoles.out`` (overwritten): one header line, then ``site L M Q_LM`` per line, :math:`Q_{LM}` in electrons*bohr\ :sup:`L` for the real harmonics :math:`R_{LM}` (M > 0 cos, M < 0 sin).
- With ``dipole_electrostatics = .true.``, ``potential%q10`` takes :math:`\sqrt{4\pi/3}\,Q_{10}/S` (S the Wigner-Seitz radius), mixed with ``dipole_mix`` as in the ``'q10'`` route. See :ref:`theory/orbital_basis`.

**Related code:** ``source/electrostatics_multipole.f90`` (``compute_multipoles``), ``source/green_chebyshev.f90`` (``chebyshev_onsite_energy_moments``)

nonlocal_dm
-----------

**Type:** Logical

**Purpose:** Run the Chebyshev recursion of the real-space SCF as the one-shot inter-site recursion, so that one recursion per recursion site gives the on-site moments the SCF needs and the inter-site moments :math:`P_j^\dagger T_m P_i` to every neighbour :math:`j` of the site.

**Default:** ``.false.``

**Notes:**

- Requires ``recur = 'chebyshev'`` and ``scf_integration = 'moments'`` (fatal otherwise).
- The on-site moments are computed in double precision with the Hamiltonian action of the one-shot route, whatever ``cheb_backend`` is, and cost about :math:`2\times` the product-identity recursion.
- Off, the SCF is bit-identical to a run without the keyword.
- Each SCF step, after :math:`E_F` is set, writes ``nonlocal_dm.out`` (overwritten): one header line, then for every recursion site ``i`` and each pair ``(i, j)`` (self pair first, then the structure-constant neighbours) a line ``site_i atom_j Rx Ry Rz`` (:math:`R` = ``cr(:,j) - cr(:,i)`` in units of ``alat``) followed by ``nb`` lines of ``nb`` (Re, Im) values of the density matrix block :math:`n_{ij}` (electrons).
- :math:`n_{ij} = -\frac{1}{\pi}\sum_m \mathrm{Im}(c_m)\, M_{ji}(m)^\dagger`, with :math:`M_{ji}(m) = \langle j|T_m(\tilde H)|i\rangle` the one-shot moments and :math:`c_m = \int f(E)\, g_m(E)\, dE` the linear weights of the moment-space SCF (same window, kernel, :math:`E_F` and temperature).

**Related code:** ``source/nonlocal_density.f90``, ``source/recursion_chebyshev.f90``, ``source/recursion_transport.f90``

nonlocal_dm_order
-----------------

**Type:** Integer

**Purpose:** Highest raw energy moment :math:`p` of the non-local density matrix that is computed and written: :math:`n_{ij}^{(p)} = \int f(E)\,E^p\,[\ldots]\,dE = -\frac{1}{\pi}\sum_m \mathrm{Im}(c^{(p)}_m)\, M_{ji}(m)^\dagger`, with :math:`c^{(p)}_m = \int f(E)\,E^p\, g_m(E)\, dE`.

**Allowed values:** ``0``, ``1``, ``2``

**Default:** ``0``

**Notes:**

- Fatal outside 0..2, and fatal if larger than 0 without ``nonlocal_dm = .true.``.
- The moments are taken about :math:`E = 0` (not about :math:`E_\nu`), in electrons\*Ry\ :sup:`p`.
- ``p = 0`` is ``nonlocal_dm.out`` (unchanged); ``p = 1`` and ``p = 2`` are written to ``nonlocal_dm_e1.out`` and ``nonlocal_dm_e2.out`` in the same format, with the moment named in the header line.

**Related code:** ``source/nonlocal_density.f90``, ``source/bands.f90``

cheb_backend
------------

**Type:** Character string

**Purpose:** CPU engine of the Chebyshev DOS moments and, with the Haydock block recursion, of the fast kernels.

**Allowed values:** ``'fast_dp'``, ``'fast'``, ``'legacy'``, ``'batched'``, ``'mkl_batch'``, ``'mkl_sparse'``

**Default:** ``'fast_dp'``

**Notes:**

- ``'fast_dp'`` and ``'legacy'``: double-precision Chebyshev moments (the matrix-free loops) and double-precision
  DOS reconstruction. ``'fast_dp'`` also selects the fp64 block-Lanczos kernel.
- ``'fast'`` and the other fast backends: single-precision moments and reconstruction (about 1e-7 relative);
  the single-precision rounding breaks exact symmetries such as zero transverse moments.
- The transport moments use ``transport_backend`` instead.

cheb_lp
-------

**Type:** Logical

**Purpose:** Extend the Chebyshev moments by linear prediction before the Jackson-damped reconstruction.

**Default:** ``.false.``

**Notes:**

- Requires ``recur = 'chebyshev'`` (fatal otherwise). The moments of each recursion atom and, for the exchange route, of each pair slot (four per pair) are extended from ``2*lld+2`` to ``(2*lld+2)*cheb_lp_factor`` orders by ``x(n) = -sum_i a_i x(n-i)``, after Ganahl et al., PRB 90, 045144. One real coefficient vector per block is shared by the real and imaginary parts of all ``nb*nb`` elements.
- Order ``p = N/2`` trained on all ``N`` known moments; pseudo-inverse cutoff 1e-6; roots with ``|lambda| >= 1`` are set to zero. These constants are module parameters in ``source/cheb_lp.f90``, not input keys.
- The moments of all local atoms (``chebyshev_recur``) and pair slots (``chebyshev_recur_ij``) are extended in one OpenMP loop after the moment loops; each block is independent and runs on one thread, so the ``Chebyshev LP`` log lines may appear in any order. The CUDA path calls the same loop after the GPU moments of all atoms and slots are back on the host.
- BLAS threading matters. Measured on ``example/exchange/bccFe_chebyshev_lp`` (12 pair slots, factor 2; 4 cores shared, gfortran 13, OpenBLAS pthread build), time from the first to the last ``Chebyshev LP pair`` log line: 1.54 s with 1 OpenMP thread and 4 OpenBLAS threads, 1.09 s with 1 and 1, 1.22 s (2.5 s in an earlier, more loaded run) with 4 and 4, 0.33 s with 4 OpenMP threads and ``OPENBLAS_NUM_THREADS=1``; ``J_ij`` was identical to the digits written in all cases. With the current kernel (normal equations from the lag Gram matrix, one real matrix product for the extension, Hermitian reduction) the same runs gave 0.25 s (1 and 1), 0.08 s (4 and 1) and 0.65 s (4 and 4). When running with ``cheb_lp`` and several OpenMP threads, limit the BLAS library to one thread (``OPENBLAS_NUM_THREADS=1`` or ``MKL_NUM_THREADS=1``) for this step; there is no code switch for it.
- Cost, measured in a separate run: one extension call at ``nb`` = 18 (648 real channels), 302 known moments, factor 2 took about 105 ms with the first implementation (72 ms in a quiet window) and 25 ms with the current one (one BLAS thread, shared machine), which agrees with it to 3e-11 relative on the extended moments of bcc Fe (the first implementation itself changes by 4e-12 under a 1e-16 input perturbation). That is about 0.1% of an on-site run (moment recursion about 22 s per atom) and well below 1% per slot on the triad exchange deck. The Hermitian reduction assumes the moment block is Hermitian and falls back to all channels when it is not (checked to a relative 1e-8).
- ``mu_n`` takes 16 bytes per element, ``nb*nb*N_ext`` per slot (about 3.1 MB at ``nb`` = 18, ``lld`` = 150, factor 2, i.e. 604 orders); the pair route has ``4*njij`` slots. There is no second array of kernel-weighted moments (the former ``mu_ng``); the two routes that need them build a scratch per call, and the DOS reconstruction does not copy ``mu_n``. The log gives one line ``Chebyshev moments: <n> orders x <slots> slots; mu_n = <x> MB`` and, per atom or pair slot, ``Chebyshev LP ...: p= rank= unstable=<n> (weight <w>) lam_max= fit_res= imag=``. ``weight`` is the share of the mode amplitude at the last known moment carried by the dropped roots.
- The spectral positivity projection still runs and its count is logged (``local_blocks=``) when non-zero. With ``cheb_lp`` it is followed by a rescale of each atom, orbital and spin block so that its integrated weight equals the weight before the clipping (the spin coherence by the geometric mean of the two factors); the largest ``|1-s|`` is logged as ``Chebyshev LP charge-conserving rescale``. Without the rescale the clipping adds charge, because the extended moments fix the zeroth moment, and the Fermi level moves: measured on frozen bcc Fe, E_F was 4.8 mRy from the computed reference (factor 4, clipping only) against 0.56 mRy (factor 2, rescaled; computed 604-moment run), with the integrated charge 18.0027 against 18.0026; on the Si example deck (lld=200, factor 2) E_F was 0.0656 Ry, below the gap (0.0738 to 0.2911 Ry), without the rescale and 0.205 Ry, inside it, with the rescale (lld=400 run: 0.1727 Ry). Not covered: ``chebyshev_dos`` in ``density_of_states.f90``, the orbital moments (``mu_n_orb``) and the conductivity moments. The CUDA hooks compile but were not run on a GPU.
- Limit (measured, bcc Fe, nsp=2): moments of successive open-boundary clusters agree to 1e-3 only up to n = 88 (5984 sites), 110 (9720) and 136 (16888); LP does not remove the boundary contamination beyond that. LP trained only on boundary-free moments did not predict a larger system's moments better than the smaller system's own moments (0 of 288 configurations). LP approximates a run with more moments on the same cluster, not a larger cluster.
- Factor 2 (default), measured: bulk bcc Fe (5984 sites, ``lld`` = 150, frozen potential, N_ext = 604), RMS DOS difference to a computed 604-moment Jackson run, occupied band and within +-0.1 Ry of E_F, up/down spin: 0.14/0.13/0.14/0.14 states/Ry, against 0.89/1.44/0.47/0.75 for the plain 302-moment run; factor 4 against its computed 1208-moment run: 0.69/0.83/0.54/0.68. Ni(001) surface deck ``tests/ldau/decks/ni001_soc_u3`` (LDA+U, E_F fixed at -0.076843 Ry, ``lld`` = 200, factor 2, reference ``lld`` = 401), occupied-band RMS DOS difference of plain ``lld`` = 200 and LP, up/down: Ni-S 0.746/0.866 and 0.104/0.120; Ni-S-1 0.614/0.475 and 0.120/0.143; empty sphere 0.014/0.025 and 0.002/0.002 states/Ry. DOS(E_F) for ``lld`` = 200, LP, reference: Ni-S 17.354, 17.017, 16.983; Ni-S-1 12.318, 12.064, 12.057 states/Ry. Si deck (``tests/scf``, ``real_space_chebyshev``, ``lld`` = 200, factor 2 against ``lld`` = 400): largest DOS difference 1.48 states/Ry with LP, 5.10 for plain ``lld`` = 200.
- The Jackson kernel is applied at the extended length (kernel on). Ganahl et al. do not apply it to their predicted moments, but the predicted moments here decay exponentially. Reconstruction of the extended bulk Fe moments without the kernel gave 473 negative DOS points and RMS errors 1.2/1.4/0.8/1.2 states/Ry, against 0.14/0.13/0.14/0.14 with Jackson at factor 2.
- One deck (bcc Fe bulk, 5984 sites, ``lld`` = 150, factor 4, ``channels_ldos`` = 2500) against a computed ``lld`` = 603 run on the same cluster: RMS DOS difference over the occupied band, up/down spin, 1.59/0.87 states/Ry (Jackson at 302 moments) and 0.67/0.52 (LP); within 0.1 Ry of E_F, 2.72/1.31 and 0.83/0.68; E_F -0.0683 (302 moments), -0.0743 (LP), -0.0695 (reference) Ry; 6480 local blocks projected with LP, 0 without; wall time 23.6 s (LP), 23.0 s, 77.8 s (reference). Exchange deck (1000-site periodic cluster, three shells): ``J_ij`` error against the ``lld`` = 603 reference +0.229/+0.013/-0.010 (302 moments) and +0.0206/+0.0074/-0.0074 (LP). ``cheb_backend = 'fast'`` against fp64: largest DOS difference 0.05 (302 moments) and 0.49 (LP) states/Ry for a peak of about 40. Two MPI ranks gave the same ``J_ij`` as the serial run.

**Related code:** ``source/cheb_lp.f90``, ``source/recursion_chebyshev.f90``, ``source/recursion_transport.f90``, ``source/green_chebyshev.f90``; campaign scripts archived at the tag ``cheb-lp-campaign-2026-11``; examples ``example/bulk/bccFe_chebyshev_lp``, ``example/exchange/bccFe_chebyshev_lp``

cheb_lp_factor
--------------

**Type:** Integer

**Purpose:** Extension factor of the linear prediction: the stored moment length is ``(2*lld+2)*cheb_lp_factor``.

**Default:** 2

**Allowed values:** 2 or larger (fatal otherwise)

**Notes:**

- Ignored unless ``cheb_lp`` is true.

**Related code:** ``source/cheb_lp.f90``, ``source/control.f90``

cheb_kernel
-----------

**Type:** Character string

**Purpose:** Damping kernel g_n applied to the Chebyshev moments in the DOS and Green-function reconstruction (CPU routes and the GPU plugin; the GPU path is implemented but has not been run on a GPU).

**Allowed values:** ``'jackson'``, ``'lorentz'``, ``'none'``

**Default:** ``'jackson'``

**Notes:**

- ``'jackson'`` is the kernel of all earlier versions; the default results are unchanged. ``'lorentz'`` uses ``g_n = sinh(lambda (1 - n/N)) / sinh(lambda)`` with ``cheb_kernel_lambda``. ``'none'`` uses ``g_n = 1``.
- ``'none'`` requires ``cheb_lp = .true.``: a hard cutoff of computed moments rings. Ganahl et al., PRB 90, 045144, apply no kernel to their linear-prediction moments; they predict 10^4 to 10^5 moments of a smooth spectrum whose added moments decay exponentially. The moments of the cluster Green functions here do not decay that way (bcc Fe: diagonal moments of order 0.1 at n = 1200).
- Measured (cheb_lp, factor 2, ``lld`` = 200 against plain Jackson at ``lld`` = 400, Si example deck ``tests/scf``, real_space_chebyshev): RMS DOS difference 0.31 with Jackson, 0.59 with Lorentz (lambda = 3), 1.00 with no kernel, 0.85 for plain ``lld`` = 200. On bcc Fe (factor 2, computed 604-moment reference, no positivity projection) the no-kernel reconstruction had 473 negative-DOS points and RMS errors of 1.0, 2.1, 0.4 and 0.6 states/Ry against 0.14, 0.13, 0.14, 0.14 with Jackson. A small imaginary shift (1e-5, as in the paper) did not change this. Si deck at factor 4 (``lld`` = 200, reference plain ``lld`` = 803): RMS DOS difference 1.28 (Jackson), 1.39 (Lorentz, lambda = 3), 1.36 (no kernel), i.e. the three are within 0.11 of each other; no kernel choice improved on the others for these moments.
- Requires ``recur = 'chebyshev'``. With ``gpu_plugin = .true.`` the damping vector is built on the host and passed to the CUDA reconstruction (the real-axis DOS and the intersite eta route); without it the CUDA code builds the Jackson factor itself, so ``'jackson'`` takes the previous code path. This GPU path is implemented but not yet run on a GPU: neither the Jackson bit-identity nor the fp32 behaviour of the other kernels has been measured there. The orbital moments (``mu_n_orb``) and the conductivity moments keep their own kernels.

**Related code:** ``source/chebyshev_fast.f90`` (optional ``damping`` argument of ``cheb_green_fast``, ``cheb_green_dp``, ``cheb_green_complex``), ``source/green_chebyshev.f90`` (``cheb_damping``), ``source/rsrec_cuda_plugin.f90`` and ``source/cuda/rsrec_gpu.cu`` (optional ``kernel`` argument of ``rsrec_chebyshev_dos``), ``source/math.f90`` (``jackson_kernel``, ``lorentz_kernel``)

cheb_kernel_lambda
------------------

**Type:** Real

**Default:** 3.0

**Purpose:** Parameter lambda of the Lorentz kernel (``cheb_kernel = 'lorentz'``). Must be positive; ignored for the other kernels.

transport_backend
-----------------

**Type:** Character string

**Purpose:** Engine that computes the stochastic Chebyshev moments of the
conductivity (``linear_out``/``linear_in`` transport). Independent of ``cheb_backend``.

**Allowed values:**

- ``'legacy'``: matrix-free double-precision loops (all operators, ``hoh``,
  ``ccor_2c``, impurity sites).
- ``'fast_sp'``: single-precision CPU kernel with bond and operator blocks indexed by atom type.
  Runs every ``linear_in``/``linear_out`` operator, with or without ``hoh``, except ``spin_torque`` and
  ``orbital_torque`` with ``hoh``. The other cases (``ccor_2c``, impurity sites, bond vectors that differ
  between sites of one atom type) run ``'legacy'`` and write one warning to the log.
- ``'gpu'``: reserved; fatal at run time.

**Default:** ``'legacy'``

**Related code:** ``source/recursion_transport.f90``, ``source/chebyshev_fast.f90``

cond_calctype
-------------

**Type:** Character string

**Purpose:** Selects how the trace of the stochastic transport moments is taken

**Allowed values:** ``'per_type'``, ``'random_vec'``

**Default:** ``'per_type'``

**Notes:**

- ``'per_type'``: one start vector set per atom type.
- ``'random_vec'``: ``random_vec_num`` random vectors.
- The code does not check the value; use one of the two.

cond_ll
-------

**Type:** Integer

**Purpose:** Chebyshev expansion cutoff used inside the conductivity module,
independent of the SCF recursion cutoffs ``llsp``/``lld``

**Typical range:** 300–1000

**Default:** 200

**Example:**

.. code-block:: fortran

   cond_ll = 600

**Notes:**

- Higher values give a smoother conductivity spectrum but increase cost.

**See also:** :ref:`reference/conductivity_module`

svac
----

**Type:** Logical

**Default:** ``.false.``

Read but currently not used.

blockrec
--------

**Type:** Logical

**Default:** ``.false.``

Read but currently not used. The recursion type is selected by ``recur``.

partype
-------

**Type:** Integer

**Default:** ``1``

Read but currently not used.

terminator
----------

**Type:** Integer

**Default:** ``5``

Read but currently not used.

conca / concb
-------------

**Type:** Real

**Default:** 0.0 / 0.0

Read but currently not used.

ruban
-----

**Type:** Real

**Default:** ``0.0``

Read but currently not used.

dipole_electrostatics (B6)
--------------------------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Activate l=1 (dipole) surface/interface electrostatics —
extends the monopole-only ASA Madelung problem with a per-atom dipole
moment Q₁₀, fed into the potential shift via the existing (previously
unused) dipole-monopole Madelung matrix ``dsz``.

**Example:**

.. code-block:: fortran

   dipole_electrostatics = .true.

**Notes:**

- With this ``.false.`` (the default), behavior is bit-identical to the
  code without the feature — this is the regression contract for B6.
- No new lattice-sum machinery: reuses the ``dsz``/``dzz`` Madelung
  matrices that were already computed but previously unfed. Narrower in
  scope than the original B6 blueprint (no 2D/3D Ewald to l≤2) — see
  ``docs/DECISIONS.md`` for what shipped vs.
  what was originally planned.
- No literature (Skriver–Rosengaard) work-function validation has been
  performed yet — gate G-B6-1 is open.
- The computed moment is exported per atom as ``potential%q10`` (see
  the ``&potential`` output namelist).

**Related code:** ``source/electrostatics_multipole.f90``
(``compute_dipole_moments``), ``source/charge.f90`` (``surfpot``)

dipole_mix (B6)
---------------

**Type:** Real

**Default:** ``0.5``

**Purpose:** SCF mixing fraction for the dipole moment Q₁₀ between
iterations, analogous to charge mixing (``mixing``/``alpha``): ``q10_new
= dipole_mix * q10_computed + (1 - dipole_mix) * q10_old``.

**Example:**

.. code-block:: fortran

   dipole_mix = 0.3

**Notes:**

- Only relevant when ``dipole_electrostatics = .true.``.

**Related code:** ``source/electrostatics_multipole.f90``

density_policy (WP7)
--------------------

**Type:** String

**Default:** ``'relaxed_reference'``

**Allowed values:** ``'constrained_spiral'``, ``'relaxed_reference'``

**Purpose:** Selects how the shared rotating-frame spin density
(``source/spin_density.f90``) is reduced to the variables the SCF mixes, and
what the radial up/down projection axis means. Both solvers -- the real-space
Green-function route and the k-space eigenvector route -- accumulate the full
per-site, per-l 2x2 spin matrix first and project onto radial channels only
afterwards, against the axis this policy resolves.

- ``'constrained_spiral'`` -- the per-site reference direction is imposed and
  held fixed. Only the charge and the *longitudinal* moment magnitude
  ``m . ref`` are mixed; the transverse density ``m - (m . ref) ref`` and the
  torque ``ref x m`` are reported per site as the constraint residual
  (``DENSITY_POLICY`` log lines) rather than absorbed.
- ``'relaxed_reference'`` -- the moment direction is mixed with ``&mix`` ``magbeta``,
  ``mom = (1 - magbeta) mom_old + magbeta mom_measured``, as in the main branch;
  ``magbeta = 0`` keeps the axes fixed. The real-space density is projected on
  that mixed ``mom``, which is not replaced by the measured axis.
  ``magbeta`` is ignored by ``'constrained_spiral'``. ``gbt_single_q`` needs a
  fixed axis and switches to ``'constrained_spiral'`` with a warning.

**Example:**

.. code-block:: fortran

   density_policy = 'relaxed_reference'

**Notes:**

- On a single-sublattice collinear ferromagnet the two policies coincide, as
  they must: the accumulated moment is already parallel to the reference.
- The projection axis is always stated explicitly; a stale ``potential%mom``
  is never used as an implicit projection definition.

**Related code:** ``source/spin_density.f90``, ``source/bands.f90``
(``accumulate_spin_density_rs``, ``resolve_density_axes``),
``source/reciprocal_spin_density.f90``

gbt_scf_diagnostics (WP12)
--------------------------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Emit ``gbt_scf_diagnostics.dat`` with one row per site and SCF
iteration. The record joins the incoming frame/moment and radial channels,
the eigenstate density-matrix components and moment, constraint target/field,
magnetic mixer inputs/outputs, next radial channels and ``potential%mom``,
Hamiltonian exchange markers, physical DFT energy, and separate constraint
diagnostics. In ``gbt_single_q`` the density-matrix components and moments are
reported in the primitive rotating frame; ordinary SCF rows use the active
solver frame.

This is an audit trace and is off by default because it performs file I/O.
It does not change the SCF equations.

**Related code:** ``source/self.f90``

hyperfine
---------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Compute hyperfine coupling constants at nuclear sites after SCF.

sym_term
--------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Apply symmetry-based terminator corrections to the Green's function
recursion.

do_cochg / do_comom
-------------------

**Type:** Logical

**Default:** ``.false.`` / ``.false.``

``do_cochg``: read but currently not used.

``do_comom``: read once, in the radial solver of ``source/self.f90``. When it
is ``.true.`` and the atom is not forced non-magnetic, the field term
``B_fsm`` of the radial equation is set to ``-mag_cfield(3)``, minus the
z component of the constraining field of the atom.

Provenance
==========

Control parameters defined and used in:

- **Definition:** ``source/control.f90::type control``
- **Reading:** ``source/control.f90::build_from_file()``
- **Usage:**
  
  - Calculation driver: ``source/calculation.f90::process()``
  - SCF loop: ``source/self.f90::run()``
  - Recursion: ``source/recursion.f90::recur()``
  - Output: ``source/density_of_states.f90``

Transport and GPU options (&control, &conductivity)
===================================================

.. list-table::
   :header-rows: 1
   :widths: 34 10 18 38

   * - Key
     - Type
     - Default
     - Meaning
   * - ``cpu_reconstruction_precision``
     - character
     - ``'fp64'``
     - Precision of the CPU reconstruction of the transport response: ``'fp64'`` or
       ``'fp32'``.
   * - ``gpu_precision``
     - character
     - ``'fp32'``
     - Arithmetic precision of the CUDA Chebyshev kernels: ``'fp32'`` or ``'fp64'``. The
       reconstructed output stays double precision.
   * - ``fresh_start``
     - logical
     - ``.false.``
     - Start from the named elemental files even when a ``<element>_out.nml`` restart
       is present.
   * - ``cpp_plugin``, ``gpu_moment_download``, ``gpu_stochastic_block``
     - logical, logical, integer
     - ``.false.``, ``.false.``, 1
     - Read but not used.
   * - ``response_geometry``
     - character
     - ``'bulk'``
     - ``&conductivity``: normalise the response per volume (``'bulk'``) or per area
       (``'sheet'``).
   * - ``response_measure``
     - real (Å³ or Å²)
     - 0.0
     - ``&conductivity``: cell volume or area used for the normalisation; 0 derives it
       from the lattice vectors and ``alat``.

See Also
========

- :ref:`keywords/scf_settings` - Related SCF parameters (nstep, conv_thr, mixtype, beta)
- :ref:`keywords/energy_mesh` - Related energy parameters
- :ref:`theory/scf_cycle` - Understanding SCF convergence
- :ref:`theory/recursion_method` - Understanding recursion cutoffs
- :doc:`../user_guide/input_files` - Input file syntax

response_projection
-------------------

The default 'seed' closes the response trace on the representative site's
basis. Select 'local' to place an onsite output observable on that site.
Local projection uses deterministic per_type, charge input and real-space
recursion. It accepts spin, orbital, quadrupole (including OAP), and octupole
accumulations and spin_torque_onsite. Use pol_alpha for the spin component
and v_beta for the electric-field direction.
