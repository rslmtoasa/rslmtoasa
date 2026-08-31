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
- Output verbosity
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

**Default:** 1

**Example:**

.. code-block:: fortran

   nsp = 1  ! Scalar relativistic, collinear global mode

**Related code:** ``source/control.f90``, ``source/self.f90::process()``

**See also:** :ref:`theory/lmto_asa_overview` for theory

llsp
----

**Type:** Integer

**Purpose:** Recursion cutoff for s and p electrons

**Typical range:** 50-150

**Default:** 80

**Meaning:**

- Number of Lanczos vectors generated during recursion
- Controls energy resolution and accuracy
- Higher value = better accuracy but slower

**Example:**

.. code-block:: fortran

   llsp = 100  ! More accurate

**Notes:**

- Too low -> missed spectral features
- Too high -> diminishing returns (convergence slows down)
- Convergence study recommended (see :doc:`../user_guide/examples`)

**Related code:** ``source/recursion.f90::recur()``

**See also:** :ref:`keywords/lattice_geometry`

lld
---

**Type:** Integer

**Purpose:** Recursion cutoff for d electrons

**Typical range:** 50-150

**Default:** 80

**Notes:**

- Often same as ``llsp``, but can be different
- For d-rich systems (Fe, Co, Ni), may need higher values
- Independent of ``llsp`` allows flexibility

**Example:**

.. code-block:: fortran

   lld = 100
   llsp = 60   ! Different for s/p and d electrons

**Related code:** ``source/recursion.f90``

max_iterations
--------------

**Type:** Integer

**Purpose:** Maximum number of SCF iterations

**Typical range:** 30-150

**Default:** 100

**Meaning:**

- If SCF doesn't converge by this iteration, calculation stops
- Set higher for difficult systems
- Set lower for quick preliminary runs

**Example:**

.. code-block:: fortran

   max_iterations = 50  ! Quick calculation

**Notes:**

- Typical convergence: 10-30 iterations for well-behaved systems
- Very stubborn systems: 100-200 iterations needed
- Monitor ``scf_convergence.dat`` to see if trend is toward convergence

**Related code:** ``source/self.f90::process()``

dq_tol
------

**Type:** Real (Ry)

**Purpose:** SCF convergence criterion (charge density)

**Typical range:** 1e-6 to 1e-4

**Default:** 1e-5

**Meaning:**

.. math::

   \sum_i |q_i^{(n+1)} - q_i^{(n)}| < \text{dq\_tol}

- Calculation converges when charge density change is below this value
- Smaller value = better convergence but may take more iterations

**Example:**

.. code-block:: fortran

   dq_tol = 1.0e-5  ! Reasonable accuracy

**Notes:**

- For production calculations: 1e-5 or better
- For quick tests: can relax to 1e-4
- If very small (1e-7+), may have numerical issues

**Related code:** ``source/self.f90::converge_scf()``

verbose
-------

**Type:** Logical

**Purpose:** Enable verbose output during calculation

**Allowed values:** ``.true.`` or ``.false.``

**Default:** .false.

**Example:**

.. code-block:: fortran

   verbose = .true.

**Effect:**

- Prints state of structure, control, self-consistent, charge after setup
- Useful for debugging parameter values
- Increases output volume (check disk space for long runs)

**Related code:** ``source/calculation.f90::pre_processing_*()``

random_vec_num
--------------

**Type:** Integer

**Purpose:** Number of stochastic vectors for Chebyshev moment calculation

**Typical range:** 1-100

**Default:** 10

**Meaning:**

When computing Chebyshev moments for DOS, use random vector estimation:

.. math::

   \mu_n \approx \frac{1}{N_{\text{random}}} \sum_i \langle r_i | H^n | r_i \rangle

- Higher value = more accurate but slower
- For small systems: can use 1 (deterministic)
- For large systems: use 10-50 (statistical average)

**Example:**

.. code-block:: fortran

   random_vec_num = 20  ! Accurate

**Related code:** ``source/recursion.f90::chebyshev_recur_ll()``

**See also:** :ref:`theory/recursion_method`

mext
----

**Type:** Integer

**Purpose:** Acceleration of spin rotation (for magnetic calculations)

**Allowed values:** 0 (off), 1 (linear), 2 (Broyden)

**Default:** 0

**Meaning:**

For non-collinear calculations (nsp=3 or 4), accelerate convergence of spin directions:

- 0: No acceleration
- 1: Linear extrapolation
- 2: Broyden acceleration

**Example:**

.. code-block:: fortran

   mext = 2  ! Broyden acceleration for spins

**Notes:**

- Only relevant for non-collinear magnetism
- Can significantly speed up convergence
- Similar to density mixing but for spin directions

**Related code:** ``source/control.f90``, ``source/self.f90``

lrot
----

**Type:** Logical

**Purpose:** Rotate spins to local coordinate system

**Default:** .false.

**Notes:**

- Advanced feature (beta, not recommended)
- When true, rotates basis to local spin quantization axis
- Usually leave as false

**Example:**

.. code-block:: fortran

   lrot = .false.  ! Recommended

incorb
------

**Type:** Logical

**Purpose:** Include orbital moment in local spin axis determination

**Default:** .false.

**Notes:**

- For non-collinear calculations with spin-orbit coupling
- When true: magnetic axis = spin + orbital angular momentum
- Usually false for spin-only considerations

idos
----

**Type:** Integer

**Purpose:** Output type for local density of states (LDOS)

**Allowed values:**

.. list-table::
   :widths: 10 20 60
   :header-rows: 1

   * - idos
     - Output
     - Description
   * - 0
     - None
     - No LDOS output
   * - 1
     - Atom-averaged
     - LDOS s,p,d for first atom type
   * - 2
     - Per-atom
     - LDOS s,p,d for each atom type

**Default:** 0

**Example:**

.. code-block:: fortran

   idos = 2  ! Output LDOS per atom type

**Related code:** ``source/density_of_states.f90``

nlim
----

**Type:** Integer

**Purpose:** Cluster size control (limiting cluster vs. bulk)

**Allowed values:** >= 0

**Default:** 0 (bulk)

**Meaning:**

- 0: All first neighbors included (bulk-like)
- n > 0: Limit cluster to atoms within n nearest-neighbor shells
- Used to study cluster vs. bulk behavior

**Notes:**

- For bulk calculations: use 0
- For cluster studies: increase to isolate finite-size effects

**Related code:** ``source/lattice.f90``

npold
-----

**Type:** Integer

**Purpose:** Number of old densities to store (mixing history)

**Typical range:** 5-20

**Default:** 10

**Notes:**

- Related to ``broyden_history`` in SCF mixing
- Allows Broyden mixing to use multiple previous iterations
- More history -> better convergence but more memory

orb_pol
-------

**Type:** Logical

**Purpose:** Include orbital polarization corrections

**Default:** .false.

**Notes:**

- Advanced feature for improved magnetism description
- Adds mean-field term to enhance orbital anisotropy
- Usually false for standard calculations

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
  up to the cutoff ``llsp``/``lld``.  Preferred for small clusters and
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
- An unrecognized value is non-fatal (warning + CPU fallback).

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

**Default:** ``'B'``

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
  ``bands.f90``, and ``energy.f90``; validated in
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

**Purpose:** Activate the atomistic spin dynamics (ASD) module after the SCF
calculation.

**Example:**

.. code-block:: fortran

   do_asd = .true.

**Notes:**

- ASD time-integration parameters are set in the ``&sd`` namelist
  (see :ref:`keywords/spin_dynamics`).

**See also:** :ref:`theory/spin_dynamics`, :ref:`keywords/spin_dynamics`

asd_jij
-------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Include computed exchange parameters :math:`J_{ij}` as the
effective field driving the ASD simulation.

**Notes:**

- Only meaningful when ``do_asd = .true.``.

cond_type
---------

**Type:** Character string

**Purpose:** Selects the formulation of the conductivity tensor calculation

**Allowed values:** ``'kubo'``, ``'kubo-bastin'``

**Default:** ``'kubo-bastin'``

**Related code:** ``source/conductivity.f90``

**See also:** :ref:`reference/conductivity_module`

cond_calctype
-------------

**Type:** Character string

**Purpose:** Selects the sub-type (orbital decomposition, trace method) used
inside the conductivity calculation

**Notes:**

- Allows selection between full tensor, diagonal, or off-diagonal evaluation.

cond_ll
-------

**Type:** Integer

**Purpose:** Chebyshev expansion cutoff used inside the conductivity module,
independent of the SCF recursion cutoffs ``llsp``/``lld``

**Typical range:** 300–1000

**Default:** 500

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

**Purpose:** Use spherically averaged atomic spheres (spherical vacancy
approximation) for the potential within the Wigner-Seitz spheres.

blockrec
--------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Enable block recursion, processing multiple starting vectors
simultaneously for improved parallel efficiency.

partype
-------

**Type:** Integer

**Default:** ``0``

**Purpose:** Selects the source of tight-binding potential parameters.
0 uses parameters from the ``&par`` namelist (standard).

terminator
----------

**Type:** Integer

**Default:** ``0``

**Purpose:** Selects the terminator function applied at the end of the recursion
chain to close the continued-fraction Green's function

**Allowed values:**

- ``0``: Beer-Pettifor terminator (analytic band-edge approximation, recommended)
- ``1``: Luchini-Nex terminator
- ``2``: Constant terminator

**Related code:** ``source/recursion.f90``

conca / concb
-------------

**Type:** Real

**Default:** 1.0 / 0.0

**Purpose:** Concentrations of species A and B for binary alloy calculations

**Example:**

.. code-block:: fortran

   conca = 0.75
   concb = 0.25

**Notes:**

- Only relevant when ``ruban > 0``.
- Should satisfy ``conca + concb = 1``.

ruban
-----

**Type:** Real

**Default:** ``0.0``

**Purpose:** Ruban-Abrikosov concentration parameter; a positive value activates
CPA-like alloy averaging.

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
  ``docs/dev/plans/B6_surface_electrostatics.md`` for what shipped vs.
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

**Default:** ``'constrained_spiral'``

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
- ``'relaxed_reference'`` -- the full rotating-frame Cartesian moment is mixed
  and the per-sublattice reference axis follows it, while the single-q ansatz
  is retained. The transverse residual and torque vanish identically.

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

**Default:** ``.false.``

**Purpose:** Enable co-charge (``do_cochg``) and co-moment (``do_comom``)
corrections for CPA alloy disorder.

Provenance
==========

Control parameters defined and used in:

- **Definition:** ``source/control.f90::type control``
- **Reading:** ``source/control.f90::build_from_file()``
- **Usage:**
  
  - Calculation driver: ``source/calculation.f90::process()``
  - SCF loop: ``source/self.f90::process()``
  - Recursion: ``source/recursion.f90::recur()``
  - Output: ``source/density_of_states.f90``

See Also
========

- :ref:`keywords/scf_settings` - Related SCF parameters (mixing, dq_tol)
- :ref:`keywords/energy_mesh` - Related energy parameters
- :ref:`theory/scf_cycle` - Understanding SCF convergence
- :ref:`theory/recursion_method` - Understanding recursion cutoffs
- :doc:`../user_guide/input_files` - Input file syntax
