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

   nsp = 1  ! Scalar relativistic, collinear (ferromagnetic)

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

- Too low → missed spectral features
- Too high → diminishing returns (convergence slows down)
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
- More history → better convergence but more memory

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

**Allowed values:** ``'lanczos'``, ``'chebyshev'``

**Default:** ``'lanczos'``

**Example:**

.. code-block:: fortran

   recur = 'chebyshev'

**Notes:**

- ``'lanczos'``: Tridiagonalises the Hamiltonian exactly up to the cutoff
  ``llsp``/``lld``.  Preferred for small clusters and impurity calculations.
- ``'chebyshev'``: Expands the Green's function in Chebyshev polynomials with
  stochastic trace evaluation (controlled by ``random_vec_num``).  More
  efficient for large disordered systems.

**Related code:** ``source/recursion.f90``

**See also:** :ref:`theory/recursion_method`

calctype
--------

**Type:** Character string

**Purpose:** Selects the top-level calculation to perform after SCF convergence

**Allowed values:** ``'bulk'``, ``'surface'``, ``'bands'``, ``'exchange'``,
``'conductivity'``, ``'dos'``

**Default:** ``'bulk'``

**Example:**

.. code-block:: fortran

   calctype = 'exchange'

**Notes:**

- Determines which post-SCF module is invoked.
- ``'exchange'``: compute Heisenberg exchange parameters :math:`J_{ij}`.
- ``'conductivity'``: compute the Kubo-Bastin conductivity tensor.
- ``'bands'``: compute band structure along a k-path.

**Related code:** ``source/calculation.f90``

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
