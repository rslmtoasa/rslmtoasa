.. _keywords/scf_settings:

=======================================
SCF Convergence Parameters (&self)
=======================================

Overview
========

The ``&self`` namelist controls self-consistent field (SCF) loop settings,
including the iteration limit, the convergence threshold and magnetic mixing.
The density mixing parameters ``mixtype`` and ``beta`` are in the ``&mix`` namelist.

Parameters
==========

mixtype (in ``&mix``)
---------------------

**Type:** Character string

**Purpose:** Type of density mixing algorithm

**Allowed values:** 'linear', 'broyden'

**Default:** 'linear'

**Example:**

.. code-block:: fortran

   mixtype = 'broyden'

**Linear Mixing:**

.. math::

   n^{(i+1/2)} = \alpha n_{\text{new}}^{(i+1)} + (1-\alpha) n_{\text{old}}^{(i)}

- Simple and robust
- Can be slow to converge
- Good for unstable systems

**Broyden Mixing:**

Uses information from multiple previous densities to predict better mix.
Much faster convergence in most cases.

**Guidance:**

- For first calculation: try 'linear'
- For production: 'broyden' usually much better
- If 'broyden' oscillates: switch to 'linear'

**Related code:** ``source/mix.f90::type mix``

**See also:** :ref:`theory/scf_cycle`

beta (in ``&mix``)
------------------

**Type:** Real

**Purpose:** Density mixing parameter (written α in the formulas below)

**Typical range:** 0.01-0.8

**Default:** 0.1

**Example:**

.. code-block:: fortran

   beta = 0.3  ! Conservative

**Meaning:**

- New density weight in mixing formula
- Larger α -> larger step (faster but riskier)
- Smaller α -> smaller step (slower but more stable)

**Guidance:**

- Start with α = 0.5
- If oscillating: reduce to 0.2-0.3
- If converging too slowly: increase to 0.6-0.7

**Notes:**

- For 'linear' mixing: critical parameter for convergence speed
- For 'broyden' mixing: less critical (auto-adjusted somewhat)

**Related code:** ``source/mix.f90::mixpq()``

SCF Convergence Behavior
========================

**Typical sequence:**

.. code-block:: text

   Iteration 1: ΔQ = 0.5 Ry
   Iteration 2: ΔQ = 0.1 Ry (5x better)
   Iteration 3: ΔQ = 0.02 Ry (5x better)
   ...
   Iteration 10: ΔQ = 0.5e-5 Ry (converged!)

**Exponential convergence (good):**

- ΔQ decreases by constant factor each iteration
- Typical factor: 3-10 (depends on mixing)

**Linear convergence (OK):**

- ΔQ decreases linearly
- Slower but still acceptable

**Oscillation (bad):**

- ΔQ alternates between small and large values
- Usually indicates mixing parameter too aggressive
- Solution: reduce α, switch to linear mixing

**Stagnation (very bad):**

- ΔQ plateaus at high value, doesn't decrease
- Usually indicates problem with Hamiltonian or structure
- Solution: check element/potential parameters, cluster geometry

Troubleshooting SCF Convergence
===============================

**Problem: Doesn't converge at all**

Solutions:
- Reduce α from 0.5 to 0.2
- Switch from 'broyden' to 'linear' mixing
- Check element parameters and lattice constant

**Problem: Converges slowly (100+ iterations)**

Solutions:
- For 'linear': increase α to 0.6-0.7
- Try 'broyden' mixing instead
- Check if problem inherently difficult (bad initial guess)

**Problem: Oscillates between states**

Solutions:
- Reduce α significantly (0.2-0.3)
- Use 'linear' mixing (more stable)

**Problem: Energy or forces unphysical**

This usually indicates non-convergent SCF.
Solutions:
- Ensure conv_thr is tight enough
- Check Hamiltonian parameters
- Verify cluster size adequate (increase r2 or nx,ny,nz)

Additional &self Parameters
===========================

The following parameters from the ``&self`` namelist are less commonly used but
important for specific situations.

ws
--

**Type:** Real array (one value per inequivalent atom)

**Purpose:** Wigner-Seitz sphere radii (in Bohr) for each inequivalent atom,
overriding the lattice-derived value

**Example:**

.. code-block:: fortran

   ws(1) = 2.67
   ws(2) = 2.40

**Notes:**

- If not set, WS radii are derived from the lattice constant and structure.
- Useful when atoms have different sizes (e.g. surface or impurity calculations).

**See also:** ``ws_all`` below; ``sws`` in ``&lattice``

ws_all
------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Apply a single uniform Wigner-Seitz radius (from ``&lattice``) to
all atoms, overriding per-atom values

conv_thr
--------

**Type:** Real

**Default:** ``0.5e-8``

**Purpose:** SCF convergence threshold: the iteration stops with "Converged!" when the printed ``Diff`` is below it.

soc_scale
---------

**Type:** Real

**Default:** ``1.0``

**Purpose:** Global scaling factor applied to all spin-orbit coupling (SOC)
matrix elements.  Values less than 1 reduce the SOC strength; values greater
than 1 enhance it.

**Example:**

.. code-block:: fortran

   soc_scale = 0.5   ! Half-strength SOC

**Notes:**

- Only active when ``nsp = 2`` or ``4`` (relativistic calculations).
- Useful for separating orbital and spin effects, or for scanning SOC strength.

fix_soc
-------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Fix the spin-orbit coupling matrix elements at their initial values
throughout the SCF cycle (no self-consistent update of the SOC part)

**Notes:**

- Useful for perturbative SOC studies where self-consistency of the SOC
  contribution is not required.

orbital_polarization
--------------------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Apply an orbital polarization correction (Brooks-Gunnarsson mean-field)
to enhance the orbital magnetic moment

**Notes:**

- Only relevant for magnetic calculations (``nsp > 1``).
- Can improve agreement with experiment for strongly correlated 3d/4f systems.

magnetic_mixing
---------------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Apply a separate mixing step for the magnetic (exchange-correlation)
part of the potential, decoupled from the charge mixing

**Notes:**

- Useful when the spin channels converge at different rates.
- Can stabilise calculations where the spin density oscillates but the charge
  density is already converged.

mixmag
------

**Type:** Real array (one value per inequivalent atom)

**Purpose:** Per-atom mixing parameter for the magnetic potential.  Overrides
the global ``beta`` for the magnetic channel on a per-site basis.

**Example:**

.. code-block:: fortran

   mixmag(1) = 0.2   ! Conservative magnetic mixing for atom 1

mixmag_all
----------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Apply the same magnetic mixing parameter (``mixmag(1)``) to all
inequivalent atoms

mix_all
-------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Mix all components of the potential (charge + magnetic + orbital)
simultaneously in a single mixing step

freeze
------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Freeze the potential at its current value — no SCF updates are
applied.  Useful for single-shot calculations starting from a converged potential.

all_inequivalent
----------------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Treat every atom in the cluster as symmetry-inequivalent, even if
the lattice has symmetry.  Forces independent potentials on all sites.

**Notes:**

- Substantially increases memory and CPU cost.
- Required for impurity calculations or when breaking lattice symmetry explicitly.

rigid_band
----------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Apply a rigid-band shift to the potential: shift the Fermi level
without updating the self-consistent potential.

**Notes:**

- Useful for electron/hole doping studies at fixed potential shape.

nstep
-----

**Type:** Integer

**Default:** ``1``

**Purpose:** Maximum number of SCF steps; the loop ends earlier when ``Diff`` is below ``conv_thr``.

init
----

**Type:** Integer

**Default:** ``0``

**Purpose:** Initialisation mode for the SCF potential

**Allowed values:**

- ``0``: Start from scratch (free-atom superposition)
- ``1``: Read potential from a previous run (restart)

cold
----

**Type:** Logical

**Default:** ``.false.``

**Purpose:** **Cold start** — perform an ASA (Atomic Sphere Approximation)
extraction to derive initial tight-binding potential parameters from the
all-electron atomic solution before entering the SCF loop.

**Notes:**

- Useful when starting a new element not in the parameter database.
- Increases setup time but often gives a better starting point.

Advanced SCF Options
====================

**Accelerated convergence techniques:**

Some codes support:

- Kerker mixing (k-dependent weighting)
- Pulay mixing (residual minimization)
- Anderson mixing

RS-LMTO currently supports linear and Broyden;
others may be added in future versions.

**Temperature effects:**

If non-zero temperature (``&energy``):

- SCF may converge more slowly
- Smeared Fermi surface -> fewer oscillations
- Can sometimes help or hinder

**Constraint fields:**

For magnetic impurity or fixed-magnetization calculations:

- May need different mixing/convergence
- Less stability if constraints are too tight
- See :ref:`keywords/constraints` for the ``&constraints`` namelist

Provenance
==========

SCF parameters defined in:

- **Type definition:** ``source/self.f90::type self``
- **Mixing:** ``source/mix.f90::type mix``, ``mix.f90::mixpq()``
- **SCF loop and convergence check:** ``source/self.f90::run()``

Magnetic seed and per-atom options (&self)
==========================================

.. list-table::
   :header-rows: 1
   :widths: 34 10 18 38

   * - Key
     - Type
     - Default
     - Meaning
   * - ``magnetic_seed_enable``
     - logical
     - ``.false.``
     - Add a temporary field to the exchange splitting of every atomic solve during the
       first ``magnetic_seed_steps`` SCF iterations, to break the symmetry of a nonmagnetic
       start. Cannot be combined with ``&constraints``.
   * - ``magnetic_seed_steps``
     - integer
     - 0
     - Number of SCF iterations with the seed; must be at least 1 when the seed is on and
       leave at least one iteration without it.
   * - ``magnetic_seed_field``
     - real (Ry)
     - 0.0
     - Seed field, in the same convention as the constraining field ``B_fsm``; must be
       nonzero when the seed is on.
   * - ``force_nonmagnetic(i)``
     - logical (nrec)
     - ``.false.``
     - Keep recursion atom i nonmagnetic throughout the SCF.
   * - ``r05_p_tolerance``
     - real
     - -1.0
     - Fixed-U convergence target for the occupation matrix; values ``<= 0`` use the
       built-in target (1e-7 k-space, 5e-7 real space), values above 1 stop the run.
   * - ``r05_v_tolerance``
     - real (Ry)
     - -1.0
     - Fixed-U convergence target for the U potential; values ``<= 0`` use 1e-7 Ry,
       values above 1 stop the run.
   * - ``rb(i)``
     - integer (nrec)
     - 2
     - Read but not used.

``var`` in ``&mix`` is read but not used.

See Also
========

- :ref:`keywords/control_parameters` - Related: nsp, llsp, lld
- :ref:`theory/scf_cycle` - Detailed SCF theory
- :doc:`../user_guide/examples` - Example SCF parameter sets
- :doc:`../user_guide/input_files` - Input file format
