.. _keywords/scf_settings:

=======================================
SCF Convergence Parameters (&self)
=======================================

Overview
========

The ``&self`` namelist controls self-consistent field (SCF) loop settings,
including density mixing, convergence criteria, and iteration limits.

Parameters
==========

mixing
------

**Type:** Character string

**Purpose:** Type of density mixing algorithm

**Allowed values:** 'linear', 'broyden'

**Default:** 'linear'

**Example:**

.. code-block:: fortran

   mixing = 'broyden'

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

alpha
-----

**Type:** Real

**Purpose:** Density mixing parameter

**Typical range:** 0.1-0.8

**Default:** 0.5

**Example:**

.. code-block:: fortran

   alpha = 0.3  ! Conservative

**Meaning:**

- New density weight in mixing formula
- Larger α → larger step (faster but riskier)
- Smaller α → smaller step (slower but more stable)

**Guidance:**

- Start with α = 0.5
- If oscillating: reduce to 0.2-0.3
- If converging too slowly: increase to 0.6-0.7

**Notes:**

- For 'linear' mixing: critical parameter for convergence speed
- For 'broyden' mixing: less critical (auto-adjusted somewhat)

**Related code:** ``source/mix.f90::mix_charge()``

broyden_history
---------------

**Type:** Integer

**Purpose:** Number of previous densities kept for Broyden mixing

**Typical range:** 3-20

**Default:** 10

**Example:**

.. code-block:: fortran

   broyden_history = 8

**Meaning:**

- Higher number: uses more history → better prediction
- Lower number: faster but less information

**Guidance:**

- Start with 10 (reasonable balance)
- If convergence stalls: increase to 15-20
- If very slow: reduce to 5-8

**Notes:**

- Only used if ``mixing = 'broyden'``
- Increases memory usage (~10% for large systems)
- Significant impact on convergence rate

**Related code:** ``source/mix.f90``

max_iterations (also in &control)
----------------------------------

**Type:** Integer

**Purpose:** Maximum number of SCF iterations

**Typical range:** 50-200

**Default:** 100

**Example:**

.. code-block:: fortran

   max_iterations = 80

**Notes:**

- Can be in &control or &self (check your version)
- Calculation stops when either iteration limit or convergence reached

dq_tol (also in &control)
--------------------------

**Type:** Real (Ry)

**Purpose:** Charge density convergence tolerance

**Typical range:** 1e-6 to 1e-4

**Default:** 1e-5

**Example:**

.. code-block:: fortran

   dq_tol = 1.0e-5

**Meaning:**

.. math::

   \sum_i |q_i^{(n+1)} - q_i^{(n)}| < \text{dq\_tol}

**Guidance:**

- 1e-5 or better: Production quality
- 1e-4: Reasonable for most purposes
- 1e-6+: Overkill; diminishing returns

**Notes:**

- Smaller tolerance = more iterations needed
- Check in &control; may override &self value

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
- Increase Broyden history (if using Broyden)

**Problem: Energy or forces unphysical**

This usually indicates non-convergent SCF.
Solutions:
- Ensure dq_tol is tight enough
- Check Hamiltonian parameters
- Verify cluster size adequate (increase r2 or nx,ny,nz)

Additional &self Parameters
==========================

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

**Default:** ``1.0e-5``

**Purpose:** SCF convergence threshold used internally by the ``&self`` module
(analogous to ``dq_tol`` in ``&control``)

**Notes:**

- In most input files only ``dq_tol`` (in ``&control``) needs to be set.
- If both are present, the more stringent one governs convergence.

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
the global ``alpha`` for the magnetic channel on a per-site basis.

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

**Purpose:** Number of SCF steps to perform (alternative to ``max_iterations``
in ``&control``).  If both are set, the minimum is used.

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
- Smeared Fermi surface → fewer oscillations
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
- **Mixing:** ``source/mix.f90::type mix``, ``mix.f90::mix_charge()``
- **SCF loop:** ``source/self.f90::process()``
- **Convergence check:** ``source/self.f90::converge_scf()``

See Also
========

- :ref:`keywords/control_parameters` - Related: nsp, llsp, lld
- :ref:`theory/scf_cycle` - Detailed SCF theory
- :doc:`../user_guide/examples` - Example SCF parameter sets
- :doc:`../user_guide/input_files` - Input file format
