.. _keywords/frozen_magnon:

===================================================
Frozen-Magnon Parameters (&frozen_magnon Namelist)
===================================================

Overview
========

The ``&frozen_magnon`` namelist configures ``post_processing =
'frozen_magnon'`` (B1): a sweep of ``&hamiltonian``'s spin-spiral wavevector
``q_ss`` over a list of q-points, evaluating the total/band energy at each
and reporting a magnon dispersion :math:`\omega(q)`. It relies on the B1
generalized-Bloch-theorem (GBT) fix in the spiral Hamiltonian block — see
:ref:`keywords/hamiltonian_parameters` (``q_ss``, ``theta_ss``,
``magnetic_representation='gbt_single_q'``) and
``GBT_RS_LMTO_completion_blueprint.md``.

The namelist is **optional**. It only has an effect when
``post_processing = 'frozen_magnon'`` (see :ref:`keywords/output_options`).

.. code-block:: fortran

   &calculation
     post_processing = 'frozen_magnon'
   /
   &frozen_magnon
     mode           = 'mft'
     branch_mode    = 'acoustic'
     q_coordinates  = 'cartesian'
     q_file         = 'q_path.dat'
     output_file    = 'frozen_magnon.dat'
   /

**Known limitation:** the single-sublattice (``branch_mode='acoustic'``)
route is validated. The multi-sublattice ``branch_mode='auto'`` route's
acoustic branch is **not yet gapless at** :math:`\Gamma` (a Goldstone-mode
violation of roughly 0.28 Ry on two-sublattice bcc FeCo, unexplained,
deferred to B11) — see ``tests/KNOWN_ISSUES.md``.

----

Parameters
==========

mode
----

**Type:** Character string

**Allowed values:** ``'mft'``, ``'scf'``

**Default:** ``'mft'``

**Purpose:** Selects how each q-point after the first is evaluated.

**Example:**

.. code-block:: fortran

   mode = 'mft'

**Notes:**

- ``'mft'`` (magnetic force theorem): converge SCF fully at the first
  q-point (``q_ss_list(:,1)`` or the first row of ``q_file``), then reuse
  that converged potential for a single-iteration band-energy pass at
  every other q. Cheap — one real SCF convergence for the whole sweep.
- ``'scf'``: converge SCF fully and independently at every q-point.
  Expensive but does not assume the force-theorem approximation.
- The first q-point is **always** the reference point regardless of mode;
  every other point's E(q)/:math:`\omega`\ (q) is measured relative to it.

**Related code:** ``source/calculation.f90::post_processing_frozen_magnon()``

branch_mode
-----------

**Type:** Character string

**Allowed values:** ``'acoustic'``, ``'auto'``

**Default:** ``'acoustic'``

**Purpose:** Selects the magnon-branch construction.

**Example:**

.. code-block:: fortran

   branch_mode = 'acoustic'

**Notes:**

- ``'acoustic'``: single-sublattice adiabatic dispersion,
  :math:`\omega(q) = \Delta E(q) / \Delta m_z` with :math:`\Delta m_z
  \propto M_{\text{tot}} \sin^2\theta` (small cone) — the validated route.
- ``'auto'``: multi-sublattice magnon branches via the direct GBT
  frozen-magnon method (second derivatives of the force-theorem
  band-energy surface with respect to sublattice cone angles; Essenberger
  *et al.*, PRB **84**, 174425, Eq. 26). Writes
  ``frozen_magnon_branches.dat``/``frozen_magnon_modes.dat`` in addition
  to ``output_file``. **Known not gapless at** :math:`\Gamma` — see the
  overview note above.
- ``'auto'`` uses ``theta_probe`` and ``active_moment_threshold`` below;
  ``'acoustic'`` does not.

**Related code:** ``source/calculation.f90::post_processing_frozen_magnon_acoustic()``,
``::post_processing_frozen_magnon_auto()``

q_file
------

**Type:** Character string (file path)

**Default:** ``''`` (empty — use ``q_ss_list``/``n_q_points`` instead)

**Purpose:** Path to a file listing the q-points to sweep. Preferred over
inline ``q_ss_list`` for long paths.

**Example:**

.. code-block:: fortran

   q_file = 'q_path.dat'

**File format:**

.. code-block:: text

   n_q_points
   qx qy qz
   qx qy qz
   ...

**Notes:**

- Row units follow ``q_coordinates`` below.
- If ``q_file`` is blank, q-points are read from the namelist's own
  ``q_ss_list(3, n_q_points)`` array instead, which has a fixed
  compile-time column bound of 500 (a Fortran namelist limitation — see
  the comment in ``source/include_codes/namelists/frozen_magnon.f90``).

**Related code:** ``source/frozen_magnon.f90``

q_coordinates
-------------

**Type:** Character string

**Allowed values:** ``'cartesian'``, ``'direct'``

**Default:** ``'cartesian'``

**Purpose:** Units of the q-point rows in ``q_file`` (or ``q_ss_list``).

**Example:**

.. code-block:: fortran

   q_coordinates = 'direct'

**Notes:**

- ``'cartesian'``: rows are in Cartesian units of :math:`2\pi/\text{alat}`,
  matching ``&hamiltonian``'s ``q_ss`` convention directly.
- ``'direct'``: rows are reciprocal-lattice coordinates, converted to
  Cartesian after the lattice is built.

**Related code:** ``source/frozen_magnon.f90``

n_q_points
----------

**Type:** Integer

**Default:** ``0``

**Purpose:** Number of q-points to use from the inline ``q_ss_list`` array
(ignored when ``q_file`` is set).

output_file
-----------

**Type:** Character string

**Default:** ``'frozen_magnon.dat'``

**Purpose:** Output file receiving, per q-point: total energy, band
energy, per-sublattice moment magnitude, and :math:`\omega(q)`.

theta_probe
-----------

**Type:** Real (degrees)

**Default:** ``-1.0`` (sentinel; the code substitutes 20° when
:math:`\le 0`)

**Purpose:** Small-cone tilt angle used by the ``'auto'`` branch's
finite-difference probes of the force-theorem band-energy surface.

**Notes:**

- Only used when ``branch_mode = 'auto'``.
- Harmonic-regime assumption: keep small (the 20° default is the
  documented working point; see
  ``source/calculation.f90::post_processing_frozen_magnon_auto()``).

active_moment_threshold
-----------------------

**Type:** Real (:math:`\mu_B`)

**Default:** ``1.0e-6``

**Purpose:** Minimum per-sublattice moment magnitude for a sublattice to
be treated as "magnetically active" in the ``'auto'`` branch construction.
Sublattices below threshold are excluded from the probe set.

**Notes:**

- Only used when ``branch_mode = 'auto'``.
- Fatal error if no sublattice clears the threshold.

----

See Also
========

- :ref:`keywords/hamiltonian_parameters` — ``q_ss``, ``theta_ss``,
  ``magnetic_representation='gbt_single_q'`` (the GBT representation)
- :ref:`keywords/output_options` — ``post_processing = 'frozen_magnon'``
- ``docs/dev/plans/B1_gbt_frozen_magnons.md`` — the GBT bug-fix spec and
  frozen-magnon workflow design
- ``docs/DEVELOPER_MAP.md`` — one-line summary with file/line pointers
- ``tests/KNOWN_ISSUES.md`` — the ``branch_mode='auto'`` Goldstone gap
