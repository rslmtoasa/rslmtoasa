.. _keywords/energy_mesh:

=======================================
Energy Mesh Parameters (&energy)
=======================================

Overview
========

The ``&energy`` namelist controls the energy grid used for integrating density of states
and calculating electronic properties. Proper energy mesh setup is critical for accuracy.

Parameters
==========

channels_ldos
-------------

**Type:** Integer

**Purpose:** Number of energy mesh points for DOS integration

**Typical range:** 200-1000

**Default:** 500

**Example:**

.. code-block:: fortran

   channels_ldos = 500

**Meaning:**

- Controls energy resolution of DOS
- More points = finer resolution = slower but more accurate
- Energy mesh typically spans ``energy_min`` to ``energy_max``

**Guidance:**

- Quick tests: 200-300 points
- Production: 500-1000 points
- Very fine structure: 1000+ points

**Related code:** ``source/energy.f90::e_mesh()``

energy_min
----------

**Type:** Real (Rydberg)

**Purpose:** Lower energy limit for DOS/properties calculation

**Typical range:** -3 to -1 Ry

**Default:** -3.0

**Example:**

.. code-block:: fortran

   energy_min = -2.5

**Meaning:**

- Lowest energy point in mesh
- Should be below lowest band
- Extends mesh well below Fermi level

**Notes:**

- For ground state properties: use E_min ≈ (lowest band) - 1 Ry
- For excited state or DOS visualization: extend lower
- Impact on DOS far below Fermi level is usually small

energy_max
----------

**Type:** Real (Rydberg)

**Purpose:** Upper energy limit for DOS/properties calculation

**Typical range:** 0 to +3 Ry

**Default:** +3.0

**Example:**

.. code-block:: fortran

   energy_max = 2.5

**Meaning:**

- Highest energy point in mesh
- Should be above highest band of interest
- Extends mesh above Fermi level

**Notes:**

- For band structure: need to reach top of bands
- For occupied states only: E_max ≈ E_Fermi sufficient
- For unoccupied states: extend higher (2-3 Ry above E_F)

fermi
-----

**Type:** Real (Rydberg)

**Purpose:** Initial guess or fixed value for Fermi energy

**Typical range:** -1 to +1 Ry

**Default:** 0.0

**Example:**

.. code-block:: fortran

   fermi = -0.3

**Meaning:**

- If ``fix_fermi = .false.``: initial guess for SCF
- If ``fix_fermi = .true.``: fixed throughout calculation

**Notes:**

- For metals: typically near 0 Ry
- For semiconductors: gap is usually near +-1-2 Ry
- Good initial guess accelerates convergence

fix_fermi
---------

**Type:** Logical

**Purpose:** Fix Fermi level to specified value

**Default:** .false.

**Example:**

.. code-block:: fortran

   fix_fermi = .true.
   fermi = -0.25

**Meaning:**

- .true.: Fermi level stays constant (useful for impurity problems)
- .false.: Fermi level adjusted to maintain correct electron count

**Notes:**

- For bulk: usually .false.
- For isolated impurity: often .true.
- Allows exploration of non-physical electron counts

fix_fermi_to_region (B7, calctype='L')
---------------------------------------

**Type:** Character string

**Purpose:** Pin the Fermi level to a named frozen region instead of the
default free E_F from cluster charge neutrality, for an interface
(``calctype='L'``) run.

**Default:** ``''`` (empty — free E_F)

**Example:**

.. code-block:: fortran

   fix_fermi_to_region = 'A'

**Notes:**

- In ``&charge`` namelist (not ``&energy`` — grouped here with
  ``fix_fermi`` for discoverability since they're the same kind of
  decision).
- Companion to ``fermi_a``/``fermi_b`` below when checking two-sided
  alignment consistency.

**Related code:** ``source/charge.f90``, ``docs/dev/CONTRACT_FROZEN_REGION.md``

fermi_a, fermi_b (B7, calctype='L')
-------------------------------------

**Type:** Real (Ry)

**Purpose:** Optional per-region source Fermi level for the interface
alignment consistency check (``align_regions``, gate G-B7-2).

**Default:** ``huge(1.0)`` (sentinel — not supplied)

**Example:**

.. code-block:: fortran

   fermi_a = -0.0899
   fermi_b = -0.0899

**Notes:**

- In ``&charge`` namelist. Leave unset unless you are supplying
  independently-converged source E_F values for each region to
  cross-check the alignment solver.

**Related code:** ``source/charge.f90``, ``docs/dev/CONTRACT_FROZEN_REGION.md``

compensation_profile, bias (B7, calctype='L')
------------------------------------------------

**Type:** Character string / Real (Ry)

**Purpose:** ``compensation_profile`` selects the charge-compensation
weighting scheme used by ``interfacepot`` when balancing the two-sided
electrostatics (B7 §1.5, gate G-B7-3 — tolerance calibration still open,
see ``tests/KNOWN_ISSUES.md``). ``bias`` is a B7 §6 hook for an applied
bias (``target_step = contact_potential + bias``); namelist plumbing only,
not yet wired to a physical effect.

**Allowed values (compensation_profile):** ``'nef_weighted'`` (the only
implemented profile — side-resolved :math:`N(E_F)` weighting)

**Default:** ``compensation_profile = 'nef_weighted'``, ``bias = 0.0``

**Example:**

.. code-block:: fortran

   compensation_profile = 'nef_weighted'
   bias = 0.0

**Notes:**

- In ``&charge`` namelist.
- ``bias`` is exposed now so it is a namelist choice rather than a future
  new code path, but it is inactive until a dedicated task wires it to
  the PM kernel — do not expect it to do anything yet.

**Related code:** ``source/charge.f90`` (``interfacepot``)

edel
----

**Type:** Real (Rydberg)

**Purpose:** Energy mesh spacing (automatic calculation)

**Meaning:**

- Computed from energy_min, energy_max, channels_ldos
- edel = (energy_max - energy_min) / (channels_ldos - 1)

**Example calculation:**

.. code-block:: fortran

   energy_min = -2.0
   energy_max = +2.0
   channels_ldos = 401
   ! edel = 4.0 / 400 = 0.01 Ry automatically

ene (array)
-----------

**Type:** Real array (dimension: channels_ldos)

**Purpose:** Actual energy mesh values

**Meaning:**

- Automatically generated as linearly spaced array
- ene(i) = energy_min + (i-1) * edel

**Notes:**

- Read-only; don't set manually
- Used internally for DOS/property calculations

ik1, nv1
--------

**Type:** Integer

**Purpose:** (Internal) Energy mesh indices for convergence

**Notes:**

- Automatically set; not typically user-adjustable

Typical Energy Mesh Setup Examples
==================================

**For Bulk Metal (Fe):**

.. code-block:: fortran

   &energy
      energy_min = -2.5,
      energy_max = +2.5,
      channels_ldos = 501,  ! Fine mesh
      fix_fermi = .false.   ! Adjust E_F
   /

**For Semiconductor (Si, Gap ~1 Ry):**

.. code-block:: fortran

   &energy
      energy_min = -3.0,
      energy_max = +1.0,    ! Above conduction band minimum
      channels_ldos = 401,
      fix_fermi = .false.
   /

**For Impurity Problem (Fixed E_F):**

.. code-block:: fortran

   &energy
      energy_min = -2.0,
      energy_max = +2.0,
      channels_ldos = 401,
      fix_fermi = .true.,
      fermi = -0.2          ! Fixed value
   /

Convergence with Respect to Energy Mesh
========================================

Recommended convergence test:

1. Run calculation with ``channels_ldos = 200``
2. Run with ``channels_ldos = 400``
3. Run with ``channels_ldos = 800``
4. Compare DOS, total energy, magnetic moments
5. Use first value where results plateau

**Typical behavior:**

- Energy changes: ~1-10 mRy when doubling mesh
- DOS shape stabilizes quickly (100-200 points usually adequate)
- Very fine structure (sharp peaks): may need 500+ points

Broadening Parameter
====================

**Related to energy mesh:**

The imaginary part ``η`` in Green's function:

.. math::

   G(E) = (E - H + i\eta)^{-1}

Controls spectral broadening. Smaller η -> sharper features -> slower convergence.

Typical: $\eta = 0.001$ to $0.01$ Ry.

(Set in recursion or green's function routines; not direct keyword.)

Provenance
==========

Energy mesh parameters defined in:

- **Type definition:** ``source/energy.f90::type energy``
- **Reading:** ``source/energy.f90::build_from_file()``
- **Mesh generation:** ``source/energy.f90::e_mesh()``
- **Usage:** ``source/green.f90``, ``source/density_of_states.f90``

See Also
========

- :ref:`keywords/control_parameters` - Related: llsp, lld (recursion depth)
- :ref:`theory/green_functions` - Green's function theory
- :doc:`../user_guide/examples` - Typical energy mesh values
