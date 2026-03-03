.. _keywords/element_parameters:

======================================
Element Parameters (&element namelist)
======================================

Overview
========

The ``&element`` namelist defines the chemical identity and electronic structure
of each **element type** used in the calculation.  One ``&element`` block is
required for every distinct chemical species present in the system.

The element database provides the fundamental atomic data from which the
LMTO tight-binding parameters (``&par``) are generated: atomic number, electron
counts, and principal quantum numbers for each angular-momentum channel.

.. code-block:: fortran

   &element
     symbol        = 'Fe'
     atomic_number = 26
     core          = 18
     valence       = 8
     f_core        = 0
     num_quant_s   = 4
     num_quant_p   = 4
     num_quant_d   = 3
   /

Parameters
==========

symbol
------

**Type:** Character string (len=10)

**Purpose:** Chemical symbol identifying the element (e.g. ``'Fe'``, ``'Co'``,
``'Pt'``).  Must match the symbol used in the ``&lattice`` / ``&par`` namelists
to associate structural positions with electronic parameters.

**Example:**

.. code-block:: fortran

   symbol = 'Fe'

**Notes:**

- Case-sensitive; use the standard two-character symbol (first letter uppercase,
  second lowercase).
- For empty spheres (interstitial atoms in the ASA), use a dummy symbol such as
  ``'Es'``.

atomic_number
-------------

**Type:** Real (stored as real for compatibility with alloy concentrations)

**Purpose:** Atomic number :math:`Z` of the element.

**Example:**

.. code-block:: fortran

   atomic_number = 26.0   ! Iron

**Notes:**

- Used to determine the nuclear potential :math:`V_{\text{nuc}}(r) = -2Z/r`
  inside the atomic sphere.
- For virtual crystal approximation (VCA) alloys, fractional values can be used
  (e.g. ``atomic_number = 26.5`` for a 50/50 Fe-Co virtual atom).

core
----

**Type:** Real

**Purpose:** Number of **core electrons** (those treated as frozen, not
participating in the valence self-consistency).

**Example:**

.. code-block:: fortran

   core = 18   ! [Ar] core for 3d transition metals

**Notes:**

- Core electrons are excluded from the charge density integration.
- The split between core and valence should be chosen so that the highest core
  level is well separated in energy from the lowest valence level.
- Typical choices for 3d transition metals: 18 (Ar core), leaving 3d + 4s in
  the valence.

valence
-------

**Type:** Real

**Purpose:** Number of **valence electrons** that are treated self-consistently.

**Example:**

.. code-block:: fortran

   valence = 8   ! 3d^6 4s^2 for Fe

**Notes:**

- ``core + valence`` must equal the total electron count ``atomic_number``.
- For spin-polarised calculations, ``valence`` is the total (spin-up + spin-down)
  valence count; the spin splitting is determined self-consistently.

f_core
------

**Type:** Integer

**Default:** ``0``

**Purpose:** Number of **f-core electrons** to include in the frozen core
(relevant for lanthanides and actinides where the 4f/5f shell may be kept in the
core).

**Example:**

.. code-block:: fortran

   f_core = 14   ! Full 4f shell in core for Lu

**Notes:**

- Set to ``0`` for all elements up to and including the 5d series.
- For rare-earth elements, setting ``f_core = 14`` (fully filled 4f in core)
  gives a good approximation for most structural and magnetic properties.

num_quant_s / num_quant_p / num_quant_d
----------------------------------------

**Type:** Integer

**Purpose:** Principal quantum numbers :math:`n` for the s, p, and d orbitals
included in the LMTO basis.

**Example for Fe (3d transition metal):**

.. code-block:: fortran

   num_quant_s = 4   ! 4s orbital
   num_quant_p = 4   ! 4p orbital
   num_quant_d = 3   ! 3d orbital

**Typical values by element group:**

.. list-table::
   :widths: 30 20 20 20
   :header-rows: 1

   * - Series
     - num_quant_s
     - num_quant_p
     - num_quant_d
   * - 3d transition metals (Sc–Cu)
     - 4
     - 4
     - 3
   * - 4d transition metals (Y–Ag)
     - 5
     - 5
     - 4
   * - 5d transition metals (La, Hf–Au)
     - 6
     - 6
     - 5
   * - sp metals (Al, Si, …)
     - n (valence)
     - n (valence)
     - n-1

**Notes:**

- These quantum numbers define which radial functions are used in the
  LMTO expansion and therefore which energy levels are included in the basis.
- Incorrect quantum numbers give wrong band centres and potential parameters.
- When in doubt, consult the ``source/element_data/`` directory for
  pre-defined database values.

Multiple Element Blocks
=======================

For a binary compound each element type needs its own ``&element`` block:

.. code-block:: fortran

   &element
     symbol        = 'Fe'
     atomic_number = 26
     core          = 18
     valence       = 8
     f_core        = 0
     num_quant_s   = 4
     num_quant_p   = 4
     num_quant_d   = 3
   /

   &element
     symbol        = 'Co'
     atomic_number = 27
     core          = 18
     valence       = 9
     f_core        = 0
     num_quant_s   = 4
     num_quant_p   = 4
     num_quant_d   = 3
   /

The ``symbol`` field links each ``&element`` block to the corresponding entry
in the ``&par`` (potential parameters) and ``&lattice`` namelists.

Provenance
==========

- **Namelist definition:** ``source/include_codes/namelists/element.f90``
- **Type definition:** ``source/element.f90``
- **Reading and setup:** ``source/symbolic_atom.f90``
- **Integration with potential:** ``source/potential.f90``

See Also
========

- :ref:`keywords/basis_parameters` — ``&par`` namelist (tight-binding potential
  parameters derived from the element database)
- :ref:`keywords/control_parameters` — ``nsp`` determines how many spin channels
  are used per element
- :doc:`../user_guide/input_files` — Input file structure and element ordering
