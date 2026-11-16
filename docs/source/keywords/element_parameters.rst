.. _keywords/element_parameters:

======================================
Element Parameters (&element namelist)
======================================

Overview
========

The ``&element`` namelist defines the chemical identity and electronic structure
of each **element type** used in the calculation.  It is not part of the main
input file: it is read from the file ``<label>.nml`` of the element (together
with its ``&par`` namelist), where ``<label>`` is an entry of ``label(:)`` in
``&atoms``. The file is looked up in the directory ``&atoms database`` (default
``./``) and then in the global database folder; a ``<label>_out.nml`` file from a
previous run is read first when it exists. Each file holds one ``&element``
group.

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

**Purpose:** Chemical symbol of the element (e.g. ``'Fe'``, ``'Co'``,
``'Pt'``). It names the per-atom output files, e.g. ``<symbol>_dos.out``.
The element is selected by its ``&atoms label``, not by ``symbol``.

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

Read and stored, but currently not used. The core charge comes from the radial
solver (``qv = Z - qc`` in ``source/self.f90``).

**Example:**

.. code-block:: fortran

   core = 18

valence
-------

**Type:** Real

**Purpose:** Number of **valence electrons** that are treated self-consistently.

**Example:**

.. code-block:: fortran

   valence = 8   ! 3d^6 4s^2 for Fe

**Notes:**

- ``valence`` is the electron count used in the charge-neutrality and Fermi-level routines; the code does not check it against ``atomic_number`` and ``core``.
- For spin-polarised calculations, ``valence`` is the total (spin-up + spin-down)
  valence count; the spin splitting is determined self-consistently.

f_core
------

**Type:** Integer

**Default:** ``-1`` (not ``0``). The f-core treatment of ``source/self.f90`` is applied whenever ``f_core /= 0`` and ``lmax < 3``, so an omitted ``f_core`` is not the same as ``0``: always set it.

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

num_quant_s / num_quant_p / num_quant_d / num_quant_f
-----------------------------------------------------

**Type:** Integer

**Purpose:** Principal quantum numbers :math:`n` for the s, p, d and f orbitals
included in the LMTO basis; ``num_quant_f`` is used with an f basis (``lmax = 3``).

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
- When in doubt, start from the ``example/`` deck of a similar element.

Multiple elements
=================

For a binary compound each element type has its own file, ``Fe.nml`` and
``Co.nml``, each with one ``&element`` group (and its ``&par`` group). The main
input names them in ``&atoms``:

.. code-block:: fortran

   &atoms
     database = './'
     label(1) = 'Fe'
     label(2) = 'Co'
   /

Two ``&element`` groups in one file are not read: only the first group of a file
is used.

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
