.. _keywords/lattice_geometry:

=======================================
Lattice & Geometry Parameters
=======================================

Overview
========

The ``&lattice`` namelist defines the crystal structure, lattice parameters, 
and cluster geometry for the calculation.

Parameters
==========

alat
----

**Type:** Real (Ångströms)

**Purpose:** Lattice constant

**Typical range:** 2-6 Å (system-dependent)

**Default:** 5.0

**Example:**

.. code-block:: fortran

   alat = 5.4307  ! Silicon lattice constant

**Notes:**

- In Ångströms (not atomic units)
- Critical for all length scales
- Should match experimental or first-principles values

**Related code:** ``source/lattice.f90::build_from_file()``

nbulk
-----

**Type:** Integer

**Purpose:** Number of atoms in bulk unit cell

**Typical range:** 1-4

**Default:** 1

**Meaning:**

- 1 for monatomic lattices (simple cubic, bcc, fcc)
- 2+ for compounds or complex structures
- Defines the base repeating unit before cluster expansion

**Example:**

.. code-block:: fortran

   nbulk = 2  ! Two atoms per unit cell

**Notes:**

- For diamond (Si, Ge): nbulk = 2
- For rock salt (NaCl): nbulk = 2
- For perovskite (ABO3): nbulk = 5

**Related code:** ``source/lattice.f90::type lattice``

nx, ny, nz
----------

**Type:** Integer (each)

**Purpose:** Cluster dimensions in lattice coordinates

**Typical range:** 3-11 (odd numbers for center)

**Default:** 5 (each)

**Meaning:**

- Defines number of unit cells in each Cartesian direction
- Center atom typically at (nx+1)/2, (ny+1)/2, (nz+1)/2
- Larger values -> larger cluster -> more expensive

**Example:**

.. code-block:: fortran

   nx = 7
   ny = 7
   nz = 1   ! Quasi-2D: surface-like

**Notes:**

- Usually odd values to have center atom
- Total atoms ≈ nbulk × nx × ny × nz
- Convergence study: vary nx, ny, nz to check bulk limit

**Related code:** ``source/lattice.f90``

r2
--

**Type:** Real (Ångströms²)

**Purpose:** Cluster cutoff radius squared

**Typical range:** 10-30 Å²

**Default:** 25.0

**Meaning:**

- Defines hopping range from central atom
- Atoms within sqrt(r2) distance are included
- Controls neighborhood for Hamiltonian construction

**Example:**

.. code-block:: fortran

   r2 = 13.0  ! First-neighbor shell only

**Notes:**

- r2 = 7.566 Å² ≈ distance to first neighbors in fcc (e.g., Au)
- r2 = 15.13 Å² ≈ distance to second neighbors
- Use "5 th neighs." as comment suggests (code note) for bulk limit

**Related code:** ``source/lattice_cluster.f90``

**See also:** :doc:`../code_structure` for details on cluster construction

Bravais Lattice Types
---------------------

Though not explicit in namelists, typical structures:

**Cubic (BCC, FCC, SC):**

.. code-block:: fortran

   alat = 2.87   ! For bcc iron
   nx = 7        ! Cubic cluster
   ny = 7
   nz = 7

**Hexagonal (HCP, graphene):**

.. code-block:: fortran

   nx = 7
   ny = 7
   nz = 3        ! Thinner in z for layers

**Tetragonal or lower:**

Use ``nx, ny, nz`` and ``alat`` to specify conventional cell dimensions.

Wigner-Seitz Sphere Parameters
==============================

ws_r (from &par namelist)
==========================

**Type:** Real (Ångströms)

**Purpose:** Wigner-Seitz sphere radius (atomic sphere size)

**Typical range:** 2-3 Å

**Default:** element/structure-dependent

**Meaning:**

- Defines size of non-overlapping sphere around each atom
- Larger sphere -> includes more electron density
- Related to ``sws`` in charge density calculations

**Example:**

.. code-block:: fortran

   ws_r = 2.827  ! Silicon

**Notes:**

- For cubic: WS radius ≈ alat × sqrt(π/6) for fcc, etc.
- ASA requires overlapping spheres to cover all space
- Too large -> overlaps excessive
- Too small -> misses density

**Constraint:**

Total sphere volume should be close to crystal volume:

.. math::

   N_{\text{atoms}} \times \frac{4}{3}\pi r_{\text{WS}}^3 \approx V_{\text{cell}}

Layered / interface geometry (calctype='L', B7)
================================================

These ``&lattice`` keys apply only when ``calctype='L'`` (interface mode,
``buildinterface`` — see :ref:`keywords/control_parameters`); they are
ignored otherwise.

nlay_a, nlay_b
--------------

**Type:** Integer

**Purpose:** Number of *atomic layers* in region A / region B of an
interface cluster, counted from the interface outward.

**Default:** ``0`` (both)

**Example:**

.. code-block:: fortran

   nlay_a = 4
   nlay_b = 4

**Notes:**

- These count physical atomic layers in ``&lattice`` — a different thing
  from the ``&charge`` namelist's own ``nlay_a``/``nlay_b``, which count
  *synthetic Madelung rows* for the electrostatics solve. **Raising the
  ``&charge`` pair without raising this one breaks alignment** (observed:
  V(B) off by ~0.45 Ry) — see ``tests/KNOWN_ISSUES.md`` and
  ``docs/DECISIONS.md``.
- Only the *first frozen non-vacuum* layer in each region is used as the
  gauge anchor for the alignment solver (``align_regions``,
  ``source/region_registry.f90``).

**Related code:** ``source/lattice_cluster.f90``, ``source/region_registry.f90``

region_b_kind
-------------

**Type:** Character string

**Purpose:** Physical kind of region B in an interface (``calctype='L'``)
cluster.

**Allowed values:** ``'metal'`` (region B is a second metallic reference,
loaded from ``&atoms label(:)``), ``'vacuum'`` (region B is semi-infinite
vacuum, with frozen potential parameters generated per run by
``vacuum_lead`` rather than hand-set empty spheres — B7.2/B7.6)

**Default:** ``'metal'``

**Example:**

.. code-block:: fortran

   region_b_kind = 'vacuum'

**Notes:**

- Case-insensitive; an unrecognized value is **fatal** at namelist read
  (deliberately — a misspelling like ``'vaccum'`` would otherwise silently
  fall through to the metallic path and report a plausible-looking wrong
  barrier height; see ``source/lattice_lifecycle.f90``, B7.6 comment).
- With ``'vacuum'``, ``vacuum_lead.f90`` generates the frozen empty-sphere
  potential parameters by running the code's own radial solver at constant
  V(r), checked against an analytic spherical-Bessel oracle
  (``tests/unit/test_vacuum_lead.f90``).
- See :doc:`../user_guide/examples/interface_fcccu111` for a worked
  ``A|vacuum`` example and the alignment-convergence behavior with buffer
  width.

**Related code:** ``source/vacuum_lead.f90``, ``source/lattice_lifecycle.f90``

Related Structural Parameters
=============================

Several structure-related parameters are in ``&lattice`` or ``&charge``:

**celldm** - Cell parameters (from PWSCF convention)

**gx, gy, gz** - Lattice distortion parameters (for structure optimization)

**a, b, c** - Direct/reciprocal lattice constants (for low-symmetry systems)

Structure constants (strux_lib)
===============================

.. list-table::
   :header-rows: 1
   :widths: 34 10 18 38

   * - Key
     - Type
     - Default
     - Meaning
   * - ``strux_backend``
     - character
     - ``'legacy'``
     - ``'legacy'`` or ``'strux_lib'``. ``strux_lib`` is required for ``lmax = 3`` (f basis) and
       for ``magnetic_representation = 'gbt_single_q'``.
   * - ``screening``
     - character
     - ``'default'``
     - Screening of the ``strux_lib`` structure constants: ``'default'`` and ``'manual'`` use
       ``screening_alpha``; ``'sigma'`` uses ``screening_sigma``; ``'fitted'`` fits the
       screening constants. Other values stop the run.
   * - ``screening_alpha(l)``
     - real (lmax+1)
     - built-in table
     - l-resolved screening constants for ``'default'``/``'manual'``.
   * - ``screening_sigma(l)``
     - real (lmax+1)
     - unset
     - l-resolved hard-core radii in units of the sphere radius, for ``'sigma'``.
   * - ``strux_want_sdot``
     - logical
     - ``.false.``
     - Also compute the energy derivative of the structure constants; needs
       ``screening = 'sigma'`` or ``'fitted'``.
   * - ``strux_solve_scale``
     - real
     - 9.0
     - Multiplier of the squared radius of the cluster on which the structure constants are
       solved.

Periodic cluster
================

.. list-table::
   :header-rows: 1
   :widths: 34 10 18 38

   * - Key
     - Type
     - Default
     - Meaning
   * - ``pbc``, ``b1``, ``b2``, ``b3``
     - logical
     - ``.false.``
     - Periodic boundary conditions; with all four true the cluster is a periodic supercell.
   * - ``n1``, ``n2``, ``n3``
     - integer
     - 0
     - Repeats of the primitive vectors ``a(:, 1..3)`` that span the periodic supercell.
   * - ``morton_sfc``
     - logical
     - ``.false.``
     - Reorder the bulk cluster along a Morton (Z-order) curve for memory locality in the
       recursion and Chebyshev kernels; results are unchanged up to round-off.

Cluster and bookkeeping keys
============================

These ``&lattice`` keys describe the cluster that pre-processing builds. Most are
written by the code to the output namelist and read back by later runs; they are
rarely set by hand. Where the code documents no further meaning, the table says so.

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Key
     - Meaning
   * - ``crystal_sym``
     - Crystal symmetry: ``'bcc'``, ``'fcc'``, ``'hcp'`` or ``'nsy'``.
   * - ``celldm``
     - Cell parameter written with the lattice; no further meaning documented in the code.
   * - ``primcell``
     - Coordinates of the atoms of the primitive cell.
   * - ``crd``
     - Basis (primitive-cell) coordinates, shape (3, ntot); the starting point of the cluster
       expansion.
   * - ``cr``
     - Cartesian coordinates of the expanded cluster, shape (3, kk).
   * - ``kk``
     - Number of cluster sites.
   * - ``ndim``, ``npe``
     - Size limits of the cluster before the cut (defaults 9900000 and 49); rarely changed.
   * - ``rc``
     - Cut radius of the cluster.
   * - ``ct``
     - Nearest-neighbour cut radii (Å).
   * - ``nmax``
     - Number of cluster sites with their own local Hamiltonian (impurities); 0 for
       ``calctype`` B and L.
   * - ``ntot``
     - Number of atoms in the primitive cell (1 for fcc and bcc without relaxation, 2 for hcp).
   * - ``ntype``
     - Number of inequivalent atoms, ``nbulk + nrec``.
   * - ``nbas``
     - Number of shells or atoms over which charge is distributed; ``reduced_nbas`` is the
       reduced count.
   * - ``nbulk_bulk``
     - Auxiliary count used to set ``ntype``, ``nbulk`` and ``ntot``.
   * - ``iz``
     - Atom type of each cluster site.
   * - ``izp``
     - Type label of each primitive-cell atom (species labels used by ``strux_lib``).
   * - ``num``
     - Bravais type of each atom.
   * - ``no``
     - Prototype class of each cluster site, used for the Wigner-Seitz radius lookup.
   * - ``ib``, ``iu``
     - Cluster site of each atom of the primitive cell (1 for a bulk material; chosen from
       the cluster for surfaces).
   * - ``irec``
     - Cluster sites treated self-consistently (one per inequivalent recursion atom).
   * - ``wav``
     - Wigner-Seitz radius (Å).
   * - ``njij``, ``ijpair``
     - Number of atom pairs and the pairs (i, j) for exchange interactions.
   * - ``njijk``, ``ijktrio``
     - Number of atom trios and the trios (i, j, k) for spin-lattice interactions.
   * - ``nclu``, ``inclu``
     - Number of impurity sites and their coordinates (impurity cluster build).
   * - ``acr``, ``reduced_acr``
     - Cluster bookkeeping from the impurity build; no further meaning documented in the code.
   * - ``surftype``
     - Surface orientation: ``'111'``, ``'110'`` or ``'001'``.
   * - ``nlay``, ``z``
     - Number of layers and the layer coordinates of a surface cluster.
   * - ``zmin``, ``zmax``, ``zstep``
     - Range and step of the layer coordinate in the surface cluster build.
   * - ``dx``, ``dy``, ``dz``, ``dw``
     - Surface indices of the surface cluster build.
   * - ``crsurf``, ``izsurf``, ``izpsurf``, ``nosurf``, ``izpo``
     - Surface-cluster bookkeeping; no further meaning documented in the code.

Provenance
==========

Lattice parameters defined and used in:

- **Definition:** ``source/lattice.f90::type lattice``
- **Reading:** ``source/lattice.f90::build_from_file()``
- **Cluster setup:** ``source/lattice_cluster.f90``
- **Geometry utilities:** ``source/math.f90`` (distance, angles)

See Also
========

- :ref:`keywords/basis_parameters` - Related: lmax, ws_r
- :ref:`keywords/control_parameters` - Related: nlim, verbose
- :doc:`../user_guide/input_files` - Input file format
- :doc:`../user_guide/examples` - Typical values for common materials
- :doc:`../code_structure` - Cluster construction details
