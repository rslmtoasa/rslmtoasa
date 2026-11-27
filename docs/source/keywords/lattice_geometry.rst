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

**Default:** none (not set in ``lattice%restore_to_default``); give it in every input.

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

Read but currently not used. The namelist value is not stored: the code sets
``nbulk`` itself, as the number of atoms of type bulk (0 for a bulk run; see
the member documentation in ``source/lattice.f90``). The number of atoms of
the bulk unit cell is ``nbulk_bulk`` (below).

**Related code:** ``source/lattice.f90::type lattice``,
``source/lattice_cluster.f90``

r2
--

**Type:** Real (Ångströms²)

**Purpose:** Cluster cutoff radius squared

**Typical range:** 10-30 Å²

**Default:** 0 (not given): with a built-in ``crystal_sym`` it becomes ``ct(1)**2``, where ``ct`` defaults to ``alat + 0.1``.

**Meaning:**

- Defines hopping range from central atom
- Atoms within sqrt(r2) distance are included
- Controls neighborhood for Hamiltonian construction

**Example:**

.. code-block:: fortran

   r2 = 13.0  ! First-neighbor shell only

**Notes:**

- r2 = 7.566 Å² is the squared first-neighbour distance of fcc Pd (a = 3.89 Å)
- r2 = 15.13 Å² is the squared second-neighbour distance of fcc Pd
- Use "5 th neighs." as comment suggests (code note) for bulk limit

**Related code:** ``source/lattice_cluster.f90``

**See also:** :doc:`../code_structure` for details on cluster construction

Bravais Lattice Types
---------------------

Though not explicit in namelists, typical structures:

**Cubic (BCC, FCC, SC):**

.. code-block:: fortran

   crystal_sym = 'bcc'
   alat = 2.86120   ! For bcc iron

**Hexagonal (HCP):**

.. code-block:: fortran

   crystal_sym = 'hcp'
   celldm = 1.633    ! c/a ratio

**Tetragonal or lower:**

Use ``crystal_sym = 'file'`` and give the lattice vectors ``a`` and the basis
``crd`` in ``lattice.nml`` (read by ``lattice%build_from_lattice``).

Wigner-Seitz Sphere Parameters
==============================

ws_r (from &par namelist)
==========================

**Type:** Real (Bohr)

**Purpose:** Wigner-Seitz sphere radius (atomic sphere size)

**Typical range:** 2-3 Bohr

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

**celldm** - c/a ratio for ``crystal_sym = 'hcp'`` (1.633, with a warning, when 0)

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
     - ``'legacy'`` or ``'strux_lib'``; any other value runs the legacy path. ``strux_lib`` is required for ``lmax = 3`` (f basis) and
       for ``magnetic_representation = 'gbt_single_q'``.
   * - ``screening``
     - character
     - ``'default'``
     - Screening of the ``strux_lib`` structure constants: ``'default'`` uses the built-in
       table and ignores the input ``screening_alpha``; ``'manual'`` uses
       ``screening_alpha``; ``'sigma'`` uses ``screening_sigma``; ``'fitted'`` fits the
       screening constants. Other values stop the run.
   * - ``screening_alpha(l)``
     - real (lmax+1)
     - built-in table
     - l-resolved screening constants for ``'manual'``.
   * - ``screening_sigma(l)``
     - real (lmax+1)
     - 0.7
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
     - Crystal symmetry: ``'bcc'``, ``'b2'``, ``'fcc'``, ``'hcp'``, ``'fcc2'``, ``'fcc3'`` or ``'file'``; any other value stops the run.
   * - ``celldm``
     - c/a ratio of the hcp cell (1.633, with a warning, when 0).
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
   * - ``pair_radius``
     - Real, in the units of the cluster coordinates (alat), default 0 (off). With a value above 0 the
       exchange post-processing (``exchange``, ``exchange_import``, ``exchange_p2rs``) generates its pair list instead of reading ``njij`` and
       ``ijpair`` (giving both is fatal). For every recursion atom ``irec(k)`` the list holds all cluster atoms
       j with :math:`|r_j - r_i| \le` ``pair_radius`` (+ 1e-6), the self pair first, then by distance and atom
       index. The sphere must lie inside the cluster, :math:`|r_i - r_1|` + ``pair_radius`` :math:`\le \sqrt{r2}`,
       otherwise the run stops. A plain scalar: it may stand anywhere in ``&lattice`` (only ``njij`` and
       ``ijpair`` must follow ``ntype``, ``ct`` and ``r2``).
   * - ``pair_unique``
     - Logical, default ``.false.``. With ``pair_radius``, keep one representative per symmetry star of the pairs of
       each recursion atom: (i, j') is equivalent to (i, j) when an operation of the crystal that leaves the site i
       invariant maps :math:`r_j - r_i` onto :math:`r_{j'} - r_i` and j, j' have the same type. The operations are the
       ones of the symmetry-reduced k mesh (spglib, ``setup_magnetic_symmetry``): for ``nsp = 2`` and ``4`` the magnetic
       group of the site moments as axial vectors, including the time-reversal-combined operations; for ``nsp = 1`` the crystal operations
       that map every site onto a site with the same type and moment. The moment of a site is its constraint target, else
       the unit ``potential%mom``. Needs spglib, ``calctype = 'B'`` and ``nsp`` 1, 2 or 4 (a stop otherwise).
       ``ijpair`` and ``njij`` then hold the representatives, so ``jij.out``, ``dij.out``, ``aij.out`` and the other
       exchange outputs list one pair per star; the full list is kept in ``lattice%pair_full``, the representative of
       each member in ``lattice%pair_rep`` and the rotation from the representative to the member in ``lattice%pair_rot``;
       the log prints the multiplicity of each star. The pair-loop route (``exchange``, ``exchange_import``) also writes
       ``exchange_stars.out``, one row per member of every star with its own vector :math:`r_j - r_i`, J, D and A of the
       representative carried over as follows. With R the Cartesian rotation (proper or improper) of the operation that maps the
       representative (i, j) onto (i, j'), the spins rotate with the proper part
       :math:`\tilde R = \det(R) R` (spin is an axial vector; a time reversal flips both spins of the bilinear pair energy and
       changes nothing), so the pair tensor :math:`J + [D]_\times + A` of the member is :math:`\tilde R (\cdot) \tilde R^T`:
       :math:`J' = J`, :math:`D' = \tilde R D`, :math:`A' = \tilde R A \tilde R^T`, for unitary and time-reversal-combined
       operations alike. With spin-orbit coupling D and A depend on the direction of the moments, which is why the
       operations are those of the moments and not of the lattice alone; the D vector and the anisotropic tensor of a
       representative do not carry over by the lattice group in general. Without it (``nsp = 1``) the spins are decoupled from the lattice and J, D and A of the representative are copied unrotated.
       With ``do_chi0`` the lattice sums, the type sums and the J0 sum run over the full stars (see ``do_chi0``).
   * - ``njijk``, ``ijktrio``
     - Number of atom trios and the trios (i, j, k) for spin-lattice interactions.
   * - ``nclu``, ``inclu``
     - Number of impurity sites and their coordinates (impurity cluster build).
   * - ``acr``, ``reduced_acr``
     - Cluster bookkeeping from the impurity build; no further meaning documented in the code.
   * - ``surftype``
     - Surface orientation as three integers read free-format, e.g. ``'1 1 1'`` or ``'0 0 1'``; four indices for ``crystal_sym = 'hcp'``.
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

- :ref:`keywords/basis_parameters` - Related: lmax, ws_r (``&par``)
- :ref:`keywords/control_parameters` - Related: lmax, calctype
- :doc:`../user_guide/input_files` - Input file format
- :doc:`../user_guide/examples` - Typical values for common materials
- :doc:`../code_structure` - Cluster construction details
