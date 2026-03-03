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
- Larger values → larger cluster → more expensive

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

**Related code:** ``source/lattice.f90::setup_cluster()``

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
- Larger sphere → includes more electron density
- Related to ``sws`` in charge density calculations

**Example:**

.. code-block:: fortran

   ws_r = 2.827  ! Silicon

**Notes:**

- For cubic: WS radius ≈ alat × sqrt(π/6) for fcc, etc.
- ASA requires overlapping spheres to cover all space
- Too large → overlaps excessive
- Too small → misses density

**Constraint:**

Total sphere volume should be close to crystal volume:

.. math::

   N_{\text{atoms}} \times \frac{4}{3}\pi r_{\text{WS}}^3 \approx V_{\text{cell}}

Related Structural Parameters
=============================

Several structure-related parameters are in ``&lattice`` or ``&charge``:

**celldm** - Cell parameters (from PWSCF convention)

**gx, gy, gz** - Lattice distortion parameters (for structure optimization)

**a, b, c** - Direct/reciprocal lattice constants (for low-symmetry systems)

Provenance
==========

Lattice parameters defined and used in:

- **Definition:** ``source/lattice.f90::type lattice``
- **Reading:** ``source/lattice.f90::build_from_file()``
- **Cluster setup:** ``source/lattice.f90::setup_cluster()``
- **Geometry utilities:** ``source/math.f90`` (distance, angles)

See Also
========

- :ref:`keywords/basis_parameters` - Related: lmax, ws_r
- :ref:`keywords/control_parameters` - Related: nlim, verbose
- :doc:`../user_guide/input_files` - Input file format
- :doc:`../user_guide/examples` - Typical values for common materials
- :doc:`../code_structure` - Cluster construction details
