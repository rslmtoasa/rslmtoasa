.. _keywords/basis_parameters:

=======================================
Basis & Potential Parameters (&par)
=======================================

Overview
========

The ``&par`` namelist contains the tight-binding LMTO potential parameters. 
These define the electronic structure basis and are typically taken from 
first-principles databases or fitted to ab-initio calculations.

Key Parameters
==============

lmax
----

**Type:** Integer

**Purpose:** Maximum orbital angular momentum

**Allowed values:** 1 (s,p), 2 (s,p,d), 3 (s,p,d,f)

**Default:** 2

.. note::

   The basis of the run is set by ``lmax`` in ``&control``. ``lmax`` in ``&par``
   records the potential file's own value; a file with a smaller ``lmax`` is
   promoted to the basis with a warning. ``lmax = 3`` needs
   ``strux_backend = 'strux_lib'`` (:ref:`keywords/lattice_geometry`).

**Example:**

.. code-block:: fortran

   lmax = 2  ! Include s, p, d orbitals

**Meaning:**

- lmax=0: single orbital (s-like) – rare, very approximate
- lmax=1: s and p orbitals – for light elements (C, N, O)
- lmax=2: s, p, d orbitals – **standard for transition metals**
- lmax=3: includes f orbitals – for rare earths (experimental)

**Notes:**

- Larger lmax = more basis functions = larger Hamiltonian matrix
- Most production calculations use lmax=2
- Must be consistent with potential database

**Related code:** ``source/potential.f90::type potential``

center_band
-----------

**Type:** Real array (dimension: lmax+1, nsp)

**Purpose:** Energy reference (orbital center) for each orbital

**Units:** Rydberg

**Typical range:** -2 to +1 Ry

**Example:**

.. code-block:: fortran

   ! Silicon
   center_band(:,1) = -1.070, -0.250, -0.098
   center_band(:,2) = -1.070, -0.250, -0.098

**Meaning:**

- center_band(1) = s orbital reference energy
- center_band(2) = p orbital reference energy
- center_band(3) = d orbital reference energy
- For collinear: spin-up and spin-down can differ

**Physical interpretation:**

- Negative values: electrons bound to atom (typical)
- Position determines electron affinity
- Related to orbital ionization potential

**Notes:**

- Usually same for spin-up/down (unless magnetic ground state)
- Key parameter affecting position of bands
- Small changes have large effect on band structure

**Related code:** ``source/potential.f90::type potential``

width_band
----------

**Type:** Real array (dimension: lmax+1, nsp)

**Purpose:** Bandwidth parameter (√Δ) for each orbital

**Units:** Rydberg

**Typical range:** 0.5-4 Ry

**Example:**

.. code-block:: fortran

   ! Silicon
   width_band(:,1) = 1.147, 1.938, 3.513
   width_band(:,2) = 1.147, 1.938, 3.513

**Meaning:**

- Controls width and shape of band(s)
- Larger value -> broader band
- Related to electron hopping strength

**Physical interpretation:**

- Proportional to transfer integral strength
- Small: localized electrons (ionic character)
- Large: delocalized electrons (metallic character)

**Notes:**

- Typically: d > p > s in magnitude
- Critical for reproducing band structure
- Affects magnetic moment predictions

**Related code:** ``source/hamiltonian.f90::build_bulkham()``

enu
---

**Type:** Real array (dimension: lmax+1, nsp)

**Purpose:** Energy parameter for orthogonal tight-binding (ORT) potential

**Units:** Rydberg

**Related code:** ``source/potential.f90``

**Notes:**

For ORT potential formulation (alternative to center_band/width_band).
Details depend on specific potential parameterization.

c, srdel, qpar, ppar, vl
------------------------

**Type:** Real arrays

**Purpose:** Additional potential parameters for ORT formulation

**Meaning:**

- c: potential parameter
- srdel: √Δ analog
- qpar, ppar: quadrupole/octupole parameters
- vl: local potential

**Notes:**

Implementation varies; see specific potential database documentation.

For exchange-correlation functional parameters, see :doc:`exchange_correlation`.

Complete Potential Parameter Sets
==================================

Potential parameters are usually NOT entered manually. Instead:

1. **Use database files:**

.. code-block:: fortran

   &element
      symbol = 'Si1',
      atomic_number = 14,
      core = 10,
      valence = 4,
      database = '/path/to/element_database/'
   /

2. **Or pre-computed in namelist** (for testing/development):

.. code-block:: fortran

   &par
      lmax = 2,
      center_band(:,1) = ...,
      width_band(:,1) = ...,
      ws_r = 2.827
   /

Potential Database Format
=========================

Standard RS-LMTO potentials are stored as namelists:

**File:** ``element_name.nml``

**Contains:** All ``&par`` parameters for that element

**Directories:**

- ``database/`` - (location depends on installation)
- Specify via ``database=`` keyword

**Creating custom potentials:**

1. Copy existing potential close to desired element
2. Modify center_band, width_band by small amounts
3. Test convergence and band structure
4. Refine via comparison with reference calculations

Provenance
==========

Basis/potential parameters defined in:

- **Type definition:** ``source/potential.f90::type potential``
- **Reading:** ``source/potential.f90::build_from_file()``
- **Usage in Hamiltonian:** ``source/hamiltonian.f90::build_bulkham()``
- **Database interface:** ``source/element.f90``

Stored potential state (&par)
=============================

These ``&par`` keys are written by the code to ``<label>_out.nml`` and read back on
restart; they are not normally edited. Orbital arrays are indexed (l+1, spin) with
spin 1 = up, 2 = down, unless stated otherwise. The Hubbard keys of ``&par`` are on
:ref:`keywords/lda_u`.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``pl``
     - Potential parameters :math:`P_l = 0.5 - \arctan(D_l)/\pi` per l and spin.
   * - ``ql``
     - Moments of the occupied density per l and spin (Phys. Rev. B 43, 9538 (1991), Eq. 48).
   * - ``ppar``, ``qpar``, ``srdel``, ``vl``
     - Potential parameters on the orthogonal basis per l and spin.
   * - ``gravity_center``
     - Band gravity centre per l and spin.
   * - ``shifted_band``
     - ``center_band`` minus ``gravity_center``.
   * - ``obar``
     - Potential parameter :math:`\bar{o}` per l and spin; ``obx`` is built from it.
   * - ``center_band_s_up``, ``center_band_s_dw``, ``center_band_p_up``, ``center_band_p_dw``, ``center_band_d_up``, ``center_band_d_dw``, ``width_band_s_up``, ``width_band_s_dw``, ``width_band_p_up``, ``width_band_p_dw``, ``width_band_d_up``, ``width_band_d_dw``
     - Copies of ``center_band`` and ``width_band`` for s, p, d and each spin, written for
       reading; the code does not read them back.
   * - ``cx``, ``wx``, ``cex``, ``obx``
     - Internal per-orbital potential parameters (9 orbitals × 2 spins): band centre, band
       width, centre minus gravity centre, and the ``hoh`` parameter o.
   * - ``cx0``, ``cx1``, ``wx0``, ``wx1``
     - Spin sum (``0``) and difference (``1``) of ``cx`` and ``wx``.
   * - ``cshi``, ``dw_l``
     - Band variables per orbital (18 entries); no further meaning documented in the code.
   * - ``xi_p``, ``xi_d``
     - Spin-orbit coupling parameters of the p and d channels (2 values each).
   * - ``rac``
     - Racah parameters (2 values).
   * - ``lmom``
     - Orbital moment of the atom (not normalised).
   * - ``sumec``, ``sumev``
     - Core and valence band-energy sums (Ry).
   * - ``etot``, ``utot``, ``ekin``, ``rhoeps``
     - Total, electrostatic and kinetic energy and the XC energy integral of the atom (Ry).
   * - ``vmad``
     - Madelung potential at the site (Ry).
   * - ``vrmax``
     - Sphere-boundary potential for the muffin-tin zero: spin average and up minus down.

See Also
========

- :ref:`keywords/control_parameters` - txc (XC functional)
- :ref:`keywords/lattice_geometry` - ws_r (sphere radius)
- :doc:`../theory/lmto_asa_overview` - Theoretical background
- :doc:`../user_guide/examples` - Example potential values
