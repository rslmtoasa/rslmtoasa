.. _keywords/constraints:

=============================================
Constraining Field Parameters (&constraints)
=============================================

Overview
========

The ``&constraints`` namelist activates and configures site-dependent
**constraining fields** that steer the local magnetization directions toward
prescribed reference vectors during the SCF (or exchange-interaction) calculation.
This is essential for non-collinear calculations where the self-consistent
magnetization would otherwise relax away from the desired spin configuration —
for example, when computing exchange parameters for a spin spiral, or when
modelling a metastable magnetic state.

The namelist is **optional**.  If omitted, no constraining field is applied.

.. code-block:: fortran

   &constraints
     constraints_enable      = .true.
     constraints_i_cons      = 3
     constraints_code_prefac = 1
     constraints_mom_ref(1,1) =  0.0
     constraints_mom_ref(2,1) =  0.0
     constraints_mom_ref(3,1) =  1.0
     ! ... one triplet per atom
   /

For the underlying theory see :ref:`theory/constraining_fields`.

----

Parameters
==========

constraints_enable
------------------

**Type:** Logical

**Default:** ``.false.``

**Purpose:** Master switch that enables or disables the constraining field
machinery.  When ``.false.`` the entire ``&constraints`` namelist is ignored
even if present.

**Example:**

.. code-block:: fortran

   constraints_enable = .true.

**Notes:**

- Must be set to ``.true.`` for any of the other parameters to have effect.
- Initialisation of internal arrays (PID history, :math:`\lambda_t`) is
  triggered by this flag in the ``self``, ``exchange``, and ``conductivity``
  constructors.

**Code location:** ``source/control.f90``, member ``constraints_enable``;
checked in ``source/self.f90:445``.

----

constraints_i_cons
------------------

**Type:** Integer

**Default:** ``0``

**Purpose:** Selects the algorithm used to compute the constraining field
:math:`\mathbf{B}_i^{\,\mathrm{con}}` at each SCF iteration.

**Allowed values:**

.. list-table::
   :widths: 10 60
   :header-rows: 1

   * - Value
     - Algorithm
   * - ``2``
     - **Lagrange multiplier without orthogonalization.**
       Penalty energy :math:`E_{\mathrm{con}} = \lambda_t \sum_i |\hat{\mathbf{m}}_i - \hat{\mathbf{e}}_i^{\,\mathrm{ref}}|^2`.
       Simple and robust; no constraint on the field direction.
   * - ``3``
     - **Lagrange multiplier with Gram-Schmidt orthogonalization**
       (:math:`\mathbf{B} \perp \hat{\mathbf{e}}^{\,\mathrm{ref}}`).
       Removes the component of the error vector parallel to the reference
       direction so that the constraining field is purely transverse.
       Recommended for general non-collinear calculations.
   * - ``4``
     - **PID controller** (AMN non-collinear constraints).
       Uses proportional (:math:`k_P=1.30`), integral (:math:`k_I=0.35`),
       and derivative (:math:`k_D=-0.10`) gains on the perpendicular moment
       error.  Useful when moments oscillate or converge slowly.
   * - ``5``
     - **Ma-Dudarev approach** (see [Ma2015]_).
       Direct computation of
       :math:`\mathbf{B}^{\,\mathrm{new}} = \lambda_t(\mathbf{m}_i - \hat{\mathbf{e}}_i)`
       followed by orthogonalization and field mixing.  Robust general-purpose
       method with an adaptive :math:`\lambda_t` schedule.

**Example:**

.. code-block:: fortran

   constraints_i_cons = 3   ! Lagrange + orthogonalization (recommended)

**Guidance:**

- Start with ``3`` for most non-collinear calculations.
- Use ``5`` (Ma-Dudarev) if ``3`` is slow to converge or oscillates.
- Use ``4`` (PID) when moment directions are changing rapidly across iterations.
- Value ``0`` (default) means no field is computed even if
  ``constraints_enable = .true.``; this is likely a misconfiguration.

**Code location:** ``source/include_codes/abspinlib/constrain.f90``,
module variable ``i_cons``.

.. [Ma2015] P.-W. Ma and S. L. Dudarev, *Phys. Rev. B* **91**, 054420 (2015).
   DOI: `10.1103/PhysRevB.91.054420 <https://doi.org/10.1103/PhysRevB.91.054420>`_

----

constraints_code_prefac
-----------------------

**Type:** Integer

**Default:** ``1``

**Purpose:** Integer prefactor applied to the constraining field before it
is added to the exchange-correlation potential.  Scales the overall
magnitude of :math:`\mathbf{B}_i^{\,\mathrm{con}}` entering the Hamiltonian.

**Example:**

.. code-block:: fortran

   constraints_code_prefac = 1   ! default; no additional scaling

**Notes:**

- In most cases this should be left at its default value of ``1``.
- The value is stored as ``cfd_prefac`` inside the ``cfd`` module and is
  used for unit-conversion output (e.g. when printing field values in Tesla).
- Changing this parameter does **not** affect the iterative update of
  :math:`\lambda_t`.

**Code location:** ``source/include_codes/abspinlib/constrain.f90``,
module variable ``cfd_prefac``.

----

constraints_mom_ref
-------------------

**Type:** Real array, shape ``(3, natoms)``

**Default:** Not allocated (reference moments are taken from the initial
spin configuration if not provided)

**Purpose:** Reference magnetization vectors :math:`\hat{\mathbf{e}}_i^{\,\mathrm{ref}}`
toward which the local moments are constrained.  Each column ``(:, i)``
gives the :math:`(x, y, z)` components of the target moment direction
for atom :math:`i`.

The vectors need not be unit vectors; the code normalises them internally.

**Example (two-atom unit cell, moments along z and along x):**

.. code-block:: fortran

   &constraints
     constraints_enable = .true.
     constraints_i_cons = 3
     ! Atom 1: moment along +z
     constraints_mom_ref(1,1) = 0.0
     constraints_mom_ref(2,1) = 0.0
     constraints_mom_ref(3,1) = 1.0
     ! Atom 2: moment along +x
     constraints_mom_ref(1,2) = 1.0
     constraints_mom_ref(2,2) = 0.0
     constraints_mom_ref(3,2) = 0.0
   /

**Notes:**

- The number of columns must equal the number of symmetry-inequivalent
  atoms ``nrec`` (not the full cluster size).
- If ``constraints_mom_ref`` is not set, the code falls back to using the
  moments from the initial self-consistent solution (``mom0``) as references.
- For collinear calculations it is sufficient to set only the
  :math:`z`-component.

**Code location:** ``source/control.f90``, member ``constraints_mom_ref``;
read from namelist at lines 340–360.

----

constraints_bfield
------------------

**Type:** Real array, shape ``(3, natoms)``

**Default:** Not allocated (constraining field initialized to zero)

**Purpose:** Optional initial values for the constraining field
:math:`\mathbf{B}_i^{\,\mathrm{con}}` at the start of the calculation (in
internal Rydberg units).  Providing a good starting field can accelerate
convergence when restarting from a previous run.

**Example:**

.. code-block:: fortran

   constraints_bfield(1,1) = 0.0
   constraints_bfield(2,1) = 0.0
   constraints_bfield(3,1) = 0.05   ! small initial field along z on atom 1

**Notes:**

- In most cases this can be omitted (the code starts from
  :math:`\mathbf{B}_i^{\,\mathrm{con}} = 0`).
- Useful when the calculation is restarted from a partially converged
  configuration and a meaningful estimate of the field is available.

**Code location:** ``source/control.f90``, member ``constraints_bfield``.

----

Module-Level Constants (not user-settable)
==========================================

The following parameters are hard-coded in the ``cfd`` module and are
documented here for reference:

.. list-table::
   :widths: 25 15 45
   :header-rows: 1

   * - Variable
     - Value
     - Meaning
   * - ``lambda``
     - 10
     - Maximum penalty strength :math:`\lambda_{\max}`
   * - ``lambda_t``
     - ``lambda`` (init.)
     - Current (adaptive) penalty strength :math:`\lambda_t`
   * - ``induced_mom_thresh``
     - 0.5
     - Minimum :math:`|\mathbf{m}_i|` below which the constraint is skipped
       (avoids constraining induced or near-zero moments)
   * - ``bfield_beta``
     - 1.0
     - Field mixing fraction :math:`\beta` in the Ma-Dudarev step

----

Minimal Working Example
=======================

Constraining all moments in a 2-atom bcc-Fe unit cell along the global
:math:`z`-axis using the orthogonal Lagrange method:

.. code-block:: fortran

   &constraints
     constraints_enable   = .true.
     constraints_i_cons   = 3
     constraints_code_prefac = 1
     constraints_mom_ref(1,1) = 0.0
     constraints_mom_ref(2,1) = 0.0
     constraints_mom_ref(3,1) = 1.0
     constraints_mom_ref(1,2) = 0.0
     constraints_mom_ref(2,2) = 0.0
     constraints_mom_ref(3,2) = 1.0
   /

Place this block at the end of the standard ``.nml`` input file alongside the
``&control``, ``&lattice``, etc. namelists.

----

Convergence Tips
================

- If the SCF oscillates after adding constraints, reduce the overall mixing
  parameter ``alpha`` (in ``&self``) to 0.2–0.3.
- Monitor the ``lambda_t`` printed to stdout (methods 4 and 5): it should
  grow toward ``lambda_max = 10`` and stabilise.
- If a moment never reaches its reference direction, check that
  ``constraints_mom_ref`` is correctly indexed (Fortran column-major order:
  first index is the Cartesian component, second is the atom index).
- For heavy elements with spin-orbit coupling, orthogonalized methods (3 or 5)
  are strongly preferred over method 2.

----

Provenance
==========

- **Namelist definition:** ``source/include_codes/namelists/constraints.f90``
- **Control type:** ``source/control.f90`` (lines 192–196, 340–360, 418–420)
- **Algorithm implementation:** ``source/include_codes/abspinlib/constrain.f90``
- **SCF integration:** ``source/self.f90`` (lines 445–449, 862–881)
- **Exchange integration:** ``source/exchange.f90`` (lines 120–124, 1471–1492)
- **Transport integration:** ``source/conductivity.f90`` (lines 102–106)

See Also
========

- :ref:`theory/constraining_fields` — Mathematical derivations and algorithm details
- :ref:`theory/scf_cycle` — SCF loop where the constraint is applied
- :ref:`keywords/scf_settings` — SCF mixing and convergence parameters
- :ref:`keywords/control_parameters` — Related: ``nsp`` (spin polarisation type)
