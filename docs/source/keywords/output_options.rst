.. _keywords/output_options:

=======================================
Output Control Parameters
=======================================

Overview
========

Several parameters control what output files are generated and their verbosity level.

Parameters
==========

verbose
-------

**Type:** Logical

**Purpose:** Enable verbose output during calculation

**Default:** .false.

**Example:**

.. code-block:: fortran

   verbose = .true.

**Effect:**

- Prints detailed information after each setup phase
- Calls ``print_state()`` routines for structure, control, etc.
- Useful for debugging or understanding parameter values
- Increases screen output significantly

**Where to find:**

- In ``&control`` or ``&calculation`` namelist

**Related code:** ``source/calculation.f90::pre_processing_*()``

idos
----

**Type:** Integer

**Purpose:** Local density of states (LDOS) output option

**Allowed values:**

.. table::
   :align: left

   +-----+------------------+--------------------------------+
   | idos| Output           | Description                    |
   +=====+==================+================================+
   | 0   | None             | No LDOS output (default)       |
   +-----+------------------+--------------------------------+
   | 1   | First atom type  | LDOS for s, p, d of first type |
   +-----+------------------+--------------------------------+
   | 2   | Per atom type    | LDOS for s, p, d of each type  |
   +-----+------------------+--------------------------------+

**Default:** 0

**Example:**

.. code-block:: fortran

   idos = 2  ! Output LDOS per atom type

**Output files:**

With idos=1 or idos=2: creates files like ``ldos_atom1_s.dat``, ``ldos_atom1_p.dat``, etc.

Format: Two columns

.. code-block:: text

   Energy(Ry)   LDOS(states/Ry)
   -2.0         0.0234
   -1.99        0.0245
   ...

**Where to find:**

- In ``&control`` or ``&dos`` namelist

**Related code:** ``source/density_of_states.f90``

post_processing
---------------

**Type:** Character string

**Purpose:** Type of post-processing/property calculation

**Allowed values:**

.. table::
   :align: left

   +--------------------+---------------------------------------+
   | post_processing    | Action                                |
   +====================+=======================================+
   | 'none'             | No post-processing (default)          |
   +--------------------+---------------------------------------+
   | 'density_of_states'| Calculate density of states           |
   +--------------------+---------------------------------------+
   | 'band_structure'   | Calculate band structure              |
   +--------------------+---------------------------------------+
   | 'paoflow2rs'       | PAOFLOW-to-RS conversion workflow     |
   +--------------------+---------------------------------------+
   | 'exchange'         | Calculate exchange interactions       |
   +--------------------+---------------------------------------+
   | 'conductivity'     | Calculate transport properties        |
   +--------------------+---------------------------------------+
   | 'exchange_p2rs'    | Exchange in p2rs workflow             |
   +--------------------+---------------------------------------+
   | 'conductivity_p2rs'| Conductivity in p2rs workflow         |
   +--------------------+---------------------------------------+
   | 'orbital_modern'   | Orbital-moment post-processing        |
   +--------------------+---------------------------------------+
   | 'bsf'              | Bloch spectral function A(k,E) (B3)   |
   +--------------------+---------------------------------------+
   | 'kspace_green'     | k-space Green's-function driver (B2)  |
   +--------------------+---------------------------------------+
   | 'frozen_magnon'    | Frozen-magnon E(q)/omega(q) sweep (B1)|
   +--------------------+---------------------------------------+

**Default:** 'none'

**Example:**

.. code-block:: fortran

   post_processing = 'exchange'

**Where to find:**

- In ``&calculation`` namelist

**Related code:** ``source/calculation.f90::process()`` - selects post-processing routine

**Output files generated:**

- 'density_of_states': ``dos_kspace.dat`` (k-space mode) and/or DOS outputs
- 'band_structure': ``band_structure.dat``
- 'exchange': ``exchange.dat``, ``exchange_expanded.dat``
- 'conductivity': ``conductivity_tensor.dat``, ``conductivity_scalar.dat``
- 'bsf': ``bsf.dat`` (see ``source/reciprocal_bsf.f90::calculate_bsf()``)
- 'frozen_magnon': the file named by ``&frozen_magnon output_file``
  (default ``frozen_magnon.dat``), plus ``frozen_magnon_branches.dat``/
  ``frozen_magnon_modes.dat`` when ``branch_mode='auto'`` — see
  :ref:`keywords/frozen_magnon`.

gf_route
--------

**Type:** Character string

**Purpose:** Selects which engine fills the ``green%gij``/``gji`` (and
``gij_eta`` torque-ladder) arrays consumed by exchange, damping, and
conductivity post-processing (B2/B5) — a route-agnostic dispatch key so
those consumers do not need to know which engine produced the Green's
function.

**Allowed values:** ``'recursion'`` (real-space Haydock/Chebyshev, the
original route), ``'lehmann'`` (k-space, strict Lehmann representation,
:math:`\Sigma=0`, B2 backend E), ``'dyson'`` (k-space, direct Dyson
inversion :math:`[zS-H(k)-\Sigma]^{-1}`, B2 backend D)

**Default:** ``'recursion'``

**Example:**

.. code-block:: fortran

   gf_route = 'lehmann'

**Notes:**

- In ``&calculation`` namelist.
- At ``gf_route='recursion'`` this is a pure no-op relative to pre-B2
  behavior — bit-identical output.
- At ``'lehmann'``/``'dyson'`` the exchange/damping/conductivity pipeline
  is unchanged; only the Green's-function producer differs. See
  ``green_backend``/``green_eta`` (:ref:`keywords/hamiltonian_parameters`)
  for the k-space engine's own settings, and
  ``docs/dev/plans/B2_reciprocal_green.md`` /
  ``docs/dev/route_agnostic_estimators.md`` for the accepted convergence
  envelopes (Lehmann-route J_ij/damping are broadening-defined at the
  documented ``green_eta`` default, not identical to the recursion route).

**Related code:** ``source/calculation.f90::check_gf_route()``

do_damping
----------

**Type:** Logical

**Purpose:** Enable Gilbert-damping (:math:`\alpha`) evaluation in the
exchange post-processing pipeline (B5.3), via the SOC-derivative
torque-correlation route already present in the ``gij_eta`` ladder.

**Default:** ``.false.``

**Example:**

.. code-block:: fortran

   do_damping = .true.

**Notes:**

- In ``&calculation`` namelist.
- Runs on whichever Green's function ``gf_route`` filled — recursion,
  Lehmann, or Dyson — with no separate wiring per route.

**Related code:** ``source/calculation.f90`` (search ``do_damping``)

green_backend
-------------

**Type:** Character string

**Purpose:** Selects the k-space engine ``gf_route='lehmann'``/``'dyson'``
dispatches to when filling ``green%gij``/``gji`` (B2, ``reciprocal_green.f90``).

**Allowed values:** ``'lehmann'`` (strict Lehmann representation,
:math:`\Sigma=0`; one diagonalization per k-point, all energies and pairs
amortized), ``'dyson'`` (direct Dyson inversion
:math:`[zS-H(k)-\Sigma(z)]^{-1}` per (k, z); required once :math:`\Sigma \neq 0`,
e.g. a future CPA/DMFT provider)

**Default:** ``'lehmann'``

**Example:**

.. code-block:: fortran

   green_backend = 'dyson'

**Notes:**

- In ``&reciprocal`` namelist (not ``&calculation`` — ``gf_route`` is set
  in ``&calculation`` and internally assigns ``reciprocal%green_backend``
  to match; see ``source/calculation.f90`` lines ~1217, ~1281).
- With :math:`\Sigma=0`, ``'dyson'`` reproduces ``'lehmann'`` to solver
  tolerance — a permanent consistency check (gate G-B2-1, see
  ``docs/dev/plans/B2_reciprocal_green.md``).

**Related code:** ``source/reciprocal_green.f90``, ``source/lehmann_kernel.f90``,
``source/dyson_kernel.f90``

green_eta
---------

**Type:** Real (Ry)

**Purpose:** Retarded broadening :math:`\eta` for the k-space Green's
function :math:`z = E + i\eta` (B2). Larger :math:`\eta` gives faster
k-mesh convergence of intersite quantities (J_ij, damping) at the cost of
a broadening-defined (not recursion-matching) estimator.

**Default:** ``0.02`` Ry

**Example:**

.. code-block:: fortran

   green_eta = 0.02

**Notes:**

- In ``&reciprocal`` namelist.
- ``0.02`` Ry is the accepted working point from gate **G-B2-2**
  (signed 2026-07-16): documented as giving ~1% k-mesh convergence of
  J_ij at a 16³ mesh. Supersedes an earlier ``0.01`` Ry placeholder from
  gate G-B2-1. See ``docs/dev/B2_GATE_G-B2-2.md`` and
  ``docs/dev/reciprocal_green_convergence.md`` for the convergence study
  behind this number.
- The :math:`\eta \to 0`, :math:`N_k \to \infty` limit is deliberately not
  the shipped default — it is expensive and currently gated by the
  backend-E full-unreduced-BZ requirement (no symmetry reduction / k-mesh
  parallelism yet; attaches with B4, not started).

**Related code:** ``source/reciprocal_lifecycle.f90`` (``restore_to_default``)

pre_processing
---------------

**Type:** Character string

**Purpose:** Pre-processing/geometry generation

**Allowed values:**

.. table::
   :align: left

   +--------------------+---------------------------------------+
   | pre_processing     | Action                                |
   +====================+=======================================+
   | 'none'             | No pre-processing (default)           |
   +--------------------+---------------------------------------+
   | 'bravais'          | Generate bulk cluster                 |
   +--------------------+---------------------------------------+
   | 'buildsurf'        | Generate surface from bulk            |
   +--------------------+---------------------------------------+
   | 'newclusurf'       | Insert impurity into surface          |
   +--------------------+---------------------------------------+
   | 'newclubulk'       | Insert impurity into bulk             |
   +--------------------+---------------------------------------+

**Default:** 'none'

**Example:**

.. code-block:: fortran

   pre_processing = 'buildsurf'

**Where to find:**

- In ``&calculation`` namelist

**Related code:** ``source/calculation.f90::pre_processing_*()``

processing
-----------

**Type:** Character string

**Purpose:** Main calculation type

**Allowed values:**

.. table::
   :align: left

   +--------------------+---------------------------------------+
   | processing         | Type                                  |
   +====================+=======================================+
   | 'none'             | SCF only (default)                    |
   +--------------------+---------------------------------------+
   | 'sd'               | Spin dynamics                         |
   +--------------------+---------------------------------------+

**Default:** 'none'

**Where to find:**

- In ``&calculation`` namelist

**Related code:** ``source/calculation.f90::processing_*()``

Typical Output File Workflow
=============================

**Example 1: Simple SCF + DOS**

.. code-block:: fortran

   &calculation
      post_processing = 'dos'
      verbose = .true.
   /

Creates:

- ``input_out.nml`` - Parameter echo
- ``scf_convergence.dat`` - Iteration data
- ``dos.dat`` - Density of states

**Example 2: Bulk + Surface + Impurity**

.. code-block:: fortran

   &calculation
      pre_processing = 'buildsurf'
   /
   ! (then run once to build surface)

   &calculation
      pre_processing = 'newclusurf'
      post_processing = 'exchange'
   /
   ! (run again on impurity system)

**Example 3: Complete magnetic study**

.. code-block:: fortran

   &control
      nsp = 2  ! Collinear FR with SOC
   /
   &calculation
      post_processing = 'exchange'
      verbose = .true.
   /
   &output
      idos = 2
   /

Creates:

- SCF convergence data
- Exchange interactions
- LDOS files
- Magnetic moments

Controlling Output Size
=======================

**Large output:**

- Long ``channels_ldos`` -> large DOS and LDOS files
- Long recursion cutoff -> slower but not larger files
- Large cluster -> more exchange pairs

**Disk space:**

Typical example calculations:

- Si (5×5×5 cluster, dos): ~100 KB - 1 MB
- Fe with exchange (7×7×7): ~1-5 MB
- Multiple energy windows: scale with number of runs

**Managing output:**

.. code-block:: bash

   # Compress after calculation
   gzip dos.dat
   gzip exchange.dat
   
   # View compressed files
   gunzip -c dos.dat.gz | head -20

Parsing Output Files
====================

**Python example:**

.. code-block:: python

   import numpy as np
   
   # Read DOS
   dos = np.loadtxt('dos.dat')
   E = dos[:, 0]
   rho = dos[:, 1]
   
   # Read convergence
   conv = np.loadtxt('scf_convergence.dat', skiprows=1)
   iters = conv[:, 0]
   dq = conv[:, 1]
   
   # Read exchange
   ex = np.loadtxt('exchange.dat', skiprows=1)
   pairs = ex[:, :2].astype(int)
   J_values = ex[:, 3]

Provenance
==========

Output control parameters found in:

- **verbose:** ``source/control.f90::type control``
- **idos:** ``source/control.f90::type control``
- **post_processing:** ``source/calculation.f90::type calculation``
- **pre_processing:** ``source/calculation.f90::type calculation``
- **Output writers:**
  
  - ``source/density_of_states.f90`` (DOS files)
  - ``source/exchange.f90`` (exchange files)
  - ``source/conductivity.f90`` (conductivity files)
  - ``source/self.f90`` (convergence file)

See Also
========

- :ref:`keywords/control_parameters` - Related control options
- :doc:`../user_guide/output_files` - Output file formats and interpretation
- :doc:`../user_guide/examples` - Example output configurations
