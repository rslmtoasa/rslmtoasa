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
   | 'dos'              | Calculate density of states           |
   +--------------------+---------------------------------------+
   | 'bands'            | Calculate band structure              |
   +--------------------+---------------------------------------+
   | 'exchange'         | Calculate exchange interactions       |
   +--------------------+---------------------------------------+
   | 'conductivity'     | Calculate transport properties        |
   +--------------------+---------------------------------------+

**Default:** 'none'

**Example:**

.. code-block:: fortran

   post_processing = 'exchange'

**Where to find:**

- In ``&calculation`` namelist

**Related code:** ``source/calculation.f90::process()`` - selects post-processing routine

**Output files generated:**

- 'dos': ``dos.dat``, possibly ``ldos_*.dat``
- 'bands': ``bands.dat``
- 'exchange': ``exchange.dat``, ``exchange_expanded.dat``
- 'conductivity': ``conductivity_tensor.dat``, ``conductivity_scalar.dat``

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

- Long ``channels_ldos`` → large DOS and LDOS files
- Long recursion cutoff → slower but not larger files
- Large cluster → more exchange pairs

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
