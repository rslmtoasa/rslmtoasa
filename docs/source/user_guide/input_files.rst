.. _user_guide/input_files:

=======================================
Input Files
=======================================

Overview
========

RS-LMTO-ASA uses **Fortran namelist format** for input. Each calculation requires at least one 
namelist file (typically ``.nml`` extension) containing all calculation parameters, structure 
information, and control settings.

File Format
===========

Fortran namelists use the following syntax:

.. code-block:: fortran

   &namelist_name
      parameter_1 = value_1,
      parameter_2 = value_2,
      ...
   /

**Key Rules:**

- Names are case-insensitive
- Commas separate parameters (trailing comma is OK)
- Parameters can span multiple lines
- Values can be integers, reals, strings, or arrays
- Comments start with ``!`` or ``*``

Example
-------

.. code-block:: fortran

   &control
      max_iterations = 100,
      nsp = 1,          ! Collinear scalar-relativistic
      verbose = .true.
   /
   
   &lattice
      alat = 5.4,       ! Lattice constant in Angstrom
      nbulk = 1,        ! Number of bulk atoms
      nx = 5, ny = 5, nz = 1  ! Cluster geometry
   /

Structure and Content
=====================

A complete RS-LMTO input file typically includes:

1. **Element** namelist - Atomic parameters (from database or custom)
2. **Lattice** namelist - Crystal structure and geometry
3. **Par** namelist - Potential parameters (LMTO tight-binding)
4. **Control** namelist - Calculation settings
5. **Energy** namelist - Energy mesh and Fermi level settings
6. **Self** namelist - SCF options

See :ref:`keywords/index` for complete parameter reference.

Essential Parameters
====================

**Minimum required:**

- ``nbulk`` - Number of bulk atoms in cluster
- ``alat`` - Lattice constant (Å)
- ``nsp`` - Type of calculation (1=SR collinear, 2=collinear FR, 3=NC SR, 4=NC FR)
- ``max_iterations`` - Maximum SCF iterations

**Recommended to specify:**

- ``dq_tol`` - SCF convergence criterion
- ``llsp`` - Recursion cutoff for sp electrons
- ``lld`` - Recursion cutoff for d electrons
- ``mixing`` and ``alpha`` - Density mixing parameters

Running a Calculation
=====================

**Basic command:**

.. code-block:: bash

   ./rslmto.x

**With output redirection:**

.. code-block:: bash

   ./rslmto.x > output.log 2>&1

**Parallel execution (OpenMP):**

.. code-block:: bash

   export OMP_NUM_THREADS=4
   ./rslmto.x

**Parallel execution (MPI):**

.. code-block:: bash

   mpirun -np 8 ./rslmto.x

Example Input Files
===================

The code reads ``input.nml`` (or the specified input file) and produces output files depending on ``post_processing`` setting:

- ``input_out.nml`` - Echo of all parameters actually used
- ``scf_convergence.dat`` - Iteration-by-iteration convergence data

Located in ``example/`` directory:

**Bulk Si:**

``example/bulk/Si/Si1.nml`` - Diamond structure, single iteration example

.. code-block:: bash

   cd example/bulk/Si
   ../../build/bin/rslmto.x Si1.nml

**Bulk bcc Fe (magnetic):**

``example/bulk/bccFe/`` - Iron with spin-polarized calculation

**Surface (if available):**

``example/surface/`` - Surface geometry and calculations

**Impurity:**

``example/impurity/`` - Single impurity in bulk host

**Exchange interactions:**

``example/exchange/`` - Magnetic exchange parameter calculations

**Transport:**

``example/conductivity/`` - Electronic transport calculations

Template: Minimal Working Example
=================================

Create a simple `test.nml` file for a new system:

.. code-block:: fortran

   &element
      symbol='X1',
      atomic_number=14,
      core=10,
      valence=4,
      f_core=0,
      num_quant_s=3,
      num_quant_p=3,
      num_quant_d=3
   /

   &par
      ! Tight-binding parameters (from database or fitted)
      center_band(:, 1) = -1.07, -0.25, -0.098
      center_band(:, 2) = -1.07, -0.25, -0.098
      width_band(:, 1) = 1.15, 1.94, 3.51
      width_band(:, 2) = 1.15, 1.94, 3.51
      ws_r = 2.83
      lmax = 2
   /

   &lattice
      alat = 5.4,
      nbulk = 1,
      nx = 5, ny = 5, nz = 1
   /

   &control
      nsp = 1,
      llsp = 60,
      lld = 60,
      max_iterations = 50,
      verbose = .true.
   /

   &energy
      channels_ldos = 500,
      energy_min = -2.0,
      energy_max = 1.0
   /

   &self
      mixing = 'broyden',
      alpha = 0.5
   /

Common Pitfalls
===============

**1. Missing parameters**

Error: "Parameter not found in namelist"

**Solution:** Check spelling and ensure parameter is in correct namelist.

**2. Wrong data types**

Error: "Unexpected character in integer"

**Solution:** Use correct format (integers: 5, reals: 5.0, strings: 'text')

**3. Database files not found**

Error: "Element parameters file not found"

**Solution:** Specify ``database`` path or ensure element files are in expected location.

**4. Lattice parameters incompatible**

Error: "Sphere overlap exceeds limit"

**Solution:** Check ``alat`` and Wigner-Seitz radius ``ws_r`` are consistent.

**5. SCF fails to converge**

**Solution:** 
- Reduce mixing parameter ``alpha``
- Increase Broyden history
- Check recursion cutoff (``llsp``, ``lld``)
- Reduce energy mesh spacing

Advanced Topics
===============

**Array parameters:**

Many parameters are 1D or 2D arrays. Syntax:

.. code-block:: fortran

   ! 1D array
   array_1d(:) = 1.0, 2.0, 3.0
   
   ! 2D array (row-major)
   array_2d(:, :) = 1.0, 2.0, 3.0, 4.0

**Element/potential database:**

For systems with multiple atom types, include separate namelists:

.. code-block:: fortran

   &element
      symbol='Atom1', ...
   /
   
   &par
      lmax=2, ...
   /
   
   ! Repeat for additional atom types
   &element
      symbol='Atom2', ...
   /
   
   &par
      lmax=3, ...
   /

**Geometry specification:**

Cluster geometry is built via:

- ``nx, ny, nz`` - Grid dimensions
- ``alat`` - Lattice constant
- ``r2`` - Cutoff radius for neighbors

More details in :ref:`keywords/lattice_geometry`.

Provenance
==========

Input file processing implemented in:

- ``source/namelist_generator.f90`` - Namelist utilities
- ``source/control.f90::build_from_file()`` - Control parameters
- ``source/lattice.f90::build_from_file()`` - Structure parameters
- ``source/potential.f90::build_from_file()`` - Potential parameters
- ``source/element.f90::build_from_file()`` - Element database

See Also
========

- :ref:`keywords/index` - Complete parameter listing
- :doc:`output_files` - Output file format
- :doc:`examples` - Worked examples
- :doc:`../getting_started` - Running calculations
