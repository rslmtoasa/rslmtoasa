.. _getting_started:

=======================================
Getting Started
=======================================

Overview
========

This section guides you through installing, compiling, and running the RS-LMTO-ASA code.

**RS-LMTO-ASA** is a Fortran 2008 code that performs first-principles electronic structure 
calculations using the Linear Muffin-Tin Orbital method with the Atomic Sphere Approximation. 
It supports both scalar-relativistic and fully relativistic calculations with spin-orbit coupling, 
collinear and non-collinear magnetism, and magnetic impurity calculations.

Prerequisites
=============

**Required:**

- Fortran 2008 compiler (Intel ifort, GNU gfortran, or equivalent)
- CMake 3.5 or later
- Make

**Optional:**

- OpenMP (for shared-memory parallelization, recommended)
- MPI (for distributed-memory parallelization)
- Python 3.7+ (for auxiliary scripts)

Building RS-LMTO-ASA
====================

Clone the Repository
--------------------

.. code-block:: bash

   git clone https://github.com/rslmto/rslmto_devel.git
   cd rslmto_devel

Create a Build Directory
------------------------

.. code-block:: bash

   mkdir build
   cd build

Configure with CMake
--------------------

**Default configuration (OpenMP enabled, MPI disabled):**

.. code-block:: bash

   cmake ..

**With MPI support:**

.. code-block:: bash

   cmake -DENABLE_MPI=ON ..

**With custom compiler:**

.. code-block:: bash

   cmake -DCMAKE_Fortran_COMPILER=ifort ..

**View available CMake options:**

.. code-block:: bash

   cmake -LAH ..

Key CMake Options
~~~~~~~~~~~~~~~~~

.. table::
   :align: left

   +---------------------------+-----------------------------------+---------+
   | Option                    | Description                       | Default |
   +===========================+===================================+=========+
   | ENABLE_OPENMP             | Enable OpenMP parallelization     | ON      |
   +---------------------------+-----------------------------------+---------+
   | ENABLE_MPI                | Enable MPI parallelization        | OFF     |
   +---------------------------+-----------------------------------+---------+
   | COLOR                     | Enable colored terminal output    | ON      |
   +---------------------------+-----------------------------------+---------+
   | ENABLE_FLUSH              | Flush print output to files       | OFF     |
   +---------------------------+-----------------------------------+---------+
   | RUN_REG_TESTS             | Run regression tests              | OFF     |
   +---------------------------+-----------------------------------+---------+

Compile
-------

.. code-block:: bash

   make -j 4

The executable will be created at ``build/bin/rslmto.x``.

Verify Installation
-------------------

.. code-block:: bash

   ./build/bin/rslmto.x --help

Running a Calculation
=====================

Basic Syntax
------------

The code reads from ``input.nml`` in the current directory by default. Optionally, you can 
specify a different input filename as the first argument.

**Default usage (reads input.nml):**

.. code-block:: bash

   ./rslmto.x

**Specify a different input file:**

.. code-block:: bash

   ./rslmto.x my_calculation.nml

Example Calculation
-------------------

Navigate to an example directory and copy the example file to ``input.nml``:

.. code-block:: bash

   cd example/bulk/Si
   cp Si1.nml input.nml
   
Run the calculation:

.. code-block:: bash

   ../../../build/bin/rslmto.x

Alternatively, specify the input file directly:

.. code-block:: bash

   ../../../build/bin/rslmto.x Si1.nml

With OpenMP parallelization:

.. code-block:: bash

   OMP_NUM_THREADS=4 ../../../build/bin/rslmto.x

With MPI (if compiled with MPI support):

.. code-block:: bash

   mpirun -np 4 ../../../build/bin/rslmto.x

Input and Output
================

**Input files** are specified in Fortran namelist format (``.nml``), containing all control 
parameters, structure information, and calculation settings.

See :ref:`user_guide/input_files` for detailed information on input file structure.

Output files are written to the current directory and include:

- ``*_out.nml``: Output parameters file
- Energy values and convergence information
- Band structure and density of states (if requested)
- Magnetic properties (if applicable)
- Additional analysis files (forces, exchange parameters, etc.)

See :ref:`user_guide/output_files` for details.

Environment Variables
======================

**OpenMP Control**

.. code-block:: bash

   export OMP_NUM_THREADS=4          # Number of threads
   export OMP_NESTED=TRUE            # Nested parallelism
   export OMP_STACKSIZE=256M         # Stack size per thread

**MPI Control**

.. code-block:: bash

   export OMP_NUM_THREADS=2          # Threads per MPI rank (hybrid)
   mpirun -np 8 ./rslmto.x

**I/O Tuning**

.. code-block:: bash

   export OMP_PROC_BIND=CLOSE        # CPU affinity

Troubleshooting
===============

**Compilation fails with "module not found"**

This usually indicates a build directory issue. Try cleaning and rebuilding:

.. code-block:: bash

   rm -rf build
   mkdir build && cd build
   cmake .. -DCMAKE_Fortran_COMPILER=gfortran
   make -j 4

**Executable crashes with MPI**

Ensure MPI libraries are linked correctly:

.. code-block:: bash

   cmake -DENABLE_MPI=ON -DMPI_Fortran_COMPILER=mpif90 ..

**Performance issues with OpenMP**

Check thread affinity and pinning. Use ``OMP_PROC_BIND=CLOSE`` for better locality.

Provenance
==========

This section documents the codebase locations used:

- **Compilation system**: ``CMakeLists.txt`` (root), ``source/CMakeLists.txt``
- **Entry point**: ``source/main.f90`` (program entry, MPI initialization, timer)
- **Build options**: ``cmake/SetFortranFlags.cmake`` (compiler flags)
- **Example inputs**: ``example/bulk/Si/Si1.nml``, ``example/bulk/bccFe/``, etc.

Next Steps
==========

After successful installation:

1. Read :doc:`code_structure` to understand the repository layout
2. Review :ref:`theory/lmto_asa_overview` for theoretical background
3. Check :ref:`user_guide/examples` for worked examples
4. Consult :ref:`keywords/index` for input parameter reference
