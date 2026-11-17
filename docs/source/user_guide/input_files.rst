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

   &calculation
      verbose = .true.
   /

   &control
      nsp = 1,          ! Collinear scalar-relativistic
   /

   &lattice
      alat = 5.4,       ! Lattice constant in Angstrom
      n1 = 5, n2 = 5, n3 = 1  ! Repeats of the periodic real-space cell
   /

Structure and Content
=====================

The main input file (``input.nml`` or the file given as first command-line argument) typically includes:

1. **Calculation** namelist - ``pre_processing``, ``post_processing``, ``verbose``
2. **Lattice** namelist - Crystal structure and geometry
3. **Atoms** namelist - ``database`` folder and one ``label`` per atom type
4. **Self** namelist - SCF options
5. **Energy** namelist - Energy mesh and Fermi level settings
6. **Control** namelist - ``calctype``, ``nsp``, recursion settings
7. **Mix** namelist - Density mixing

The ``&element`` and ``&par`` namelists (atomic data and LMTO potential parameters) are not read
from the main input file. They are read per label from ``<label>.nml`` in the ``database`` folder
(for example ``Fe.nml`` for ``label(1) = 'Fe'``), or from ``<label>_out.nml`` of a previous run
when it exists.

See :ref:`keywords/index` for complete parameter reference.

Essential Parameters
====================

**Minimum required:**

- ``calctype`` (``&control``) - 'B' bulk, 'S' surface, 'I' impurity or 'L' layered; the run stops if it is not one of these
- ``pre_processing`` (``&calculation``) - Cluster construction, e.g. 'bravais'
- ``alat`` - Lattice constant (Å)
- ``crystal_sym`` - Crystal structure (or the file ``lattice.nml``, ``crystal_sym = 'file'``)
- ``ntype`` and ``&atoms`` - Number of atom types and their labels
- ``nsp`` (``&control``) - Type of calculation (1=SR collinear, 2=collinear FR, 3=NC SR, 4=NC FR)

**Recommended to specify:**

- ``conv_thr`` (``&self``) - SCF convergence criterion
- ``llsp`` - Recursion cutoff for sp electrons
- ``lld`` - Recursion cutoff for d electrons
- ``mixtype`` and ``beta`` (``&mix``) - Density mixing parameters

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

- ``<label>_out.nml`` - Per atom type, the element and potential parameters at the end of the run
- ``report.out`` - Total energy, Fermi energy, spin moments and the final SCF difference

Located in ``example/`` directory:

**Bulk Si:**

``example/bulk/Si/input.nml`` - Diamond structure (``Si1.nml`` and ``Si2.nml`` are the
per-label parameter files of the two Si sites)

.. code-block:: bash

   cd example/bulk/Si
   ../../../build/bin/rslmto.x

**Bulk bcc Fe (magnetic):**

``example/bulk/bccFe/`` - Iron with spin-polarized calculation

**Surface (if available):**

``example/surface/`` - Surface geometry and calculations

**Impurity:**

``example/impurity/`` - Single impurity in bulk host


``hubbard_radial_moments='auto'`` resolves to ``local`` for a nonzero fixed U/J (``hubbard_u_general``/``hubbard_j_general`` or ``hubbard_u_impurity``/``hubbard_j_impurity`` in ``&hamiltonian``) where ``local`` is
implemented (see ``docs/source/theory/lda_u_workflow.rst``); an explicit ``spectral`` or ``local`` is never changed. A frozen
host written in the other mode has its density matrix reset with a warning. The impurity energy routing is described in
the same page.

**Exchange interactions:**

``example/exchange/`` - Magnetic exchange parameter calculations

**Transport:**

``example/conductivity/`` - Electronic transport calculations

Template: Minimal Working Example
=================================

A new system needs two kinds of files: the main input file and one ``<label>.nml`` per atom type.
The main input file below is ``example/bulk/bccFe/input.nml``:

.. code-block:: fortran

   &calculation
      pre_processing = 'bravais'
      verbose = T
   /

   &lattice
      rc = 080
      alat = 2.86120
      crystal_sym = 'bcc'
      wav = 1.40880
      ntype = 1
      ct(1) = 3.0d0
      r2 = 9.00d0
   /

   &atoms
      database = './'
      label(1) = 'Fe'
   /

   &self
      nstep = 100
   /

   &energy
      fermi = -0.069282
      energy_min = -1.0
      energy_max = 0.5
      channels_ldos = 2500
   /

   &control
      calctype = 'B'
      nsp = 2
      lld = 21
      recur = 'block'
   /

   &mix
      beta = 0.01
      mixtype = 'broyden'
   /

The file ``Fe.nml`` in the same folder (``database = './'``) holds the ``&element`` and ``&par``
namelists of the atom type (the full ``&par`` block is in ``example/bulk/bccFe/Fe.nml``):

.. code-block:: fortran

   &element
      symbol='Fe',
      atomic_number=26,
      core=18,
      valence=8,
      f_core=0,
      num_quant_s=4,
      num_quant_p=4,
      num_quant_d=3
   /

   &par
      lmax = 2
      pl(:, 1) = 4.6656807311, 4.4101846972, 3.8747773735
      center_band(:, 1) = -0.404970091, 0.342156352, -0.213346791
      width_band(:, 1) = 0.399582709, 0.260624565, 0.117621554
      ws_r = 2.662200
      ! (abridged: ql, spin-down columns and energies are in Fe.nml)
   /

Common Pitfalls
===============

**1. Missing parameters**

Error: ``Invalid &<name> namelist; iostatus = ...``

**Solution:** Check spelling and ensure each parameter is in its own namelist (for example ``verbose`` is in ``&calculation``, not ``&control``).

**2. Wrong data types**

The Fortran runtime reports a read error in the namelist.

**Solution:** Use correct format (integers: 5, reals: 5.0, strings: 'text')

**3. Database files not found**

Error: ``Element <label> not found in any database``

**Solution:** Specify ``database`` (``&atoms``) or ensure ``<label>.nml`` is in the expected location.

**4. Lattice parameters incompatible**

Inconsistent ``alat``, ``wav`` and ``ws_r`` give strongly overlapping spheres and poor results.

**Solution:** Check ``alat`` and Wigner-Seitz radius ``ws_r`` are consistent.

**5. SCF fails to converge**

**Solution:** 
- Reduce mixing parameter ``beta`` (``&mix``)
- Check recursion cutoff (``llsp``, ``lld``)
- Reduce energy mesh spacing

Advanced Topics
===============

Reciprocal-space mode selection
-------------------------------

When using reciprocal-space post-processing (for example ``post_processing='band_structure'`` or
``post_processing='density_of_states'``), the ``&reciprocal`` namelist supports:

.. code-block:: fortran

   &reciprocal
      reciprocal_mode = 'ham_only'  ! default
      kanpur_diagnostics = .true.
      gamma_bounds_diagnostics = .false.
      hall_diag_experimental = .false.
   /

Allowed ``reciprocal_mode`` values:

- ``'ham_only'``:
  
  - Standard Hamiltonian-only k-space diagonalization.

- ``'generalized_overlap_proxy'``:
  
  - Experimental generalized-eigen workflow using a proxy overlap matrix.
  - Intended for diagnostics and method development.

- ``'generalized_overlap_kanpur'``:
  
  - Reserved for a future formal Kanpur-consistent overlap implementation.
  - Not implemented; the run stops with an error.

**Array parameters:**

Many parameters are 1D or 2D arrays. Syntax:

.. code-block:: fortran

   ! 1D array
   array_1d(:) = 1.0, 2.0, 3.0
   
   ! 2D array (column-major: first index varies fastest)
   array_2d(:, :) = 1.0, 2.0, 3.0, 4.0

**Element/potential database:**

For systems with multiple atom types, list the labels in ``&atoms`` of the main input file and
put each atom type's namelists in its own ``<label>.nml`` file:

.. code-block:: fortran

   ! Atom1.nml
   &element
      symbol='Atom1', ...
   /
   
   &par
      lmax=2, ...
   /
   
   ! Atom2.nml
   &element
      symbol='Atom2', ...
   /
   
   &par
      lmax=3, ...
   /

**Geometry specification:**

Cluster geometry is built via:

- ``crystal_sym``, ``n1, n2, n3`` - Crystal structure and number of repeats of the periodic real-space cell
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
