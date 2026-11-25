.. _examples/bulk_si:

====================================
Bulk Silicon (Diamond Structure)
====================================

**Location:** ``example/bulk/Si/``

**System:** Non-magnetic diamond cubic silicon

**Physics:** Demonstrates basic SCF convergence for a simple semiconductor with covalent bonding.

Overview
========

Silicon crystallizes in the diamond cubic structure (FCC with 2-atom basis). This example shows how to set up a periodic bulk calculation with the RS-LMTO-ASA method. The calculation uses:

- ``pre_processing = 'bravais'`` to generate the crystal structure
- ``crystal_sym = 'file'`` to read structure from ``lattice.nml``
- Periodic boundary conditions (PBC) in all directions
- Collinear treatment with spin-orbit coupling (``nsp = 2``)

Files
=====

The example directory contains:

- ``input.nml`` - Main calculation parameters
- ``lattice.nml`` - Crystal structure (input file written by cif2cell; read because ``crystal_sym = 'file'``)
- ``Si1.nml`` - Tight-binding parameters for first Si atom
- ``Si2.nml`` - Tight-binding parameters for second Si atom

Input File: input.nml
=====================

The main input file specifies all calculation parameters:

.. code-block:: fortran

   &calculation
   pre_processing = 'bravais'    ! Generate Bravais lattice
   verbose = T                   ! Print detailed output
   /
   
   &lattice
   rc = 50                       ! Cluster radius (a.u.)
   alat = 5.443702               ! Lattice constant (Å)
   crystal_sym = 'file'          ! Read structure from lattice.nml
   wav = 1.4962757958            ! Wigner-Seitz radius (Å)
   pbc = .true.                  ! Use periodic boundary conditions
   b1 = .true.                   ! PBC along lattice vector a1
   b2 = .true.                   ! PBC along lattice vector a2
   b3 = .true.                   ! PBC along lattice vector a3
   n1 = 20                       ! Repeats of the periodic cell along a1
   n2 = 20
   n3 = 20
   ntype = 2                     ! Number of atom types
   ct(:) = 2*4.0d0               ! Neighbour cut radius (Å), one per atom type
   r2 = 16.0d0                   ! Squared neighbour cut radius (Å²), r2 = ct²
   /
   
   &atoms
   database = './'               ! Atom parameters in current directory
   label(1) = 'Si1'              ! First Si atom in basis
   label(2) = 'Si2'              ! Second Si atom in basis
   /
   
   &self
   nstep = 25                    ! Maximum SCF iterations
   cold = .true.                 ! Cold start (no restart)
   /
   
   &energy
   fermi = -0.168353             ! Initial Fermi energy (Ry)
   energy_min = -1.5             ! Energy window minimum (Ry)
   energy_max = 1.0              ! Energy window maximum (Ry)
   channels_ldos = 2000          ! Energy points for DOS
   /
   
   &control
   calctype = 'B'                ! Bulk calculation
   nsp = 2                       ! Collinear, with spin-orbit coupling
   lld = 250                     ! Recursion levels
   recur = 'chebyshev'           ! Chebyshev recursion method
   /
   
   &mix
   beta = 0.1                    ! Mixing parameter (conservative)
   mixtype = 'broyden'           ! Broyden mixing scheme
   /
   
   &hamiltonian
   hoh = .false.                 ! No higher-order hopping
   /

Key Parameters Explained
------------------------

**Lattice Parameters:**

- ``rc = 50``: Cluster radius in atomic units. Should include 5-10 neighbor shells. For Si with alat ≈ 5.44 Å ≈ 10.3 a.u., rc=50 includes neighbors up to ~5 lattice constants.

- ``alat = 5.443702``: Experimental lattice constant of Si at room temperature in Ångströms.

- ``crystal_sym = 'file'``: Instructs code to read full crystal structure from ``lattice.nml``. Other accepted values: 'bcc', 'b2', 'fcc', 'hcp', 'fcc2', 'fcc3' for standard lattices.

- ``wav = 1.4962757958``: Wigner-Seitz radius in Å. If ``wav`` is not set, the code computes it from the cell volume per atom.

- ``pbc = .true.`` with ``b1 = b2 = b3 = .true.``: Enables 3D periodic boundary conditions (essential for bulk calculations).

- ``n1 = n2 = n3 = 20``: number of repeats of the primitive cell along ``a1``, ``a2``, ``a3`` that form the periodic real-space cell. The k-mesh of the k-space code is ``nk1, nk2, nk3`` in ``&reciprocal``.

- ``ntype = 2``: Two inequivalent Si atoms in diamond structure.

- ``ct(:) = 2*4.0d0``: Neighbour cut radius in Å for each atom type (4.0 Å).

- ``r2 = 16.0``: Squared neighbour cut radius in Å², ``r2 = ct**2``.

**SCF Parameters:**

- ``nstep = 25``: Maximum number of SCF iterations. Si typically converges in 10-15 iterations.

- ``cold = .true.``: Start from scratch (no restart from previous calculation).

**Energy Parameters:**

- ``fermi = -0.168353``: Initial guess for Fermi energy. Will be updated during SCF.

- ``energy_min/max``: Energy window for Green's function calculation. Should bracket all valence states.

- ``channels_ldos = 2000``: Number of energy points for DOS. More points -> smoother DOS.

**Control Parameters:**

- ``calctype = 'B'``: **Bulk** calculation type. Other options: 'S' (surface), 'I' (impurity), 'L' (layered).

- ``nsp = 2``: Collinear calculation with spin-orbit coupling (``nsp = 1`` is collinear scalar-relativistic).

- ``lld = 250``: Number of recursion levels. Higher -> better convergence but slower.

- ``recur = 'chebyshev'``: Use Chebyshev recursion (Kernel Polynomial Method). Alternatives: 'lanczos', 'block'.

**Mixing Parameters:**

- ``beta = 0.1``: Mixing parameter for SCF convergence. Lower values (0.01-0.1) -> more stable but slower. Higher values (0.3-0.5) -> faster but may oscillate.

- ``mixtype = 'broyden'``: Broyden mixing for accelerated convergence.

Lattice File: lattice.nml
==========================

The ``lattice.nml`` file defines the crystal structure. It is an input file (here written by cif2cell) and is read when ``crystal_sym = 'file'``:

.. code-block:: fortran

   &lattice
   nbulk_bulk = 2                ! Number of bulk atoms
   ntot = 2                      ! Total number of atoms in unit cell
   nbas = 2                      ! Number of basis atoms
   nrec = 2                      ! Number of recursive atoms
   
   ! Lattice vectors (FCC, fractional coordinates)
   a(:, 1) = 0.5, 0.5, 0.0       ! First lattice vector
   a(:, 2) = 0.5, 0.0, 0.5       ! Second lattice vector
   a(:, 3) = 0.0, 0.5, 0.5       ! Third lattice vector
   
   ! Atomic positions (fractional coordinates)
   crd(:, 1) = 0.0, 0.0, 0.0     ! Si1 at origin
   crd(:, 2) = 0.25, 0.25, 0.25  ! Si2 at (1/4, 1/4, 1/4)
   
   ! Atom type indices
   izp(1) = 1                    ! Atom 1 is type 1 (Si1)
   izp(2) = 2                    ! Atom 2 is type 2 (Si2)
   
   ! Atom numbering
   no(1) = 1
   no(2) = 2
   
   ! Unique atom indices (for symmetry)
   iu(1) = 1
   iu(2) = 2
   
   ! Basis indices
   ib(1) = 1
   ib(2) = 2
   
   ! Recursion indices
   irec(1) = 1
   irec(2) = 2
   /

Structure Explanation
---------------------

**Diamond Structure:**

- Diamond = FCC lattice with 2-atom basis
- Basis atoms at (0,0,0) and (1/4,1/4,1/4)
- Each Si has 4 nearest neighbors in tetrahedral coordination
- Nearest neighbor distance: a√3/4 ≈ 2.35 Å

**Lattice Vectors:**

The three vectors define an FCC unit cell:

.. code-block:: text

   a₁ = (a/2, a/2, 0)
   a₂ = (a/2, 0, a/2)
   a₃ = (0, a/2, a/2)

These are fractional coordinates; multiply by ``alat`` to get Cartesian coordinates.

**Atomic Positions:**

.. code-block:: text

   Si1: (0, 0, 0)           - at cube corner
   Si2: (a/4, a/4, a/4)     - at body center

Atom Parameters: Si1.nml
=========================

Each atom type requires a separate namelist file with tight-binding parameters:

.. code-block:: fortran

   &element
   symbol = 'Si1'                ! Atom label
   atomic_number = 14            ! Z = 14 for silicon
   core = 10                     ! [Ne] core: 10 electrons
   valence = 4                   ! 3s²3p² valence electrons
   f_core = 0                    ! No f-electrons in core
   num_quant_s = 3               ! Principal quantum number for s
   num_quant_p = 3               ! Principal quantum number for p
   num_quant_d = 3               ! Principal quantum number for d
   /
   
   &par
   lmax = 2
   ! Charge moments Q0 and Q2 (s, p, d; spin-up, spin-down)
   ql(1, :, 1) = 0.7500000000, 0.9900000000, 0.2450000000
   ql(1, :, 2) = 0.7500000000, 1.0100000000, 0.2550000000
   ql(3, :, 1) = 0.0037500000, 0.0049500000, 0.0012250000
   ql(3, :, 2) = 0.0037500000, 0.0050500000, 0.0012750000
   ! Potential functions
   pl(:, 1) = 3.88, 3.68, 3.23
   pl(:, 2) = 3.88, 3.68, 3.23
   ! Band centers and widths (Ry)
   center_band(:, 1) = -1.0701613521419, -0.24985484360886, -0.097585280247979
   center_band(:, 2) = -1.0701613521419, -0.24985484360886, -0.097585280247979
   width_band(:, 1) = 1.1472986747192, 1.9376890489319, 3.5134553841396
   width_band(:, 2) = 1.1472986747192, 1.9376890489319, 3.5134553841396
   ! Wigner-Seitz radius (a.u.)
   ws_r = 2.827551
   ! (the file also holds the ORT parameters enu, c, srdel, qpar, ppar, vl)
   /

Parameters Explained
--------------------

**Element Block:**

- ``core = 10``: Si has [Ne] = 1s²2s²2p⁶ core (10 electrons). These are frozen in the calculation.

- ``valence = 4``: Si has 4 valence electrons (3s²3p²) that participate in bonding.

**Tight-Binding Parameters:**

- ``pl(:, spin)``: Potential parameters for each orbital (s, p, d). Related to radial extent of atomic orbitals.

- ``center_band``: Energy center of each orbital band. Note s and p have negative energies (bound states).

- ``width_band``: Width of each orbital band. Related to hopping strength.

- ``ql(1, :, spin)``: Charge moments (Q₀) of the s, p and d channels for each spin.

**Wigner-Seitz Radius:**

- ``ws_r = 2.827551``: Atomic sphere radius in a.u. (1.496 Å, the same radius as ``wav`` in ``input.nml``).

**Si2.nml:** The second Si atom (Si2.nml) differs from Si1.nml in ``symbol`` and in the last digits of ``center_band`` and ``width_band`` (differences up to about 3e-11), because both are Si atoms in equivalent chemical environments (diamond symmetry).

Running the Calculation
========================

Step 1: Copy the Example to a Run Directory
-------------------------------------------

.. code-block:: bash

   mkdir -p runs
   cp -r example/bulk/Si runs/Si
   cd runs/Si

Step 2: Run the Calculation
----------------------------

The code reads ``input.nml`` from the current directory by default:

.. code-block:: bash

   ../../build/bin/rslmto.x

**Alternative:** Redirect output to a log file:

.. code-block:: bash

   ../../build/bin/rslmto.x > si_scf.log 2>&1

Step 3: Monitor Progress
-------------------------

During execution, the code prints SCF iteration information:

.. code-block:: text

   Not converged! Diff= <delta>
   Not converged! Diff= <delta>
   ...
   Converged!<delta>

``<delta>`` is the mixing difference of that iteration.

**Convergence indicators:**

- Fermi energy stabilizes (varies by < 0.001 Ry)
- Total charge equals expected value (8 electrons = 4 per Si × 2 atoms)
- Energy components stop changing

Expected Output
===============

Output Files
------------

After successful completion, the following files are created:

- ``report.out`` - Total energy, Fermi energy, moments and charge transfer
- ``Si1_out.nml`` - Converged tight-binding parameters for Si1
- ``Si2_out.nml`` - Converged tight-binding parameters for Si2
- ``totaldos.out``, ``Si1_dos.out``, ``Si2_dos.out`` - Density of states
- ``band_structure.dat`` - Band structure (only with ``post_processing = 'band_structure'``)

Convergence Metrics
-------------------

**Check Fermi energy:**

.. code-block:: bash

   grep "Fermi energy" report.out

Expected: fermi ≈ -0.168 +- 0.005 Ry

**Check total charge:**

.. code-block:: bash

   grep "Occupation at atom" report.out

Expected: Each Si should have ~4 valence electrons.

**Check SCF convergence:**

.. code-block:: bash

   grep "Total RMS Diff" report.out

The final mixing difference should be below the SCF threshold.

Physical Results
----------------

**Band Gap:**

Si is a semiconductor with an indirect band gap of ~1.1 eV (experimental). The LMTO-ASA typically reproduces this within ~10-20%.

**Density of States:**

If you plot ``totaldos.out``, you should see:

- Deep minimum (band gap) near Fermi energy
- Valence band maximum below Fermi level
- Conduction band minimum above Fermi level

**Total Energy:**

.. code-block:: bash

   grep "etot" Si1_out.nml

Typical values: etot ≈ -579 Ry per Si atom.

Performance
-----------

**Expected runtime:** 20-40 seconds on a single CPU core.

**Memory usage:** ~50-100 MB.

**Parallelization:** Use OpenMP for faster execution:

.. code-block:: bash

   export OMP_NUM_THREADS=4
   ../../build/bin/rslmto.x

Troubleshooting
===============

SCF Not Converging
------------------

If the calculation doesn't converge in 25 iterations:

1. **Reduce mixing parameter:**

   .. code-block:: fortran
   
      beta = 0.05    ! Instead of 0.1

2. **Increase iteration limit:**

   .. code-block:: fortran
   
      nstep = 50

3. **Check cluster size:** Ensure ``rc`` is large enough (at least 40 a.u. for Si).

Wrong Total Charge
------------------

If total charge deviates significantly from 8 electrons:

1. **Adjust the neighbour cut radius:**

   .. code-block:: fortran
   
      ct(:) = 2*4.5d0    ! Try different values
      r2 = 20.25d0       ! r2 = ct**2

2. **Increase cluster radius:**

   .. code-block:: fortran
   
      rc = 60

3. **Check atom parameters:** Verify Si1.nml and Si2.nml have correct ``valence = 4``.

Slow Convergence
----------------

To speed up the calculation:

1. **Reduce recursion levels:**

   .. code-block:: fortran
   
      lld = 150    ! Instead of 250

2. **Reduce energy points:**

   .. code-block:: fortran
   
      channels_ldos = 1000

3. **Use OpenMP:** Set ``OMP_NUM_THREADS`` environment variable.

Extensions
==========

Calculating Band Structure
---------------------------

To calculate bands along high-symmetry directions, set ``post_processing`` in ``&calculation``
(the k-path uses the ``&reciprocal`` and ``&kpath`` namelists):

.. code-block:: fortran

   &calculation
   pre_processing = 'bravais'
   post_processing = 'band_structure'
   /

Calculating Density of States
------------------------------

DOS is calculated by default if ``channels_ldos > 0``. Increase for smoother DOS:

.. code-block:: fortran

   channels_ldos = 5000    ! Very smooth DOS

Pressure/Volume Study
---------------------

Vary ``alat`` to study pressure dependence:

.. code-block:: bash

   for alat in 5.30 5.35 5.40 5.45 5.50; do
       sed "s/alat = .*/alat = $alat/" input.nml > input_$alat.nml
       ../../build/bin/rslmto.x input_$alat.nml
       grep "etot" Si1_out.nml >> energy_vs_volume.dat
   done

Then fit E(V) to Birch-Murnaghan equation of state to extract bulk modulus.

See Also
========

- :doc:`bulk_bccfe` - Magnetic bulk system
- :doc:`bulk_mn3sn` - Non-collinear magnetism
- :doc:`../input_files` - Complete input file reference
- :doc:`../../keywords/index` - All input keywords
- :doc:`../../theory/lmto_asa_overview` - LMTO theory background