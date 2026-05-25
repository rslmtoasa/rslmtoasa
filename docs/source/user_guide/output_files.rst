.. _user_guide/output_files:

=======================================
Output Files
=======================================

Overview
========

After running RS-LMTO-ASA, various output files are generated depending on your calculation 
settings. This section describes the main output files and how to interpret them.

Standard Output Files
=====================

**input_out.nml**

Echo of all parameters used in the calculation. Useful for:

- Verifying that input was read correctly
- Reproducing a calculation
- Documenting what parameters were used

Format: Same as input file (Fortran namelist)

Example:

.. code-block:: fortran

   &control
      nsp = 1,
      max_iterations = 50,
      dq_tol = 0.000100,
      ...
   /

**scf_convergence.dat** (or similar)

SCF iteration-by-iteration convergence data. Columns typically:

.. code-block:: text

   Iteration    ΔQ (Ry)    ΔE (Ry)    Wall Time (s)    Fermi Energy (Ry)
   -----------  -------    -------    -----------      -----------------
   1            0.523      -0.287     2.1              -0.125
   2            0.156       0.089     1.9              -0.124
   3            0.042       0.024     1.8              -0.123
   ...

Use this to:

- Monitor convergence rate
- Detect oscillations or divergence
- Estimate total time needed
- Extract final Fermi energy

Density of States Output
========================

If ``post_processing = 'dos'`` or similar:

**dos.dat** (format: energy, total DOS, atom-resolved DOS)

Columns:

.. code-block:: text

   Energy(Ry)   Total DOS      Atom 1 DOS   Atom 2 DOS   ...
   -----------  -----------    -------      -------
   -2.000       0.00000        0.00000      0.00000
   -1.995       0.00234        0.00123      0.00111
   -1.990       0.00289        0.00145      0.00144
   ...

**Usage:**

- Plot with gnuplot, matplotlib, etc.
- Integrate to verify electron count
- Identify band structure features

**Example plot (matplotlib):**

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   data = np.loadtxt('dos.dat')
   E = data[:, 0]
   dos = data[:, 1]
   
   plt.plot(E, dos)
   plt.xlabel('Energy (Ry)')
   plt.ylabel('DOS')
   plt.show()

Band Structure Output
=====================

If ``post_processing = 'bands'``:

**bands.dat** (or similar)

Usually contains energy bands along high-symmetry k-paths.

Format depends on implementation, but typically:

.. code-block:: text

   k-point     Band 1      Band 2      Band 3      ...
   -------     ------      ------      ------
   0.000       -1.234      -0.567       0.891
   0.010       -1.230      -0.560       0.895
   ...

**Plotting:**

Use standard band structure visualization tools.

Magnetic Properties
===================

If ``nsp > 1`` (magnetic calculation):

**magnetization.dat** or similar

Atom-resolved magnetic moments:

.. code-block:: text

   Atom    m_x        m_y        m_z       |m|
   ----    -----      -----      -----     -----
   1       0.0000     0.0000     2.234     2.234
   2       0.0000     0.0000    -2.201     2.201
   ...

**Interpretation:**

- $m_z$ is the z-component of magnetization
- $|m|$ is total magnetic moment
- Negative values indicate opposite spin direction

**For non-collinear magnetism (nsp=3,4):**

All three components ($m_x, m_y, m_z$) have meaning; they define the spin direction.

Exchange Interactions
=====================

If ``post_processing = 'exchange'``:

**exchange.dat** or similar

Magnetic exchange interactions:

.. code-block:: text

   Atom1  Atom2   Distance(Å)   J (mRy)    J/2 (mRy)
   -----  -----   -----------   -------    ---------
   1      2       2.834         -1.234     -0.617
   1      3       4.006         -0.123     -0.062
   ...

**Interpretation:**

- $J > 0$ → Ferromagnetic coupling
- $J < 0$ → Antiferromagnetic coupling
- Magnitude indicates strength

**Use for:**

- Mapping of exchange networks
- Estimation of ordering temperature (Curie/Néel)
- Comparison with experiment

Transport Properties
====================

If ``post_processing = 'conductivity'``:

**conductivity_tensor.dat** or similar

Electrical conductivity tensor as function of energy:

.. code-block:: text

   Energy(Ry)   σ_xx      σ_yy      σ_zz     σ_xy (etc.)
   -----------  ------    ------    ------   ------
   -2.0         0.00      0.00      0.00     0.00
   -1.9         0.023     0.023     0.021    0.001
   ...

**Interpretation:**

- Diagonal elements: longitudinal conductivity
- Off-diagonal elements: Hall conductivity, anomalous Hall effect
- Sum over energies (with Fermi distribution) gives transport properties at temperature

Charge Density Output
=====================

If charge density output is requested:

**charge_density.dat** or **charge.dat**

Local charge density per atom or per orbital.

Format:

.. code-block:: text

   Atom    Orbital    Charge (e)
   ----    -------    ----------
   1       s          0.750
   1       p          0.990
   1       d          0.245
   2       s          0.750
   ...

Convergence Diagnostics
=======================

**Typical signs of good convergence:**

- $\Delta Q$ decreases monotonically
- Fermi energy stabilizes
- Energy change becomes very small (< $10^{-6}$ Ry)
- Usually achieved in 10-100 iterations

**Warning signs (non-convergence):**

- Charge density oscillates
- Energy increases in some iterations
- $\Delta Q$ plateaus at high value
- Forces/torques don't decrease

**Action if non-convergent:**

- Reduce mixing parameter ``alpha`` (from 0.5 to 0.3)
- Increase Broyden history length
- Increase recursion cutoff (``llsp``, ``lld``)
- Use coarser energy mesh initially, then refine
- Check Hamiltonian parameters are reasonable

Reading and Plotting Output
============================

**Python script example:**

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   # Read DOS
   dos_data = np.loadtxt('dos.dat')
   energy = dos_data[:, 0]
   dos = dos_data[:, 1]
   
   # Read convergence
   scf_data = np.loadtxt('scf_convergence.dat', skiprows=1)
   iterations = scf_data[:, 0]
   dq = scf_data[:, 1]
   
   # Plot
   fig, axes = plt.subplots(1, 2, figsize=(12, 4))
   
   axes[0].plot(energy, dos)
   axes[0].set_xlabel('Energy (Ry)')
   axes[0].set_ylabel('DOS')
   axes[0].set_title('Density of States')
   
   axes[1].semilogy(iterations, dq)
   axes[1].set_xlabel('Iteration')
   axes[1].set_ylabel('ΔQ (Ry)')
   axes[1].set_title('SCF Convergence')
   
   plt.tight_layout()
   plt.show()

Analyzing Magnetic Structure
============================

**Extract magnetization:**

.. code-block:: python

   mag_data = np.loadtxt('magnetization.dat')
   m_z = mag_data[:, 3]  # z-component
   m_mag = mag_data[:, 4] # magnitude
   
   # Check if system is ferromagnetic or antiferromagnetic
   if np.all(np.sign(m_z[:-1]) == np.sign(m_z[1:])):
       print("Ferromagnetic")
   else:
       print("Antiferromagnetic")
   
   # Total magnetic moment
   m_total = np.sum(m_mag)
   print(f"Total moment: {m_total:.3f} μB")

File Size and Storage
=====================

**Typical sizes for a single calculation:**

- ``input_out.nml`` - 1-10 KB
- ``scf_convergence.dat`` - 1-100 KB (depending on iterations)
- ``dos.dat`` - 100 KB - 1 MB (depends on energy mesh)
- ``bands.dat`` - varies widely
- Other files - 1-10 MB typical

**For multiple energy windows or temperatures:**

Output files can be significantly larger; plan disk space accordingly.

Provenance
==========

Output is generated by:

- **Main output:** ``source/calculation.f90`` (overall workflow)
- **SCF convergence:** ``source/self.f90`` (each iteration)
- **DOS output:** ``source/density_of_states.f90``
- **Band structure:** ``source/bands.f90``
- **Exchange:** ``source/exchange.f90``
- **Conductivity:** ``source/conductivity.f90``
- **Magnetic properties:** ``source/hamiltonian.f90``, related modules

See Also
========

- :ref:`user_guide/input_files` - Input file format
- :doc:`examples` - Worked examples with sample output
- :ref:`keywords/output_options` - Parameters controlling output
