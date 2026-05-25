.. _theory/lmto_asa_overview:

=======================================
LMTO-ASA Formalism Overview
=======================================

Introduction
============

The Linear Muffin-Tin Orbital (LMTO) method is a first-principles electronic structure method
based on density functional theory (DFT) that provides an efficient, physical description
of crystal band structures and local properties.

The **Atomic Sphere Approximation (ASA)** is the simplification at the heart of RS-LMTO,
making the method fast and scalable while maintaining good accuracy for bulk materials,
surfaces, and clusters.

Key Concepts
============

Linear Muffin-Tin Orbital (LMTO) Basis
--------------------------------------

The LMTO method expands the electronic wavefunction in terms of *muffin-tin orbitals*:

.. math::

   \psi(\mathbf{r}) = \sum_{R,n} c_{R,n} \chi_{R,n}(\mathbf{r})

where:

- $\mathbf{R}$ indexes atomic sites
- $n$ indexes orbital quantum numbers (s, p, d)
- $\chi_{R,n}(\mathbf{r})$ are the muffin-tin orbitals (MTOs)

**Each MTO** consists of:

1. A solution to the radial Schrödinger equation inside the **atomic sphere** (Wigner-Seitz)
2. An energy-dependent *envelope function* in the **interstitial** region

This hybrid construction ensures:

- Variational flexibility (energy-dependent)
- Physical transparency (atomic character inside sphere)
- Computational efficiency (no planewave expansions)

Atomic Sphere Approximation (ASA)
---------------------------------

The ASA partitions real space into **overlapping atomic spheres**, one per atom, such that
the sum of sphere volumes slightly exceeds the total crystal volume. Inside each sphere,
the potential and kinetic energy operator are *assumed to be spherically symmetric about the atomic nucleus*.

**Advantages:**

- Reduces to solving radial differential equations (1D instead of 3D)
- Potential represented by a few parameters per atom
- Compatible with efficient Green's function methods
- Good for **bulk materials**, **surfaces**, and **clusters**

**Accuracy:**

- Excellent for s, p, d band metals and semiconductors
- Less accurate for systems with significant d-character mixing or f-electrons
- Well-suited for magnetic systems

Hamiltonian in LMTO
-------------------

In the LMTO basis, the effective Hamiltonian matrix elements are constructed from:

1. **One-electron energies** (orbital center $C_l$ and band width $\sqrt{\Delta_l}$)
2. **Madelung/Coulomb corrections** (monopole, dipole, etc.)
3. **Spin-orbit coupling** (if included)
4. **Exchange-correlation potential** (LDA, GGA, etc.)

The **Hamiltonian matrix** $H_{ij}$ includes:

- **On-site (diagonal)** terms: orbital energies and local potential
- **Hopping (off-diagonal)** terms: LMTO hopping integrals between neighboring atoms

**Tight-Binding Representation:**

The LMTO Hamiltonian can be cast as a tight-binding model:

.. math::

   H_{ij} = \epsilon_i \delta_{ij} + t_{ij} (1 - \delta_{ij})

where hopping integrals $t_{ij}$ depend on inter-atomic distances and directions.

ASA and Potential Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Within the ASA, the effective crystal potential is represented by:

- **Spherical part** inside each atomic sphere (solved for radial wavefunctions)
- **Constant part** in the interstitial region
- **Madelung corrections** for charge neutrality and long-range effects

The potential is characterized by few parameters:

- **Orbital centers** $C_l$ (energy reference for each $l$)
- **Band widths** $\sqrt{\Delta_l}$ (bandwidth)
- **Wigner-Seitz radius** $s$ (defines sphere size)

These are pre-computed or fitted to first-principles or ab-initio data.

Green's Function Formalism
--------------------------

Instead of diagonalizing the Hamiltonian, RS-LMTO computes the **single-particle Green's function**:

.. math::

   G_{ij}(E) = \langle i | (E - H + i\eta)^{-1} | j \rangle

where $E$ is energy and $\eta$ is a broadening parameter.

**Key Properties:**

- Encodes all spectral information (eigenvalues, eigenvectors)
- Efficient recursion schemes available (see :ref:`theory/recursion_method`)
- Real-space formulation: direct access to local and inter-site quantities
- Imaginary part gives density of states: $\rho(E) = -\frac{1}{\pi} \text{Im} G(E)$

Charge Density and Self-Consistency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The charge density is obtained by integrating the imaginary part of the Green's function:

.. math::

   n_i = \int_{-\infty}^{E_F} dE \, \text{Im} G_{ii}(E) \, f(E)

where $E_F$ is the Fermi energy and $f(E)$ is the Fermi-Dirac distribution.

The self-consistent field (SCF) procedure:

1. Start with initial charge density
2. Construct Hamiltonian using exchange-correlation potential
3. Calculate new charge density from Green's function
4. Mix old and new densities
5. Check convergence; repeat if needed

Reciprocal vs. Real-Space Representation
-----------------------------------------

**Traditional k-space LMTO:**

- Computes electronic bands $E_n(\mathbf{k})$ on a k-point mesh
- Requires Brillouin zone integration
- Best for perfect periodic crystals

**Real-space LMTO (RS-LMTO):**

- Works directly in real space with cluster geometry
- Avoids k-sampling; efficient for disorder, surfaces, impurities
- Green's function encodes all spatial and energetic information
- On-site and inter-site Green's functions directly computed

Collinear vs. Non-Collinear Magnetism
--------------------------------------

**Collinear (NSP = 1,2):**

- Spin-up and spin-down electrons handled separately
- Simpler magnetic structures (ferromagnetism, antiferromagnetism)
- Exchange splitting enters via different potentials for each spin

**Non-collinear (NSP = 3,4):**

- Spins can point in arbitrary directions (canted, helical, etc.)
- Requires vector spin density representation
- More computationally expensive; richer physics

Relativistic Effects
--------------------

RS-LMTO supports:

1. **Scalar relativistic (SR)** - Darwin and mass-velocity corrections
2. **Fully relativistic (FR)** - Dirac equation; includes spin-orbit coupling
3. **Spin-orbit coupling (SOC)** - Off-diagonal terms in Hamiltonian relating spin and orbital angular momentum

For heavy elements (Au, W, etc.), relativistic effects are crucial.

Typical Workflow
================

.. code-block:: text

   1. Define crystal structure → lattice.f90
   2. Load element and potential parameters → element.f90, potential.f90
   3. Build Hamiltonian matrix → hamiltonian.f90
   4. Compute recursion coefficients → recursion.f90
   5. Calculate Green's functions → green.f90
   6. Integrate DOS for charge density → density_of_states.f90, charge.f90
   7. Update potential and check convergence → self.f90
   8. Repeat 3-7 until convergence
   9. Calculate properties (DOS, bands, exchange, etc.) → various modules

Advantages and Limitations
===========================

**Advantages:**

✓ Efficient for large clusters (scaling better than k-space methods)
✓ Natural treatment of disorder, surfaces, impurities
✓ Direct access to real-space quantities (local moments, magnetization density)
✓ Strong analytical and semi-empirical foundation
✓ Excellent for transition metals and magnetism

**Limitations:**

✗ ASA less accurate for systems with rough potentials (f-electrons)
✗ Cluster size effects (need to reach bulk limit)
✗ No systematic way to improve beyond ASA
✗ Fewer automated fitting procedures than plane-wave codes

Provenance
==========

This section is based on:

- **Main code files:**
  
  - ``source/calculation.f90`` - Calculation driver and workflow
  - ``source/hamiltonian.f90`` - Hamiltonian construction
  - ``source/green.f90`` - Green's function calculations
  - ``source/recursion.f90`` - Recursion method
  - ``source/self.f90`` - SCF loop

- **Reference material:**
  
  - Andersen & Jepsen, PRL 53, 2571 (1984) - Original LMTO paper
  - Andersen, Phys. Rev. B 12, 3060 (1975) - ASA formalism
  - Skriver, The LMTO Method (Springer, 1992) - Comprehensive review

See Also
========

- :ref:`theory/green_functions` - Detailed Green's function formalism
- :ref:`theory/recursion_method` - Efficient recursion algorithms
- :ref:`theory/scf_cycle` - Self-consistent field procedure
- :ref:`reference/module_overview` - Code structure and modules
