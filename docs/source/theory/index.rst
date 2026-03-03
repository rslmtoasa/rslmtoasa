.. _theory_index:

=======================================
Theory and Methods
=======================================

This section provides an overview of the theoretical foundations of the RS-LMTO-ASA method,
including the Linear Muffin-Tin Orbital (LMTO) formalism, the Atomic Sphere Approximation (ASA),
and the computational techniques employed.

Table of Contents
=================

.. toctree::
   :maxdepth: 2

   lmto_asa_overview
   green_functions
   recursion_method
   scf_cycle
   spin_dynamics
   constraining_fields

Overview
========

**RS-LMTO-ASA** is a density functional theory (DFT) code based on:

1. **Linear Muffin-Tin Orbital (LMTO) Method**
   
   - Tight-binding expansion in the local orbital basis
   - Canonical transformation to atomic-like Bloch sums
   - Efficient representation of band structure

2. **Atomic Sphere Approximation (ASA)**
   
   - Partitions real space into overlapping atomic spheres
   - Simplification of the potential and kinetic energy operator
   - Significant computational advantages for bulk and cluster calculations

3. **Real-Space Green's Function Approach**
   
   - Avoids Brillouin zone sampling
   - Direct access to local and inter-site properties
   - Efficient for disorder and surface calculations

4. **Recursion Method**
   
   - Lanczos recursion for Green's function moments
   - Chebyshev expansion for spectral quantities
   - Linear scaling in system size (within clusters)

Key Features
============

- **Spin-orbit coupling** and fully relativistic treatments
- **Non-collinear** and **collinear** magnetic calculations
- **Electronic transport** properties (conductivity, magnetotransport)
- **Magnetic impurities** and defects
- **Hybrid OpenMP/MPI parallelization** for HPC systems
- **Density mixing** strategies (linear, Broyden)

Further Reading
===============

For detailed derivations and applications, see the references listed in each theory section.

Connection to Code
==================

Each theory subsection includes a **Provenance** section mapping to relevant source files:

- ``source/calculation.f90`` - Main calculation driver
- ``source/self.f90`` - SCF loop implementation
- ``source/green.f90`` - Green's function calculations
- ``source/recursion.f90`` - Recursion method implementation
- And others...

These map textbook theory to working Fortran code.
