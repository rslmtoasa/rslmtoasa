.. _code_features:

=======================================
Code Features
=======================================

This page is an inventory of what the code can do — the computational
frameworks it offers, the geometries it can build, and the physical
quantities it can compute — without going into the theory behind any of
them. See :ref:`theory_index` for the underlying formalism, and
:ref:`keywords/index` for the exact namelist keys that select each option.

Computational frameworks
=========================

The Hamiltonian is always assembled in a real-space tight-binding LMTO-ASA
basis; from there, three distinct methods are available for evaluating the
electronic structure, selected via ``control%recur``, plus a separate
diagonalization-based k-space route for periodic systems.

Recursion methods (real space)
-------------------------------

- **Lanczos recursion** (``recur = 'lanczos'``) — the original scalar
  recursion method, building a continued-fraction Green's function from a
  tridiagonal chain of recursion coefficients. One chain per orbital.
- **Block (Haydock) recursion** (``recur = 'block'``) — recursion on
  blocks of orbitals simultaneously rather than one at a time, needed for
  Hamiltonians where individual orbital channels are not decoupled (e.g.
  spin-orbit coupling, non-collinear magnetism).
- **Chebyshev / kernel polynomial method (KPM)** (``recur = 'chebyshev'``)
  — expands the Green's function in Chebyshev polynomials of the
  Hamiltonian rather than a continued fraction. Has dedicated fast CPU
  kernels (``chebyshev_fast.f90``), an optional Intel MKL-accelerated path
  (batched GEMM and inspector-executor sparse, ``ENABLE_MKL_KERNELS``),
  and an optional CUDA plugin (``ENABLE_CUDA_PLUGIN``) for GPU offload of
  the same kernels. This is also the route the Chebyshev moment
  formulation of the conductivity calculation and the Kernel Polynomial
  Method transport machinery build on.

All three recursion routes are energy-independent per site (the expensive
part — building the recursion coefficients — is done once and reused
across the whole energy mesh), which is what makes the real-space route
practical for the large, low-symmetry clusters (impurities, surfaces,
interfaces) the code is built around.

k-space route
-------------

For periodic bulk systems, a separate diagonalization-based path
(``reciprocal.f90`` and its submodules) builds :math:`H(k)` by Fourier
transform of the real-space hoppings and diagonalizes it directly on a
Monkhorst-Pack mesh, with tetrahedron, Blöchl-tetrahedron, or Gaussian
Brillouin-zone integration for the density of states. This route does not
go through the recursion machinery, and is normally faster for translationally
periodic bulk cells where diagonalization at each k-point is cheap relative
to the size of a real-space cluster that would represent the same physics.

Two further additions extend this route beyond band structure and DOS:

- A k-space Green's-function engine (``reciprocal_green.f90``) that fills
  the same Green's-function arrays the real-space recursion route fills,
  from either a strict Lehmann (eigenpair) representation or a direct
  Dyson-equation inversion. Selecting it (``gf_route = 'lehmann'`` or
  ``'dyson'``) lets exchange, damping, and conductivity post-processing
  run on k-space input instead of real-space recursion, without any
  change to those modules — the k-space and real-space routes become two
  interchangeable producers of the same downstream quantities.
- An exact Chebyshev moment generator (``reciprocal_moments.f90``) that
  computes conductivity moments directly from k-space eigenpairs, feeding
  the same moment arrays the recursion-based Chebyshev/KPM route produces,
  so the two can be compared directly on the same footing.

Both recursion and k-space routes support scalar-relativistic and fully
relativistic (spin-orbit coupled) treatments, and collinear or
non-collinear spin arrangements (``control%nsp``).

Geometries
==========

The code builds four kinds of real-space cluster, selected through
``control%calctype`` together with the corresponding pre-processing step:

- **Bulk** (``calctype = 'B'``, ``pre_processing = 'bravais'``) — an
  ordinary periodic crystal, built from its Bravais lattice and basis.
- **Surface** (``calctype = 'S'``, ``pre_processing = 'buildsurf'``) — a
  one-sided slab terminated by empty spheres representing vacuum, built
  from a converged bulk reference.
- **Impurity** (``calctype = 'I'``, ``pre_processing = 'newclubulk'`` /
  ``'newclusurf'``) — one or more defect atoms embedded in an otherwise
  bulk or surface host, with the host treated as a frozen reference medium
  around a self-consistent interaction zone.
- **Interface / layered** (``calctype = 'L'``, ``pre_processing =
  'buildinterface'``) — two frozen reference regions (A and B) joined at
  an interface, with a self-consistent buffer of active layers in
  between. Region B can be a second metallic reference (``&lattice
  region_b_kind = 'metal'``) or a semi-infinite vacuum whose frozen
  potential parameters are generated per run rather than hand-set
  (``region_b_kind = 'vacuum'``). A two-sided electrostatics and
  potential-alignment solve keeps the two regions' Fermi levels and
  potential zeros consistent across the interface.

All four geometries are built on the same embedding idea: a finite cluster
with a self-consistent interior and one or more frozen exterior regions
whose potential parameters are imported rather than converged. This is why
features like exchange coupling or Gilbert damping work identically
regardless of which geometry produced the underlying Green's function.

Substitutional disorder within a geometry (alloy-like concentration
mixing on a site or sublattice) is available through the Ruban-Abrikosov
concentration treatment (``control%ruban``, ``conca``/``concb``).

Magnetic and electronic structure
==================================

- **Collinear and non-collinear magnetism**, with or without spin-orbit
  coupling (four combinations via ``control%nsp``).
- **Spin spirals**: a generalized-Bloch-theorem construction lets a single
  chemical cell represent a spin spiral at an arbitrary (including
  incommensurate) wavevector ``q_ss``, rather than requiring an explicit
  supercell.
- **Frozen-magnon post-processing**: sweeps ``q_ss`` over a path and
  reports a magnon dispersion :math:`\omega(q)`, either from a single
  force-theorem SCF convergence (``mode = 'mft'``) or from independent SCF
  convergence at every q-point (``mode = 'scf'``); a multi-sublattice
  branch construction is available alongside the validated
  single-sublattice route.
- **Constraining fields**: site-dependent fields that steer the local
  magnetization toward a prescribed reference direction during SCF,
  needed to hold a non-collinear configuration (e.g. a spin spiral, or a
  metastable state) that would otherwise relax away under self-consistency.
- **LDA+U / DFT+U**: on-site Coulomb and exchange corrections for
  correlated (typically d or f) orbitals, with the double-counting
  potential built either from the standard Liechtenstein form or from the
  ACBN0 (self-consistently determined Hubbard parameters) construction.
- **Hyperfine coupling constants** at nuclear sites, computed after SCF
  convergence.
- **Surface and interface dipole electrostatics**: an optional l=1
  (dipole) correction to the otherwise monopole-only ASA Madelung
  potential, relevant to work functions and potential steps at broken-symmetry
  boundaries.

Measurable quantities (post-processing)
=========================================

Once a geometry has converged self-consistently (or, for the k-space
route, been diagonalized), the following properties can be computed:

- **Density of states**: total and site/orbital-resolved, via real-space
  recursion or k-space tetrahedron/Gaussian integration.
- **Band structure**: along a user-specified k-path (k-space route), or
  as a Bloch spectral function :math:`A(k,E)` built from either k-space
  Green's-function backend — the natural generalization of a band
  structure to a system with disorder- or correlation-broadened states.
- **Magnetic exchange interactions**: isotropic Heisenberg couplings
  :math:`J_{ij}`, Dzyaloshinskii-Moriya vectors :math:`\mathbf{D}_{ij}`,
  magnetic anisotropy tensors, and spin-lattice coupling parameters
  :math:`J_{ijk}`, from the real-space LKAG formalism; computable from
  real-space recursion or from either k-space Green's-function backend.
- **Gilbert damping**: the same torque-correlation machinery that
  produces exchange parameters also yields a damping tensor
  :math:`\alpha_{ij}`, on whichever route (recursion, k-space Lehmann, or
  k-space Dyson) produced the underlying Green's function.
- **Electrical conductivity**: the full :math:`3\times3` conductivity
  tensor via the Kubo-Bastin formalism and Chebyshev expansion, including
  DC and frequency-dependent (optical) conductivity, orbital-resolved
  contributions, and a stochastic-trace evaluation option for large
  systems. Moments can come from real-space Chebyshev recursion or from
  the exact k-space moment generator.
- **Orbital magnetic moments**, including a moment-of-inertia tensor
  formulation.
- **Atomistic spin dynamics**: real-time Landau-Lifshitz-Gilbert
  integration of the classical spin degrees of freedom, using the
  exchange parameters and (optionally) damping computed above as inputs.

Interoperability
=================

- **PAOFLOW**: Hamiltonians can be exported to the PAOFLOW pseudo-atomic
  orbital format for use in external post-processing, and (separately)
  imported from a PAOFLOW-derived tight-binding Hamiltonian to run the
  code's own exchange or conductivity post-processing on it.

Parallelization and hardware acceleration
============================================

- **OpenMP** shared-memory parallelization over the natural loop
  structures (atoms, Chebyshev sweeps, k-points).
- **MPI** distributed-memory parallelization, usable together with OpenMP.
- **Optional GPU acceleration** of the Chebyshev recursion kernels via a
  CUDA plugin, and optional Intel MKL-accelerated batched/sparse kernels
  for the same recursion route on CPU. Both are compile-time options and
  fall back cleanly to the portable CPU path when not enabled.

See Also
========

- :ref:`theory_index` — the formalism behind each item above
- :ref:`keywords/index` — the namelist keys that select each option
- :doc:`code_structure` — where each feature lives in the source tree
- :doc:`user_guide/examples` — worked examples for several of the geometries and properties above
