.. _keywords/index:

=======================================
Input Keywords Reference
=======================================

Complete alphabetical and categorical listing of all RS-LMTO-ASA input parameters.

Quick Navigation
================

.. toctree::
   :maxdepth: 2

   control_parameters
   hamiltonian_parameters
   lattice_geometry
   basis_parameters
   energy_mesh
   scf_settings
   exchange_correlation
   output_options
   constraints
   spin_dynamics
   element_parameters
   frozen_magnon
   reciprocal
   spin_response
   lda_u

Parameters by Namelist
======================

**&control**

Main calculation control parameters. See :ref:`keywords/control_parameters`.

**&lattice**

Crystal structure and geometry. See :ref:`keywords/lattice_geometry`.

**&energy**

Energy mesh and Fermi level settings. See :ref:`keywords/energy_mesh`.

**&self**

Self-consistent field convergence. See :ref:`keywords/scf_settings`.

**&element**

Atomic parameters (element database).

**&par**

LMTO tight-binding potential parameters. See :ref:`keywords/basis_parameters`.

**&reciprocal**, **&kpath**

k-space SCF, band structure, DOS, Fermi surface and k-space Green functions. See :ref:`keywords/reciprocal`.

**&spin_response**

Transverse spin susceptibility on the k-space route. See :ref:`keywords/spin_response`.

**LDA+U** (``&hamiltonian`` Hubbard keys, ``&mix ldm_beta``, ``&par`` Hubbard data, ``&density_matrix``)

See :ref:`keywords/lda_u`.

**&constraints**

Constraining fields for non-collinear magnetism. See :ref:`keywords/constraints`.

**&sd**

Atomistic spin dynamics time-integration parameters. See :ref:`keywords/spin_dynamics`.

**&element**

Element database (atomic number, core/valence split, quantum numbers). See :ref:`keywords/element_parameters`.

**&conductivity**

Chebyshev expansion depth for transport calculations. See :ref:`reference/conductivity_module`.

**&frozen_magnon**

Frozen-magnon E(q)/omega(q) sweep configuration (B1). See :ref:`keywords/frozen_magnon`.

Quick Reference Table
=====================

Most commonly used parameters:

.. list-table::
   :widths: 20 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Typical Value
   * - nsp
     - int
     - 1 (scalar-rel), 2 (collinear)
   * - nstep
     - int
     - 50-100
   * - llsp
     - int
     - 60-100
   * - lld
     - int
     - 60-100
   * - alat
     - real
     - system-dependent (Å)
   * - conv_thr
     - real
     - 1e-9 to 1e-5
   * - mixtype
     - str
     - 'linear' or 'broyden'
   * - beta
     - real
     - 0.05-0.5
   * - channels_ldos
     - int
     - 300-500

Alphabetical Listing (All Parameters)
=====================================

**A**

- ``alat`` - Lattice constant

**A**

- ``alpha`` (sd) - Gilbert damping parameter (ASD)
- ``asd_jij`` - Include J_ij in ASD effective field
- ``atomic_number`` - Element atomic number Z

**B**

- ``beta`` - Density mixing parameter (``&mix``)
- ``blockrec`` - Block recursion flag
- ``b_stochastic`` - Stochastic field array (ASD)

**C**

- ``calctype`` - Cluster/embedding topology: 'B'/'S'/'I'/'L' (bulk/surface/impurity/interface) —
  **not** the post-SCF property, see ``post_processing``
- ``center_band`` - LMTO orbital center energy
- ``channels_ldos`` - Number of energy mesh points
- ``cold`` - Cold start (ASA potential extraction)
- ``cond_calctype`` - Conductivity sub-type
- ``cond_ll`` - Conductivity Chebyshev cutoff (in &control)
- ``conca``, ``concb`` - Alloy concentrations
- ``constraints_bfield`` - Initial constraining field values per atom
- ``constraints_code_prefac`` - Scaling prefactor for constraining field
- ``constraints_diagnostics`` - Enable detailed constraint diagnostics
- ``constraints_enable`` - Enable constraining field functionality
- ``constraints_i_cons`` - Constraint algorithm selector (2–5)
- ``constraints_mom_ref`` - Reference magnetization directions per atom
- ``constraints_tolerance`` - RMS angular constraint convergence tolerance
- ``conv_thr`` - SCF convergence threshold (self module)
- ``core`` - Number of core electrons

**D**

- ``dipole_electrostatics`` - Enable l=1 dipole surface/interface electrostatics (B6)
- ``dipole_mix`` - SCF mixing fraction for the dipole moment Q10 (B6)
- ``density_policy`` - Rotating-frame SCF density policy (WP7)
- ``do_asd`` - Enable atomistic spin dynamics
- ``do_cochg`` - Enable co-charge correction
- ``do_comom`` - Enable co-moment correction
- ``exchange_delta`` - Exchange splitting source for imported Hamiltonians
- ``do_chi0`` - Static transverse chi0_ij in exchange post-processing
- ``do_damping`` - Enable Gilbert-damping evaluation in exchange post-processing (B5)
- ``do_sigma`` - One-shot G_ii(z) for a site-local self-energy read from sigma.dat (gloc_sigma.dat)
- ``dt`` - ASD time step

**E**

- ``energy_max`` - Upper energy limit
- ``energy_min`` - Lower energy limit

**F**

- ``f_core`` - Number of f-core electrons
- ``fermi`` - Fermi energy (initial guess or fixed)
- ``fermi_a``, ``fermi_b`` - Per-region source Fermi level for interface alignment check (B7)
- ``fix_fermi`` - Fix Fermi level to value
- ``fix_fermi_to_region`` - Pin Fermi level to a named frozen region (B7, calctype='L')
- ``fix_soc`` - Fix spin-orbit coupling magnitude
- ``freeze`` - Freeze SCF potential

**G**

- ``gf_route`` - Route selector for Green's-function producer: recursion/lehmann/dyson (B2/B5)
- ``green_backend`` - k-space Green's-function backend: lehmann/dyson (B2)
- ``green_eta`` - Retarded broadening for the k-space Green's function (B2)

**H**

- ``hyperfine`` - Compute hyperfine coupling constants

**I**

- ``init`` - SCF initialisation mode
- ``integrator`` - ASD integration scheme (depon/heun)

**L**

- ``linear_in``, ``linear_out`` - Linear-response input and output operators (in &control)
- ``ll_cond`` - Read but not used; the conductivity cutoff is ``cond_ll`` (in &control)
- ``lld`` - Recursion cutoff (number of recursion levels)
- ``llsp`` - Enters only the array size ``max(llsp, lld)``
- ``lmax`` - Maximum orbital angular momentum

**M**

- ``magnetic_mixing`` - Separate magnetic mixing in SCF
- ``mix_all`` - Mix all potential components together
- ``mixmag`` - Per-atom magnetic mixing parameter
- ``mixtype`` - Density mixing type, 'linear' or 'broyden' (``&mix``)
- ``mixmag_all`` - Apply same magnetic mixing to all atoms

**N**

- ``nbulk`` - Read but not used; the code sets it internally
- ``nlay_a``, ``nlay_b`` - Layer counts for interface regions A/B (B7, calctype='L'; also a
  separate, differently-scoped &charge pair of the same name — see :ref:`keywords/lattice_geometry`)
- ``nlim`` - Cluster size control
- ``nsp`` - Type of calculation (relativistic, collinear/non-collinear)
- ``nstep`` - Number of SCF steps
- ``nt`` - Number of ASD time steps
- ``num_quant_d`` - Principal quantum number for d channel
- ``num_quant_p`` - Principal quantum number for p channel
- ``num_quant_s`` - Principal quantum number for s channel

**O**

- ``idos`` - LDOS output type
- ``orbital_polarization`` - Enable orbital polarization correction
- ``orb_pol`` - Include orbital polarization

**Q**

- ``q_ss`` - Spin-spiral wavevector (Cartesian units of 2*pi/alat)
- ``q10`` - Computed l=1 dipole moment per atom, output only (B6)

**R**

- ``random_vec_num`` - Number of stochastic vectors for moment calculation
- ``r2`` - Cluster cutoff radius squared
- ``recur`` - Recursion algorithm (lanczos/chebyshev/block)
- ``region_b_kind`` - Physical kind of interface region B: metal/vacuum (B7, calctype='L')
- ``rigid_band`` - Rigid band shift
- ``ruban`` - Ruban-Abrikosov alloy concentration parameter

**S**

- ``sd_temp`` - ASD simulation temperature
- ``sigma_tol`` - Tail bound that marks a do_sigma G_ii(z) block converged
- ``soc_scale`` - Spin-orbit coupling scaling factor
- ``svac`` - Spherically averaged atomic spheres
- ``symbol`` - Element chemical symbol
- ``sws`` - Wigner-Seitz radius (atomic sphere)
- ``sym_term`` - Symmetry terminator correction

**T**

- ``t_f`` - ASD final simulation time
- ``t_i`` - ASD initial simulation time
- ``temperature`` - Temperature (Kelvin) of the exchange-interaction energy integration (``&energy``)
- ``terminator`` - Recursion terminator type
- ``theta_ss`` - Spin-spiral cone angle (degrees)
- ``txc`` - Exchange-correlation functional type

**V**

- ``valence`` - Number of valence electrons
- ``verbose`` - Enable verbose output

**W**

- ``width_band`` - LMTO bandwidth parameter
- ``ws`` - Per-atom Wigner-Seitz radii
- ``ws_all`` - Apply uniform WS radius to all atoms
- ``ws_r`` - Wigner-Seitz radius (input format)

See Full Reference Pages
========================

For detailed descriptions, constraints, and examples, see:

- :ref:`keywords/control_parameters` - Calculation control (nsp, lld, etc.)
- :ref:`keywords/lattice_geometry` - Structure (alat, lattice type, cluster size)
- :ref:`keywords/basis_parameters` - LMTO basis (lmax, center_band, width_band)
- :ref:`keywords/energy_mesh` - Energy integration (channels_ldos, energy range)
- :ref:`keywords/scf_settings` - SCF convergence (nstep, conv_thr, mixtype, beta)
- :ref:`keywords/exchange_correlation` - XC functional (txc)
- :ref:`keywords/output_options` - Output control (post_processing, idos)
- :ref:`keywords/constraints` - Constraining fields (constraints_enable, constraints_i_cons, etc.)
- :ref:`keywords/spin_dynamics` - Spin dynamics (t_i, t_f, dt, alpha, integrator, etc.)
- :ref:`keywords/element_parameters` - Element database (symbol, atomic_number, core, valence, etc.)

Cross-References
================

Each parameter page includes:

- **Purpose** - Physical meaning and use
- **Allowed values** - Valid ranges, defaults
- **Units** - What units are expected
- **Example** - Typical usage
- **Related parameters** - Cross-links to related options
- **Code location** - Where in source code it's used

Provenance
==========

Parameter definitions and defaults come from:

- ``source/control.f90`` - Control parameters type and defaults
- ``source/lattice.f90`` - Lattice parameters
- ``source/energy.f90`` - Energy mesh parameters
- ``source/self.f90`` - SCF parameters
- ``source/namelist_generator.f90`` - Namelist generation utilities
- Example files: ``example/*/*.nml`` - Typical values

See Also
========

- :ref:`user_guide/input_files` - Input file syntax and format
- :doc:`../user_guide/examples` - Worked examples with parameter settings
- :doc:`../theory/scf_cycle` - Understanding SCF-related parameters
- :doc:`../theory/recursion_method` - Recursion cutoff parameters
