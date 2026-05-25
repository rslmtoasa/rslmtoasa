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
   lattice_geometry
   basis_parameters
   energy_mesh
   scf_settings
   exchange_correlation
   output_options
   constraints
   spin_dynamics
   element_parameters

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

**&constraints**

Constraining fields for non-collinear magnetism. See :ref:`keywords/constraints`.

**&sd**

Atomistic spin dynamics time-integration parameters. See :ref:`keywords/spin_dynamics`.

**&element**

Element database (atomic number, core/valence split, quantum numbers). See :ref:`keywords/element_parameters`.

**&conductivity**

Chebyshev expansion depth for transport calculations. See :ref:`reference/conductivity_module`.

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
   * - max_iterations
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
   * - dq_tol
     - real
     - 1e-5 to 1e-4
   * - mixing
     - str
     - 'linear' or 'broyden'
   * - alpha
     - real
     - 0.3-0.5
   * - channels_ldos
     - int
     - 300-500

Alphabetical Listing (All Parameters)
=====================================

**A**

- ``alpha`` - Mixing parameter (density mixing)
- ``alat`` - Lattice constant

**B**

- ``broyden_history`` - Number of previous densities to keep (Broyden mixing)

**A**

- ``alpha`` (sd) - Gilbert damping parameter (ASD)
- ``asd_jij`` - Include J_ij in ASD effective field
- ``atomic_number`` - Element atomic number Z

**B**

- ``blockrec`` - Block recursion flag
- ``broyden_history`` - Number of previous densities to keep (Broyden mixing)
- ``b_stochastic`` - Stochastic field array (ASD)

**C**

- ``calctype`` - Top-level calculation type
- ``center_band`` - LMTO orbital center energy
- ``channels_ldos`` - Number of energy mesh points
- ``cold`` - Cold start (ASA potential extraction)
- ``cond_calctype`` - Conductivity sub-type
- ``cond_ll`` - Conductivity Chebyshev cutoff (in &control)
- ``cond_type`` - Conductivity formulation
- ``conca``, ``concb`` - Alloy concentrations
- ``constraints_bfield`` - Initial constraining field values per atom
- ``constraints_code_prefac`` - Scaling prefactor for constraining field
- ``constraints_enable`` - Enable constraining field functionality
- ``constraints_i_cons`` - Constraint algorithm selector (2–5)
- ``constraints_mom_ref`` - Reference magnetization directions per atom
- ``conv_thr`` - SCF convergence threshold (self module)
- ``core`` - Number of core electrons

**D**

- ``do_asd`` - Enable atomistic spin dynamics
- ``do_cochg`` - Enable co-charge correction
- ``do_comom`` - Enable co-moment correction
- ``dq_tol`` - Charge density convergence tolerance
- ``dt`` - ASD time step

**E**

- ``energy_max`` - Upper energy limit
- ``energy_min`` - Lower energy limit

**F**

- ``f_core`` - Number of f-core electrons
- ``fermi`` - Fermi energy (initial guess or fixed)
- ``fix_fermi`` - Fix Fermi level to value
- ``fix_soc`` - Fix spin-orbit coupling magnitude
- ``freeze`` - Freeze SCF potential
- ``hyperfine`` - Compute hyperfine coupling constants

**I**

- ``init`` - SCF initialisation mode
- ``integrator`` - ASD integration scheme (depon/heun)

**L**

- ``ll_cond`` - Conductivity Chebyshev cutoff (in &conductivity)
- ``lld`` - Recursion cutoff for d electrons
- ``llsp`` - Recursion cutoff for s/p electrons
- ``lmax`` - Maximum orbital angular momentum

**M**

- ``magnetic_mixing`` - Separate magnetic mixing in SCF
- ``max_iterations`` - Maximum SCF iterations
- ``mix_all`` - Mix all potential components together
- ``mixmag`` - Per-atom magnetic mixing parameter
- ``mixmag_all`` - Apply same magnetic mixing to all atoms
- ``mixing`` - Density mixing type

**N**

- ``nbulk`` - Number of bulk atoms in cluster
- ``nlim`` - Cluster size control
- ``nsp`` - Type of calculation (relativistic, collinear/non-collinear)
- ``nstep`` - Number of SCF steps
- ``nt`` - Number of ASD time steps
- ``num_quant_d`` - Principal quantum number for d channel
- ``num_quant_p`` - Principal quantum number for p channel
- ``num_quant_s`` - Principal quantum number for s channel
- ``nx, ny, nz`` - Cluster dimensions in lattice coordinates

**O**

- ``idos`` - LDOS output type
- ``orbital_polarization`` - Enable orbital polarization correction
- ``orb_pol`` - Include orbital polarization

**R**

- ``random_vec_num`` - Number of stochastic vectors for moment calculation
- ``r2`` - Cluster cutoff radius squared
- ``recur`` - Recursion algorithm (lanczos/chebyshev)
- ``rigid_band`` - Rigid band shift
- ``ruban`` - Ruban-Abrikosov alloy concentration parameter

**S**

- ``sd_temp`` - ASD simulation temperature
- ``soc_scale`` - Spin-orbit coupling scaling factor
- ``svac`` - Spherically averaged atomic spheres
- ``symbol`` - Element chemical symbol
- ``sws`` - Wigner-Seitz radius (atomic sphere)
- ``sym_term`` - Symmetry terminator correction

**T**

- ``t_f`` - ASD final simulation time
- ``t_i`` - ASD initial simulation time
- ``temperature`` - Electronic temperature (Kelvin)
- ``terminator`` - Recursion terminator type
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

- :ref:`keywords/control_parameters` - Calculation control (nsp, max_iterations, etc.)
- :ref:`keywords/lattice_geometry` - Structure (alat, lattice type, cluster size)
- :ref:`keywords/basis_parameters` - LMTO basis (lmax, center_band, width_band)
- :ref:`keywords/energy_mesh` - Energy integration (channels_ldos, energy range)
- :ref:`keywords/scf_settings` - SCF convergence (mixing, dq_tol, max_iterations)
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
