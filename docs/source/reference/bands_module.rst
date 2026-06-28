.. _reference/bands_module:

=======================================
Bands and Moments Module
=======================================

Overview
========

The **bands module** (``bands.f90``) calculates electronic structure properties from Green's functions in the LMTO-ASA framework:

- **Density of states (DOS)** - Total and projected
- **Spectral functions** - Energy-resolved electronic states
- **Band energy** - Electron kinetic energy contribution
- **Magnetic moments** - Spin and orbital moments
- **Magnetic torques** - Forces on magnetic moments
- **Fermi level** - Chemical potential determination
- **Occupation numbers** - Electron filling per orbital/atom

**Key features:**

- Real-space formulation (no explicit band structure, but DOS equivalent)
- Energy integration via Gauss-Legendre and Chebyshev methods
- Magnetic moment projections (spin + orbital angular momentum)
- Support for non-collinear magnetism
- MPI parallelization over atoms

Module Structure
================

Type: bands
-----------

.. code-block:: fortran

   type bands
      ! Pointers to system components
      class(green), pointer :: green
      class(lattice), pointer :: lattice
      class(energy), pointer :: en
      class(recursion), pointer :: recursion
      class(dos), pointer :: dos
      class(control), pointer :: control
      
      ! Electronic structure quantities
      real(rp) :: qqv                        ! Total valence electrons
      real(rp) :: eband                      ! Band energy
      real(rp) :: e1                         ! Energy integration cutoff
      integer :: nv1                         ! Energy mesh index
      
      ! Density of states arrays
      real(rp), allocatable :: dtot(:)       ! Total DOS
      real(rp), allocatable :: dx(:,:)       ! Projected DOS (x-component)
      real(rp), allocatable :: dy(:,:)       ! Projected DOS (y-component)
      real(rp), allocatable :: dz(:,:)       ! Projected DOS (z-component)
      real(rp), allocatable :: dup(:,:)      ! Spin-up DOS
      real(rp), allocatable :: ddw(:,:)      ! Spin-down DOS
      real(rp), allocatable :: d_orb(:,:,:,:,:) ! Orbital-resolved DOS
      
      ! Magnetic properties
      real(rp), allocatable :: mag_for(:,:)  ! Magnetic forces/torques
   end type

Theory: Density of States from Green's Functions
=================================================

Spectral Function and DOS
--------------------------

The **density of states (DOS)** is obtained from the imaginary part of the Green's function:

.. math::

   \rho(E) = -\frac{1}{\pi} \text{Im}\left[ \text{Tr} G^+(E) \right]

where $G^+(E) = (E - H + i\eta)^{-1}$ is the retarded Green's function.

For the **local DOS (LDOS)** at atom $i$:

.. math::

   \rho_i(E) = -\frac{1}{\pi} \sum_L \text{Im}\left[ G_{iL,iL}^+(E) \right]

where $L = (l, m, \sigma)$ combines orbital and spin indices.

Projected Density of States
----------------------------

**Orbital-resolved DOS:**

.. math::

   \rho_{i,l}(E) = -\frac{1}{\pi} \sum_{m=-l}^{l} \sum_{\sigma} 
   \text{Im}\left[ G_{i(lm\sigma), i(lm\sigma)}^+(E) \right]

giving contributions from $s, p, d, f$ orbitals.

**Spin-resolved DOS:**

.. math::

   \rho_i^{\uparrow}(E) &= -\frac{1}{\pi} \sum_{Lm} \text{Im}\left[ G_{i(Lm,\uparrow), i(Lm,\uparrow)}^+(E) \right] \\
   \rho_i^{\downarrow}(E) &= -\frac{1}{\pi} \sum_{Lm} \text{Im}\left[ G_{i(Lm,\downarrow), i(Lm,\downarrow)}^+(E) \right]

**Magnetic moment:**

.. math::

   m_i = \int_{-\infty}^{E_F} dE \left[ \rho_i^{\uparrow}(E) - \rho_i^{\downarrow}(E) \right]

Non-Collinear DOS Projections
------------------------------

For **non-collinear magnetism**, the DOS is projected along local magnetic moment direction $\hat{\mathbf{m}}_i = (m_x, m_y, m_z)$:

**Spin projection operator:**

.. math::

   P_i^{\pm} = \frac{1}{2}(1 \pm \hat{\mathbf{m}}_i \cdot \boldsymbol{\sigma})

where $\boldsymbol{\sigma} = (\sigma_x, \sigma_y, \sigma_z)$ are Pauli matrices.

**Projected DOS along $\hat{\mathbf{m}}_i$:**

.. math::

   \rho_i^{\parallel/\perp}(E) = -\frac{1}{\pi} \text{Im}\left[ 
   \text{Tr} P_i^{\pm} G_{ii}^+(E) \right]

**Implementation** (from ``calculate_moments``):

.. code-block:: text

   ! Project onto local moment direction
   do l = 1, 3  ! s, p, d orbitals
      do m = 1, 2*l - 1  ! magnetic quantum numbers
         o = (l-1)**2 + m  ! orbital index
         
         ! DOS with moment projection
         dspd(l, ie, na) = -Im[G_↑↑(o,o) + G_↓↓(o,o)]  ! Total
                          - m_z * Im[G_↑↑(o,o) - G_↓↓(o,o)]  ! z-component
                          - m_y * Im[i(G_↑↓(o,o) - G_↓↑(o,o))]  ! y-component
                          - m_x * Re[G_↑↓(o,o) + G_↓↑(o,o)]  ! x-component
      end do
   end do

Band Energy and Electron Count
===============================

Band Energy Functional
-----------------------

The **band energy** (kinetic + potential energy of electrons):

.. math::

   E_{\text{band}} = \int_{-\infty}^{E_F} dE \, E \, \rho_{\text{tot}}(E)

This is the **energy integral of DOS**:

.. code-block:: fortran

   call simpson_m(eband, en%edel, en%fermi, nv1, dtot, e1, 1, en%ene)

where ``simpson_m`` performs energy-weighted Simpson integration.

**Physical meaning:** Electronic contribution to total energy (excluding double-counting corrections).

Total Valence Electrons
------------------------

The **number of valence electrons** $Q_v$:

.. math::

   Q_v = \int_{-\infty}^{E_F} dE \, \rho_{\text{tot}}(E)

**Fermi level** $E_F$ is determined self-consistently by:

.. math::

   \int_{-\infty}^{E_F} dE \, \rho(E) = N_{\text{electrons}}

where $N_{\text{electrons}}$ is the target electron count.

Fermi Level Determination
--------------------------

**Procedure** ``fermi()``:

1. **Simpson integration** of DOS from $E_{\min}$ upward:

.. math::

   Q(E) = \int_{E_{\min}}^{E} dE' \, \rho(E')

2. **Find $E_F$** where $Q(E_F) = Q_v$:

   - If exact match found: $E_F = E_i$
   - Otherwise: Linear interpolation between $E_i$ and $E_{i+1}$

.. math::

   E_F = E_i + \frac{Q_v - Q(E_i)}{\rho(E_i)} \quad \text{(linear approx)}

3. **Set integration cutoff** ``e1`` and mesh index ``nv1``

Magnetic Moments
================

Spin Magnetic Moment
--------------------

The **spin magnetic moment** $\mu_s^i$ at atom $i$:

.. math::

   \mu_s^i = \int_{-\infty}^{E_F} dE \left[ \rho_i^{\uparrow}(E) - \rho_i^{\downarrow}(E) \right]

In **non-collinear case**, project onto local quantization axis:

.. math::

   \mu_s^{i,\alpha} = \int_{-\infty}^{E_F} dE \, 
   m_i^{\alpha} \left[ \rho_i^{\uparrow}(E) - \rho_i^{\downarrow}(E) \right]

for $\alpha \in \{x, y, z\}$.

**Implementation** (``calculate_magnetic_moments``):

.. code-block:: text

   do na = 1, lattice%nrec
      do l = 0, lmax  ! Sum over orbital angular momenta
         ! Spin moment from charge difference
         spin_moment(na) += ql(l, spin_up) - ql(l, spin_down)
      end do
      
      ! Project onto moment direction
      moment_vector(:, na) = spin_moment(na) * mom(:) / |mom|
   end do

Orbital Magnetic Moment
------------------------

The **orbital magnetic moment** $\mu_L^i$ arises from spin-orbit coupling:

.. math::

   \mu_L^i = -\mu_B \langle \mathbf{L} \rangle_i

where $\mathbf{L}$ is the orbital angular momentum operator and $\mu_B$ is the Bohr magneton.

**In LMTO**, calculated from off-diagonal Green's function elements:

.. math::

   \langle L_z \rangle = \sum_{lm} m \int_{-\infty}^{E_F} dE \, \rho_{lm}(E)

where $m$ is the magnetic quantum number.

**General direction** (non-collinear):

.. math::

   \langle \mathbf{L} \rangle = \int_{-\infty}^{E_F} dE \, 
   \text{Tr}\left[ \hat{\mathbf{L}} \cdot \text{Im}[G^+(E)] \right]

**Implementation** (``calculate_orbital_moments``):

.. code-block:: fortran

   do l = 1, 3  ! p, d, f orbitals
      do m = -l, l
         ! Orbital moment from L operator expectation value
         orbital_moment(na) += m * charge_lm(l, m, na)
      end do
   end do

Magnetic Torques
----------------

The **magnetic torque** $\mathbf{T}_i = \mathbf{m}_i \times \mathbf{B}_i^{\text{eff}}$ represents the force trying to rotate the magnetic moment.

**Effective field** from energy gradient:

.. math::

   \mathbf{B}_i^{\text{eff}} = -\frac{1}{\mu_i} \frac{\partial E}{\partial \hat{\mathbf{m}}_i}

**Torque components** (``calculate_magnetic_torques``):

.. math::

   T_i^{\alpha} = \sum_{\beta\gamma} \epsilon_{\alpha\beta\gamma} m_i^{\beta} B_i^{\gamma}

where $\epsilon$ is the Levi-Civita tensor.

**Used for:**

- Non-collinear spin dynamics
- Constrained DFT calculations
- Magnetic anisotropy energy (MAE) via torque method

Integration Methods
===================

Simpson's Rule
--------------

**Standard Simpson integration** for smooth integrands:

.. math::

   \int_a^b f(x) dx \approx \frac{h}{3} \sum_{i=0}^{N/2-1} 
   \left[ f(x_{2i}) + 4f(x_{2i+1}) + f(x_{2i+2}) \right]

where $h = (b-a)/N$ and $N$ is even.

**Implementation:** ``simpson_f(result, x, x_upper, npts, y, ...)``

**Accuracy:** $\mathcal{O}(h^4)$ for smooth functions.

Gauss-Legendre Quadrature
--------------------------

For **higher accuracy** with fewer points, use **Gauss-Legendre**:

.. math::

   \int_{-1}^{1} f(x) dx \approx \sum_{i=1}^{N} w_i f(x_i)

where $x_i$ are roots of Legendre polynomial $P_N(x)$ and $w_i$ are weights.

**Advantages:**

- Exact for polynomials up to degree $2N-1$
- Fewer points needed for given accuracy
- Optimal for smooth integrands

**Implementation:** ``calculate_moments_gauss_legendre()``

**Transform** to arbitrary interval $[a, b]$:

.. math::

   \int_a^b f(x) dx = \frac{b-a}{2} \int_{-1}^{1} f\left(\frac{b-a}{2}t + \frac{a+b}{2}\right) dt

Chebyshev-Gauss Quadrature
---------------------------

For **oscillatory** or **nearly singular** functions, use **Chebyshev-Gauss**:

.. math::

   \int_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}} dx \approx \sum_{i=1}^{N} w_i f(x_i)

where $x_i = \cos\frac{(2i-1)\pi}{2N}$ and $w_i = \frac{\pi}{N}$.

**Implementation:** ``calculate_moments_chebgauss()``

**Advantage:** Natural for Chebyshev expansion (Kernel Polynomial Method).

Main Procedures
===============

calculate_projected_dos
-----------------------

**Purpose:** Calculate orbital and site-projected DOS.

**Algorithm:**

1. Loop over all atoms and energy points
2. For each atom $i$ and energy $E$:

.. code-block:: text

   do na = 1, lattice%nrec
      do ie = 1, en%channels_ldos
         ! Total DOS at this atom
         dtot(ie) += -Im[Tr(G_ii(E))] / pi
         
         ! Projected onto orbitals
         do l = 0, lmax
            d_orb(l, ie, na) = -Im[sum_m G_ii(lm, lm, E)] / pi
         end do
         
         ! Spin-resolved
         dup(ie, na) = -Im[Tr(G_ii^↑↑(E))] / pi
         ddw(ie, na) = -Im[Tr(G_ii^↓↓(E))] / pi
      end do
   end do

3. **Output** to DOS files

calculate_moments
-----------------

**Purpose:** Calculate electronic moments $m^{(q)}$ with $q = 0, 1, 2$.

**Definition:**

.. math::

   m^{(q)} = \int_{-\infty}^{E_F} dE \, E^q \rho(E)

- $m^{(0)} = Q$ - Total electron count
- $m^{(1)} = E_{\text{band}}$ - Band energy
- $m^{(2)}$ - Second moment (related to kinetic energy)

**Implementation:**

.. code-block:: fortran

   do q = 0, 2
      moment(q) = 0.0
      do ie = 1, nv1
         moment(q) += dtot(ie) * (en%ene(ie))**q * weight(ie)
      end do
   end do

calculate_magnetic_moments
---------------------------

**Purpose:** Calculate spin and orbital magnetic moments for each atom.

**Algorithm:**

1. **Charge integration** per orbital:

.. math::

   Q_{l,\sigma}^i = \int_{-\infty}^{E_F} dE \, \rho_{i,l,\sigma}(E)

2. **Spin moment:**

.. math::

   \mu_s^i = \sum_l \left( Q_{l,\uparrow}^i - Q_{l,\downarrow}^i \right)

3. **Orbital moment:** From $L$-operator expectation value

4. **Total moment:**

.. math::

   \mathbf{\mu}_i = \mathbf{\mu}_s^i + \mathbf{\mu}_L^i

5. **Output** to ``moments.out``

calculate_band_energy
---------------------

**Purpose:** Calculate total band energy contribution to DFT total energy.

**Formula:**

.. math::

   E_{\text{band}} = \int_{-\infty}^{E_F} dE \, E \, \rho_{\text{tot}}(E)

**Implementation:**

.. code-block:: fortran

   call simpson_m(eband, en%edel, en%fermi, nv1, dtot, e1, 1, en%ene)

where ``simpson_m`` is energy-weighted Simpson integration.

**Usage:** Part of total energy in SCF:

.. math::

   E_{\text{tot}} = E_{\text{band}} + E_{\text{double-count}} + E_{\text{xc}} + E_{\text{Madelung}} + \ldots

calculate_fermi_gauss
---------------------

**Purpose:** Determine Fermi level using Gauss-Legendre integration.

**Algorithm:**

1. Transform DOS to Gauss-Legendre points $\{x_i, w_i\}$
2. Integrate charge:

.. math::

   Q(E) = \sum_{i: E_i < E} w_i \rho(E_i)

3. Find $E_F$ where $Q(E_F) = Q_v$ by interpolation
4. More accurate than Simpson for given number of points

Output Files
============

dos_total.dat
-------------

Total DOS as function of energy:

.. code-block:: text

   # E(eV)   DOS(states/eV)
     -10.0   5.234
     -9.98   5.456
     ...
      0.00   8.123  ← Fermi level
     ...

dos_atom_<i>.dat
----------------

DOS projected onto atom $i$:

.. code-block:: text

   # E(eV)   DOS_tot   DOS_s   DOS_p   DOS_d   DOS_up   DOS_down
     -10.0   1.234     0.123   0.456   0.655   0.678    0.556
     ...

**Columns:**

- Total LDOS at atom $i$
- $s$, $p$, $d$ orbital projections
- Spin-up and spin-down contributions

moments.out
-----------

Magnetic moments per atom:

.. code-block:: text

   # Atom   μ_spin(μB)   μ_orb(μB)   μ_tot(μB)   m_x   m_y   m_z
      1      2.234        0.123       2.357      0.0   0.0   1.0
      2     -2.156       -0.098      -2.254      0.0   0.0  -1.0
      ...

**Columns:**

- Atom index
- Spin magnetic moment
- Orbital magnetic moment
- Total magnetic moment
- Direction cosines $(m_x, m_y, m_z)$

band_energy.out
---------------

Band energy contribution:

.. code-block:: text

   Band energy: -245.6723 eV

Integration and Usage
=====================

Workflow
--------

1. **SCF convergence** must complete first
2. **DOS calculation** automatic during SCF or post-processing
3. **Access results:**

   - Read DOS files for plotting
   - Extract Fermi level from output
   - Use magnetic moments for spin dynamics

**Typical tasks:**

- Plot density of states
- Analyze electronic structure
- Extract magnetic properties
- Calculate derived quantities (effective mass, etc.)

Connection to Other Modules
----------------------------

**Green's function** $\to$ **Bands module** $\to$ **Properties**:

- ``green.f90`` provides $G^+(E)$ at each energy
- ``bands.f90`` extracts DOS: $\rho(E) = -\text{Im}[G^+]/\pi$
- Properties calculated from DOS integrals

**Self-consistency:**

.. math::

   \rho(E) \to n(\mathbf{r}) \to V(\mathbf{r}) \to H \to G(E) \to \rho(E)

**Exchange calculation:**

- Magnetic moments from ``bands`` used in ``exchange`` module
- Spin directions for non-collinear exchange

Performance Considerations
==========================

Computational Cost
------------------

**DOS calculation:**

- $\mathcal{O}(N_{\text{atoms}} \times N_E \times N_L^2)$
- $N_E \sim 200-500$: Energy mesh points
- $N_L \sim 9-16$: Orbital basis size
- **Typical:** Few seconds to minutes

**Moment integration:**

- Dominated by SCF, not post-processing
- Simpson vs Gauss-Legendre: Similar cost, GL more accurate

Accuracy Tips
-------------

**For converged results:**

1. **Dense energy mesh:** ``channels_ldos`` $\geq 300$
2. **Proper Fermi smearing:** ``temp`` $\sim 0.001-0.01$ eV (if using)
3. **Converged SCF:** $\Delta Q < 10^{-5}$
4. **Sufficient recursion:** ``ll_gf`` $\geq 50$

Provenance
==========

**Implementation:**

- ``source/bands.f90`` - Main bands module (1193 lines)
- ``source/density_of_states.f90`` - DOS utilities
- ``source/green.f90`` - Green's functions

**Key procedures:**

- ``calculate_projected_dos()`` - Orbital/site-projected DOS
- ``calculate_moments()`` - Electronic moments
- ``calculate_magnetic_moments()`` - Spin/orbital moments
- ``calculate_band_energy()`` - Band energy integral
- ``fermi()`` - Fermi level determination

**Authors:**

- Angela Klautau, Ramon Cardias, Lucas P. Campagna
- S. Frota-Pessôa, Pascoal R. Peduto
- Anders Bergman, Ivan P. Miranda
- S. B. Legoas, H. M. Petrilli

References
==========

1. **DOS from Green's functions:** A. P. Sutton and R. W. Balluffi,
   *Interfaces in Crystalline Materials*, Oxford University Press (1995)

2. **Magnetic moments in DFT:** R. M. Martin,
   *Electronic Structure: Basic Theory and Practical Methods*, Cambridge (2004)

3. **Numerical integration:** W. H. Press et al.,
   *Numerical Recipes*, Cambridge University Press (2007)

4. **Orbital magnetism in LMTO:** O. Eriksson et al.,
   *Phys. Rev. B* **41**, 7311 (1990)

5. **Magnetic torques:** P. H. Dederichs et al.,
   *Phys. Rev. Lett.* **53**, 2512 (1984)

See Also
========

- :ref:`theory/green_functions` - Green's function theory
- :ref:`theory/scf_cycle` - Self-consistent field procedure
- :ref:`reference/exchange_module` - Exchange interactions using moments
- :ref:`user_guide/output_files` - DOS and moment output files
