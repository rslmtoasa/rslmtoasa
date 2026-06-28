.. _reference/exchange_module:

=======================================
Exchange Interactions Module
=======================================

Overview
========

The **exchange module** (``exchange.f90``) calculates magnetic exchange interactions using the **Liechtenstein-Katsnelson-Antropov-Guslienko (LKAG) formalism** within the LMTO-ASA framework. This module computes:

- Isotropic Heisenberg exchange interactions $J_{ij}$
- Dzyaloshinskii-Moriya interactions (DMI) $\mathbf{D}_{ij}$
- Magnetic anisotropy tensors $\mathbf{A}_{ij}$
- Gilbert damping parameters $\alpha_{ij}$
- Spin-lattice coupling parameters $J_{ijk}$
- Moment of inertia tensors

**Key features:**

- Real-space formulation (no k-points)
- Green's function based (energy-resolved)
- Supports collinear and non-collinear magnetism
- Includes spin-orbit coupling effects
- MPI parallelization for large systems

Module Structure
================

Type: exchange
--------------

The ``exchange`` type stores exchange coupling data and provides calculation methods:

.. code-block:: fortran

   type exchange
      ! Pointers to other modules
      class(green), pointer :: green
      class(lattice), pointer :: lattice
      class(hamiltonian), pointer :: hamiltonian
      class(control), pointer :: control
      
      ! Exchange interaction parameters
      real(rp) :: jij, jijcd, jijsd, jijcc, jijsc  ! Heisenberg J
      real(rp), dimension(9) :: jij_aux             ! J tensor (auxiliary GF)
      real(rp) :: jij00_aux                         ! On-site J0
      real(rp), dimension(9) :: jijk                ! Spin-lattice coupling
      real(rp), dimension(3) :: dmi, dmisc, dmicc   ! DMI vectors
      real(rp), dimension(3,3) :: aij, aijsd, aijsc ! Anisotropy tensors
   end type

Theory: LKAG Formalism
======================

Liechtenstein-Katsnelson-Antropov-Guslienko Formula
----------------------------------------------------

The exchange interactions are calculated using the LKAG approach, which expresses the exchange parameters as energy integrals over Green's functions.

**Isotropic Heisenberg exchange** $J_{ij}$:

.. math::

   J_{ij} = -\frac{1}{4\pi} \text{Im} \int_{-\infty}^{E_F} dE \, 
   \text{Tr}\left[ \Delta_i G_{ij}^{\uparrow\uparrow}(E) 
   \Delta_j G_{ji}^{\downarrow\downarrow}(E) \right]

where:

- $\Delta_i = P_i^{\uparrow} - P_i^{\downarrow}$ is the spin-split potential function
- $P_i^{\sigma}$ are the potential functions for spin $\sigma$
- $G_{ij}^{\sigma\sigma'}(E)$ is the intersite Green's function from site $j$ to site $i$
- $E_F$ is the Fermi energy
- Tr denotes the trace over orbital indices

**Alternative tensor formulation:**

For non-collinear systems, the exchange can be expressed as a tensor:

.. math::

   J_{ij}^{\alpha\beta} = -\frac{1}{2\pi} \text{Im} \int_{-\infty}^{E_F} dE \, 
   \text{Tr}\left[ \Delta_i G_{ij}^{\alpha}(E) \Delta_j G_{ji}^{\beta}(E) \right]

where $\alpha, \beta \in \{x, y, z\}$ label spatial directions and:

.. math::

   G_{ij}^{\alpha}(E) = \frac{1}{2}\left( G_{ij}^{\uparrow\uparrow} + G_{ij}^{\downarrow\downarrow} \right) 
   + \frac{\mathbf{m}_i \cdot \mathbf{e}_\alpha}{2} 
   \left( G_{ij}^{\uparrow\uparrow} - G_{ij}^{\downarrow\downarrow} \right)

The isotropic exchange is then $J_{ij} = \frac{1}{3}(J_{ij}^{xx} + J_{ij}^{yy} + J_{ij}^{zz})$.

Dzyaloshinskii-Moriya Interaction (DMI)
----------------------------------------

The **DMI vector** $\mathbf{D}_{ij}$ arises from spin-orbit coupling and broken inversion symmetry:

.. math::

   \mathbf{D}_{ij} = -\frac{1}{2\pi} \text{Im} \int_{-\infty}^{E_F} dE \, 
   \text{Tr}\left[ \Delta_i G_{ij}^{x}(E) \Delta_j G_{ji}^{y}(E) - 
   \Delta_i G_{ij}^{y}(E) \Delta_j G_{ji}^{x}(E), \ldots \right]

The three components are:

.. math::

   D_{ij}^z &= -\frac{1}{2\pi} \text{Im} \int dE \, \text{Tr}\left[ 
   \Delta_i G_{ij}^{x} \Delta_j G_{ji}^{y} - \Delta_i G_{ij}^{y} \Delta_j G_{ji}^{x} \right] \\
   D_{ij}^x &= -\frac{1}{2\pi} \text{Im} \int dE \, \text{Tr}\left[ 
   \Delta_i G_{ij}^{y} \Delta_j G_{ji}^{z} - \Delta_i G_{ij}^{z} \Delta_j G_{ji}^{y} \right] \\
   D_{ij}^y &= -\frac{1}{2\pi} \text{Im} \int dE \, \text{Tr}\left[ 
   \Delta_i G_{ij}^{z} \Delta_j G_{ji}^{x} - \Delta_i G_{ij}^{z} \Delta_j G_{ji}^{x} \right]

DMI contributions to the spin Hamiltonian: $\mathcal{H}_{\text{DMI}} = \sum_{ij} \mathbf{D}_{ij} \cdot (\mathbf{S}_i \times \mathbf{S}_j)$

Magnetic Anisotropy Tensor
---------------------------

The **anisotropy tensor** $\mathbf{A}_{ij}$ describes direction-dependent exchange:

.. math::

   A_{ij}^{\alpha\beta} = -\frac{1}{2\pi} \text{Im} \int_{-\infty}^{E_F} dE \, 
   \text{Tr}\left[ \Delta_i G_{ij}^{\alpha}(E) \Delta_j G_{ji}^{\beta}(E) \right]

For $\alpha, \beta \in \{x, y, z\}$, giving a $3 \times 3$ tensor.

**Anisotropy contribution** to spin Hamiltonian:

.. math::

   \mathcal{H}_{\text{ani}} = -\sum_{ij} \sum_{\alpha\beta} A_{ij}^{\alpha\beta} S_i^{\alpha} S_j^{\beta}

On-Site Exchange J₀
-------------------

For the **on-site term** ($i = j$), a modified formula is used:

.. math::

   J_{00} = \frac{1}{\pi} \text{Im} \int_{-\infty}^{E_F} dE \, 
   \text{Tr}\left[ \Delta_i \left( G_{ii}^{\uparrow\uparrow}(E) - G_{ii}^{\downarrow\downarrow}(E) \right) 
   + \Delta_i G_{ii}^{\uparrow\uparrow}(E) \Delta_i G_{ii}^{\downarrow\downarrow}(E) \right]

This represents the energy cost of flipping the spin at site $i$.

Implementation Details
======================

Auxiliary Green's Function Method
----------------------------------

The code uses **auxiliary Green's functions** to evaluate the LKAG formulas. These are constructed from the full Green's function matrices by projecting onto specific spin channels:

**Implementation in** ``calculate_jij_auxgreen()``:

.. code-block:: fortran

   ! Calculate potential functions P for atoms i and j
   call symbolic_atom(iz(i))%p_matrix(pmatrix_i, lmaxi, en%ene)
   call symbolic_atom(iz(j))%p_matrix(pmatrix_j, lmaxj, en%ene)
   
   ! Get auxiliary Green's functions
   call green%auxiliary_gij(green%gij(:,:,:,njij), aux_gij, i, j)
   call green%auxiliary_gij(green%gji(:,:,:,njij), aux_gji, j, i)
   
   ! Construct ΔP = P↑ - P↓
   deltap_i(:,:,:) = pmatrix_i(1:nup, 1:nup, :) - 
                     pmatrix_i(nup+1:2*nup, nup+1:2*nup, :)
   deltap_j(:,:,:) = pmatrix_j(1:nup, 1:nup, :) - 
                     pmatrix_j(nup+1:2*nup, nup+1:2*nup, :)

**Energy integration** using Simpson's rule:

.. code-block:: fortran

   ! Calculate integrand at each energy point
   do nv = 1, en%channels_ldos
      temp1(:,:,nv) = matmul(deltap_i(:,:,nv), aux_gij_up(:,:,nv))
      temp2(:,:,nv) = matmul(deltap_j(:,:,nv), aux_gji_dn(:,:,nv))
      integrand(nv) = imtrace(matmul(temp1(:,:,nv), temp2(:,:,nv)))
   end do
   
   ! Simpson integration up to Fermi level
   call simpson_f(jij, en%ene, en%fermi, en%nv1, integrand, ...)
   jij = jij / (4.0 * pi)  ! Normalize

Angular Decomposition for Tensor Components
--------------------------------------------

For tensor components $J_{ij}^{\alpha\beta}$, the code uses rotation angles to project onto different directions:

**Rotation angles** (from ``calculate_jij_auxgreen``):

.. list-table::
   :header-rows: 1
   :widths: 10 15 15 15 15

   * - Component
     - $\theta$ (rad)
     - $\theta'$ (rad)
     - $\phi$ (rad)
     - $\phi'$ (rad)
   * - $J_{ij}^{xx}$
     - $\pi/2$
     - $\pi/2$
     - $0$
     - $0$
   * - $J_{ij}^{xy}$
     - $\pi/2$
     - $\pi/2$
     - $0$
     - $\pi/2$
   * - $J_{ij}^{xz}$
     - $\pi/2$
     - $0$
     - $0$
     - $0$
   * - $J_{ij}^{yy}$
     - $\pi/2$
     - $\pi/2$
     - $\pi/2$
     - $\pi/2$
   * - $J_{ij}^{zz}$
     - $0$
     - $0$
     - $0$
     - $0$

**Angular integration formula:**

.. math::

   J_{ij}^{\alpha\beta} = \sum_{k} \cos\theta \cos\theta' \, T_{\uparrow\uparrow} 
   + \sin\theta \sin\theta' e^{i(\phi'-\phi)} \, T_{\downarrow\downarrow} + \ldots

where $T$ represents the traces of Green's function products.

Main Procedures
===============

calculate_jij_auxgreen
----------------------

**Purpose:** Calculate Heisenberg exchange using auxiliary Green's function formalism.

**Algorithm:**

1. Loop over all atom pairs $(i,j)$ specified in input
2. For each pair:
   
   a. Get potential functions $P_i^{\sigma}$, $P_j^{\sigma}$
   b. Construct $\Delta P_i = P_i^{\uparrow} - P_i^{\downarrow}$
   c. Get auxiliary Green's functions $G_{ij}$, $G_{ji}$
   d. Calculate 9 tensor components $J_{ij}^{\alpha\beta}$
   e. For $i = j$: calculate on-site $J_{00}$

3. Energy integration: Simpson's rule from $-\infty$ to $E_F$
4. Output to ``exchange.out``

**Formula implemented:**

.. math::

   J_{ij}^{\alpha\beta}(\mathbf{e}_i, \mathbf{e}_j) = -\frac{1}{2\pi} \text{Im} 
   \int_{-\infty}^{E_F} dE \sum_{LL'MM'} 
   \Delta P_{iL}^{\alpha L'} G_{ij,L'M}(E) \Delta P_{jM}^{\beta M'} G_{ji,M'L}(E)

where $L, L', M, M'$ are combined orbital-spin indices.

calculate_exchange_twoindex
----------------------------

**Purpose:** Calculate all exchange interactions (J, D, A) in non-collinear formalism.

**Decomposition into spin-dependent contributions:**

- **CD:** Core-core and down-down cross terms
- **SD:** Spin-dependent terms with opposite spins
- **CC:** Core-core with same spin
- **SC:** Spin-core mixing terms

**Output quantities:**

- ``jijcd, jijsd, jijcc, jijsc`` - J components
- ``dmi(3)`` - DMI vector components
- ``aij(3,3)`` - Full anisotropy tensor

calculate_gilbert_damping
--------------------------

**Purpose:** Calculate Gilbert damping parameters using torque-correlation method.

**Formula:** Based on Kamberský's torque-correlation model (PRM 2, 013801 (2018)):

.. math::

   \alpha_{ij}^{\alpha\beta} = \frac{g}{2\pi M_i} \int_{-\infty}^{E_F} dE \, 
   \text{Tr}\left[ \mathbf{T}_i^{\alpha} A_{ij}(E) \mathbf{T}_j^{\beta\dagger} A_{ji}(E) \right]

where:

- $\mathbf{T}_i^{\alpha}$ is the torque operator on site $i$ in direction $\alpha$
- $A_{ij}(E) = G_{ij}(E) - G_{ij}^{\dagger}(E)$ is the anti-Hermitian part of Green's function
- $M_i$ is the magnetic moment on site $i$
- $g$ is the Landé g-factor

**Implementation details:**

1. Calculate torque operators for each atom type
2. Compute anti-Hermitian parts: $A_{ij} = G_{ij} - G_{ji}^{\dagger}$
3. For each pair and energy:

.. code-block:: text

   temp1 = matmul(T_i^α, A_ij)
   temp2 = matmul(T_j^β†, A_ji)
   temp3 = matmul(temp1, temp2)
   damping^αβ(E) = Re[Tr(temp3)]

4. Integrate over energy up to $E_F$
5. Output damping tensor to ``damping-energy.out``

calculate_jijk
--------------

**Purpose:** Calculate three-body spin-lattice coupling parameters.

**Formula:** Extension to include lattice distortions $\mathbf{u}$:

.. math::

   J_{ijk}^{\alpha\beta\gamma} = \frac{\partial^2 J_{ij}^{\alpha\beta}}{\partial u_k^{\gamma}}

Used for **magnetoelastic coupling** and spin-phonon interactions.

Exchange Interaction Contributions
===================================

First-Order and Second-Order Perturbation
------------------------------------------

The code decomposes exchange into perturbative orders:

**Second-order (SO):**

.. math::

   J_{ij}^{\text{SO}} = J_{ij}^{\text{SD}} + J_{ij}^{\text{SC}}

**First-order (FO):**

.. math::

   J_{ij}^{\text{FO}} = J_{ij}^{\text{SC}} - J_{ij}^{\text{CC}}

These correspond to different orders in the spin-orbit coupling perturbation expansion.

Physical Interpretation
-----------------------

**Heisenberg J:**

- $J_{ij} > 0$: Ferromagnetic coupling (parallel spins favored)
- $J_{ij} < 0$: Antiferromagnetic coupling (antiparallel spins favored)

**DMI $\mathbf{D}_{ij}$:**

- Chiral interaction favoring non-collinear spin textures
- Important for skyrmions, spin spirals
- Typically $\|\mathbf{D}\| \ll \|J\|$ unless strong SOC + inversion breaking

**Anisotropy $\mathbf{A}_{ij}$:**

- Direction-dependent exchange
- $A_{ij}^{zz} - \frac{1}{2}(A_{ij}^{xx} + A_{ij}^{yy})$ gives uniaxial anisotropy
- Important for magnetic domain orientation

Output Files
============

exchange.out
------------

**Format:**

.. code-block:: text

   # i    j    Jij(meV)   Jxx   Jxy   Jxz   Jyx   Jyy   Jyz   Jzx   Jzy   Jzz   Rij(Å)
     1    2    15.243     ...   ...   ...   ...   ...   ...   ...   ...   ...   2.866
     1    3   -2.156      ...   ...   ...   ...   ...   ...   ...   ...   ...   4.050

- Atom indices ``i, j``
- Isotropic exchange $J_{ij}$ in meV
- Full $3 \times 3$ tensor components
- Distance between atoms

dmi.out
-------

DMI vectors for each pair:

.. code-block:: text

   # i    j    Dx(meV)   Dy(meV)   Dz(meV)   |D|(meV)   Rij(Å)
     1    2    0.234     -0.156     1.234     1.267      2.866

anisotropy.out
--------------

Full anisotropy tensors:

.. code-block:: text

   # i    j    Axx   Axy   Axz   Ayx   Ayy   Ayz   Azx   Azy   Azz   (meV)
     1    1    2.3   0.0   0.0   0.0   2.3   0.0   0.0   0.0   2.8

damping-energy.out
------------------

Gilbert damping as function of energy:

.. code-block:: text

   # E-EF(eV)   α_xx   α_xy   α_xz   ...   α_zz
     -0.500     0.012  0.000  0.000  ...   0.015

Integration and Usage
=====================

Workflow
--------

1. **SCF convergence** must be completed first
2. **Set exchange calculation flag** in input:

.. code-block:: fortran

   &control
      calc_exchange = .true.
      njij = 10           ! Number of pairs
   /

3. **Specify atom pairs** in ``clust`` file or via namelist
4. **Run calculation:** Exchange computed during post-processing
5. **Output files** written to working directory

Connection to Spin Dynamics
----------------------------

Exchange parameters feed directly into spin dynamics module (``spin_dynamics.f90``):

**Classical Heisenberg Hamiltonian:**

.. math::

   \mathcal{H} = -\sum_{ij} J_{ij} \mathbf{S}_i \cdot \mathbf{S}_j 
   + \sum_{ij} \mathbf{D}_{ij} \cdot (\mathbf{S}_i \times \mathbf{S}_j)
   - \sum_{ij\alpha\beta} A_{ij}^{\alpha\beta} S_i^{\alpha} S_j^{\beta}

See :ref:`theory/spin_dynamics` for atomistic spin dynamics using these parameters.

Performance Considerations
==========================

Computational Cost
------------------

**Scaling:**

- $\mathcal{O}(N_{\text{pairs}} \times N_E \times N_L^3)$
- $N_{\text{pairs}}$: Number of exchange pairs
- $N_E$: Energy mesh points (~200-500)
- $N_L$: Orbital basis size ($\sim (l_{\max}+1)^2 \approx 9-16$)

**Typical timings:**

- 10 pairs, 400 energy points, $l_{\max}=2$: ~5-10 minutes
- 100 pairs: ~1-2 hours

Parallelization
---------------

**MPI parallelization** over atom pairs:

.. code-block:: fortran

   #ifdef USE_MPI
   do njij_glob = start_pair, end_pair
      njij = g2l_pair_map(njij_glob)
      ! Calculate exchange for this pair
   end do
   #endif

Good scaling up to $N_{\text{cores}} = N_{\text{pairs}}$.

Numerical Stability
-------------------

**Tips for accurate results:**

1. **Dense energy mesh:** Use 300-500 energy points
2. **Converged SCF:** Tight convergence ($\Delta Q < 10^{-5}$)
3. **Sufficient recursion:** ``ll_gf`` $\geq$ 50 for good real-space decay
4. **Check Im[Tr]:** Should be well-defined (not oscillating wildly)

Provenance
==========

**Implementation:**

- ``source/exchange.f90`` - Main exchange module (1917 lines)
- ``source/green.f90`` - Auxiliary Green's function construction
- ``source/hamiltonian.f90`` - Torque operators for damping

**Key procedures:**

- ``calculate_jij_auxgreen()`` - LKAG formula with auxiliary GFs
- ``calculate_exchange_twoindex()`` - Full non-collinear exchange
- ``calculate_gilbert_damping()`` - Torque-correlation damping
- ``calculate_dij()`` - DMI calculation
- ``calculate_aij()`` - Anisotropy calculation
- ``calculate_jijk()`` - Spin-lattice coupling

**Authors:**

- Angela Klautau, Ramon Cardias, Lucas P. Campagna
- S. Frota-Pessôa, Pascoal R. Peduto
- Anders Bergman, Ivan P. Miranda
- S. B. Legoas, H. M. Petrilli

References
==========

1. **LKAG formalism:** A. I. Liechtenstein, M. I. Katsnelson, V. P. Antropov, and V. A. Gubanov, 
   *J. Magn. Magn. Mater.* **67**, 65 (1987)

2. **Auxiliary Green's function:** M. Pajda, J. Kudrnovský, I. Turek, V. Drchal, and P. Bruno,
   *Phys. Rev. B* **64**, 174402 (2001)

3. **DMI from LMTO:** L. M. Sandratskii, *Phys. Rev. B* **96**, 024450 (2017)

4. **Gilbert damping:** H. Ebert, S. Mankovsky, D. Ködderitzsch, and P. J. Kelly,
   *Phys. Rev. Lett.* **107**, 066603 (2011); 
   I. P. Miranda, A. B. Klautau, et al., *Phys. Rev. Materials* **2**, 013801 (2018)

5. **Torque method:** V. Kamberský, *Czech. J. Phys. B* **26**, 1366 (1976)

See Also
========

- :ref:`theory/spin_dynamics` - Atomistic spin dynamics using exchange parameters
- :ref:`theory/green_functions` - Green's function formalism
- :ref:`user_guide/examples` - Example: Exchange calculation for Fe bcc
- :ref:`keywords/output_options` - Output control for exchange files
