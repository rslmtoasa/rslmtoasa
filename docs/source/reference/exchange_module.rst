.. _reference/exchange_module:

=======================================
Exchange Interactions Module
=======================================

Overview
========

The **exchange module** (``exchange.f90``) calculates magnetic exchange interactions using the **Liechtenstein-Katsnelson-Antropov-Guslienko (LKAG) formalism** within the LMTO-ASA framework. This module computes:

- Isotropic Heisenberg exchange interactions :math:`J_{ij}`
- Dzyaloshinskii-Moriya interactions (DMI) :math:`\mathbf{D}_{ij}`
- Magnetic anisotropy tensors :math:`\mathbf{A}_{ij}`
- Gilbert damping parameters :math:`\alpha_{ij}`
- Moment of inertia tensors

``calculate_jijk`` (spin-lattice coupling :math:`J_{ijk}`) is a type-bound
procedure with no caller in the production driver.

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

**Isotropic Heisenberg exchange** :math:`J_{ij}`:

.. math::

   J_{ij} = -\frac{1}{4\pi} \text{Im} \int_{-\infty}^{E_F} dE \, 
   \text{Tr}\left[ \Delta_i G_{ij}^{\uparrow\uparrow}(E) 
   \Delta_j G_{ji}^{\downarrow\downarrow}(E) \right]

where:

- :math:`\Delta_i = P_i^{\uparrow} - P_i^{\downarrow}` is the spin-split potential function
- :math:`P_i^{\sigma}` are the potential functions for spin :math:`\sigma`
- :math:`G_{ij}^{\sigma\sigma'}(E)` is the intersite Green's function from site :math:`j` to site :math:`i`
- :math:`E_F` is the Fermi energy
- Tr denotes the trace over orbital indices

DMI contributions to the spin Hamiltonian: :math:`\mathcal{H}_{\text{DMI}} = \sum_{ij} \mathbf{D}_{ij} \cdot (\mathbf{S}_i \times \mathbf{S}_j)`

**Anisotropy contribution** to spin Hamiltonian:

.. math::

   \mathcal{H}_{\text{ani}} = -\sum_{ij} \sum_{\alpha\beta} A_{ij}^{\alpha\beta} S_i^{\alpha} S_j^{\beta}

On-Site Exchange J₀
-------------------

For the **on-site term** (:math:`i = j`), ``calculate_jij_auxgreen`` uses:

.. math::

   J_{00} = \frac{1}{\pi} \text{Im} \int_{-\infty}^{E_F} dE \, 
   \text{Tr}\left[ \Delta_i \left( G_{ii}^{\uparrow\uparrow}(E) - G_{ii}^{\downarrow\downarrow}(E) \right) 
   + \Delta_i G_{ii}^{\uparrow\uparrow}(E) \Delta_i G_{ii}^{\downarrow\downarrow}(E) \right]

This represents the energy cost of flipping the spin at site :math:`i`.

Implemented Kernels
===================

``calculate_exchange`` evaluates three separate kernels at every energy point
of the mesh (``dGdG_Jnc``, ``dGdG_Dnc``, ``dGdG_Anc``). With
:math:`\Delta_i` the ``d_matrix`` of the type of atom :math:`i`,
:math:`G^0` the charge part (``Ginmag``, ``Gjnmag``) and :math:`G^\alpha`
(:math:`\alpha = x, y, z`) the spin parts (``Gix``... ``Gjz``) of the
intersite Green function:

.. math::

   J_{ij}(E) &= \text{Im}\,\text{Tr}\Big[ \Delta_i G_i^0 \Delta_j G_j^0
              - \sum_\alpha \Delta_i G_i^\alpha \Delta_j G_j^\alpha \Big] \\
   D_{ij}^\alpha(E) &= \text{Re}\,\text{Tr}\Big[ \Delta_i G_i^0 \Delta_j G_j^\alpha
              - \Delta_j G_j^0 \Delta_i G_i^\alpha \Big] \\
   A_{ij}^{\alpha\beta}(E) &= \text{Im}\,\text{Tr}\Big[ \tfrac{1}{2}\big(
              \Delta_i G_i^\alpha \Delta_j G_j^\beta
              + \Delta_j G_j^\alpha \Delta_i G_i^\beta \big) \Big]

Each kernel is integrated with ``simpson_f`` up to the Fermi energy with
the Fermi weight at ``en%temperature`` and multiplied by
:math:`10^3/(4\pi)`. :math:`J` is not a trace of :math:`A`.
``jtens.out`` holds :math:`J\,\mathbb{1}` plus the skew-symmetric matrix of
:math:`\mathbf{D}` plus :math:`A`.

calculate_exchange_twoindex
---------------------------

Splits each kernel into four parts built from the two-index Green functions
(``g00ij``/``g00ji``, ``g01ij``/``g01ji``, ``gx0ij``... ``gz1ji``):

- **CD:** charge-density part (``g00ij``, ``g00ji``)
- **CC:** charge-current part (``g01ij``, ``g01ji``)
- **SD:** spin-density part (``gx0ij``, ``gy0ij``, ``gz0ij`` and their ``ji`` partners, summed)
- **SC:** spin-current part (``gx1ij``, ``gy1ij``, ``gz1ij`` and their ``ji`` partners, summed)

and combines them into second-order (SO) and first-order (FO) results:

.. math::

   J^{\text{SO}} &= J^{\text{CD}} - J^{\text{SD}} + J^{\text{CC}} - J^{\text{SC}}, &
   J^{\text{FO}} &= J^{\text{CD}} + J^{\text{SD}} - J^{\text{CC}} - J^{\text{SC}} \\
   D^{\text{SO}} &= 2\,(D^{\text{SC}} + D^{\text{CC}}), &
   D^{\text{FO}} &= 2\,(D^{\text{SC}} - D^{\text{CC}}) \\
   A^{\text{SO}} &= A^{\text{SD}} + A^{\text{SC}}, &
   A^{\text{FO}} &= -A^{\text{SD}} + A^{\text{SC}}

Here the D parts use the real part of the trace, J and A the imaginary part.

calculate_jij_auxgreen
----------------------

LKAG formula with auxiliary Green functions (``green%auxiliary_gij``) and
potential-function differences :math:`\Delta P = P^\uparrow - P^\downarrow`.
It has no caller in ``post_processing='exchange'`` and writes only to
standard output. For tensor components it uses rotation angles:

.. list-table::
   :header-rows: 1
   :widths: 10 15 15 15 15

   * - Component
     - :math:`\theta` (rad)
     - :math:`\theta'` (rad)
     - :math:`\phi` (rad)
     - :math:`\phi'` (rad)
   * - :math:`J_{ij}^{xx}`
     - :math:`\pi/2`
     - :math:`\pi/2`
     - :math:`0`
     - :math:`0`
   * - :math:`J_{ij}^{xy}`
     - :math:`\pi/2`
     - :math:`\pi/2`
     - :math:`0`
     - :math:`\pi/2`
   * - :math:`J_{ij}^{xz}`
     - :math:`\pi/2`
     - :math:`0`
     - :math:`0`
     - :math:`0`
   * - :math:`J_{ij}^{yy}`
     - :math:`\pi/2`
     - :math:`\pi/2`
     - :math:`\pi/2`
     - :math:`\pi/2`
   * - :math:`J_{ij}^{zz}`
     - :math:`0`
     - :math:`0`
     - :math:`0`
     - :math:`0`

**Angular integration formula:**

.. math::

   J_{ij}^{\alpha\beta} = \sum_{k} \cos\theta \cos\theta' \, T_{\uparrow\uparrow} 
   + \sin\theta \sin\theta' e^{i(\phi'-\phi)} \, T_{\downarrow\downarrow} + \ldots

where :math:`T` represents the traces of Green's function products.

calculate_gilbert_damping
-------------------------

Torque-correlation Gilbert damping (``exchange_dynamics.f90``; see
PRM 2, 013801 (2018)). Runs only when ``control%nsp == 2``; otherwise it
does nothing. Called when ``&calculation do_damping = .true.``.

For each pair, with torque matrices :math:`T` from
``hamiltonian%torque_operator_collinear`` and
:math:`A_{ij} = G_{ij} - G_{ji}^\dagger`:

.. math::

   \alpha_{ij}^{\alpha\beta} = -\frac{1}{2\pi S_i}\,
   \text{Re}\,\text{Tr}\left[ T_i^{\alpha} A_{ij}(E) T_j^{\beta\dagger} A_{ji}(E) \right]_{E = E_k}

where :math:`E_k` is the energy-mesh point nearest to :math:`E_F` (no energy
integral) and :math:`S_i = \sum_l (q_{l}^{\uparrow} - q_{l}^{\downarrow})`
is taken from ``potential%ql`` of atom :math:`i`. No g factor is applied.

calculate_moment_of_inertia
---------------------------

Torque-correlation moment-of-inertia kernel (Sci. Rep. 7, 931 (2017)),
called when ``&calculation do_inertia = .true.``; requires ``nsp = 2``.

calculate_jijk
--------------

Spin-lattice coupling :math:`J_{ijk}`. The routine has no caller in the
production driver; with ``njijk > 0`` ``post_processing='exchange'`` runs no
exchange routine. It writes only to standard output.

Physical Interpretation
-----------------------

**Heisenberg J:**

- :math:`J_{ij} > 0`: Ferromagnetic coupling (parallel spins favored)
- :math:`J_{ij} < 0`: Antiferromagnetic coupling (antiparallel spins favored)

**DMI** :math:`\mathbf{D}_{ij}`:

- Chiral interaction favoring non-collinear spin textures
- Important for skyrmions, spin spirals
- Typically :math:`\|\mathbf{D}\| \ll \|J\|` unless strong SOC + inversion breaking

**Anisotropy** :math:`\mathbf{A}_{ij}`:

- Direction-dependent exchange
- :math:`A_{ij}^{zz} - \frac{1}{2}(A_{ij}^{xx} + A_{ij}^{yy})` gives uniaxial anisotropy
- Important for magnetic domain orientation

Output Files
============

J, D and A are :math:`10^3/(4\pi)` times an energy integral in Ry (mRy
scale); no conversion to meV is applied. The files have no header line.
Each record is:

.. code-block:: text

   type_i  type_j  dx  dy  dz  <values>  distance

where ``type_i``, ``type_j`` are the atom types (``iz``), ``dx dy dz`` is
:math:`\mathbf{r}_j - \mathbf{r}_i` and ``distance`` its norm.

- From ``calculate_exchange``: ``jij.out`` (J), ``dij.out`` (3 D
  components), ``aij.out`` (9 A components), ``jtens.out``.
- From ``calculate_exchange_twoindex``: ``jijso.out``, ``jijfo.out``,
  ``dijso.out``, ``dijfo.out``, ``aijso.out``, ``aijfo.out``,
  ``jtensso.out``, ``jtensfo.out``, ``jijparts.out``, ``dijparts.out``,
  ``aijparts.out``.
- From ``calculate_gilbert_damping``: ``alldampings.out`` (per pair, atom
  indices i, j) and ``damping-energy.out`` (energy-resolved sum; first
  column :math:`E - E_F` in Ry; sum over all listed pairs, each with its
  own prefactor, so at :math:`E_F` it equals the sum of the ``alldampings.out`` rows).

Integration and Usage
=====================

1. **SCF convergence** must be completed first
2. Select the exchange post-processing and the pairs:

.. code-block:: fortran

   &calculation
      post_processing = 'exchange'
   /
   &lattice
      njij = 2
      ijpair(1,:) = 1, 2
      ijpair(2,:) = 1, 3
   /

3. ``post_processing_exchange`` calls ``calculate_exchange`` and
   ``calculate_exchange_twoindex`` when ``njij /= 0`` and ``njijk == 0``.

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

- :math:`\mathcal{O}(N_{\text{pairs}} \times N_E \times N_L^3)`
- :math:`N_{\text{pairs}}`: Number of exchange pairs
- :math:`N_E`: Energy mesh points (~200-500)
- :math:`N_L`: Orbital basis size (:math:`\sim (l_{\max}+1)^2 \approx 9-16`)

**Typical timings:**

- 10 pairs, 400 energy points, :math:`l_{\max}=2`: ~5-10 minutes
- 100 pairs: ~1-2 hours

Parallelization
---------------

**MPI parallelization** over atom pairs:

.. code-block:: fortran

   do njij_glob = start_atom, end_atom
      njij = g2l_map(njij_glob)
      ! Calculate exchange for this pair
   end do

Results are combined with ``MPI_ALLREDUCE``.

Numerical Stability
-------------------

**Tips for accurate results:**

1. **Dense energy mesh:** Use 300-500 energy points
2. **Converged SCF:** Tight convergence (:math:`\Delta Q < 10^{-5}`)
3. **Sufficient recursion:** ``lld`` / ``llsp`` large enough for good real-space decay
4. **Check Im[Tr]:** Should be well-defined (not oscillating wildly)

Provenance
==========

**Implementation:**

- ``source/exchange.f90`` - Main exchange module
- ``source/exchange_dynamics.f90`` - Gilbert damping and moment of inertia
- ``source/green.f90`` - Auxiliary Green's function construction
- ``source/hamiltonian.f90`` - Torque operators for damping

**Key procedures:**

- ``calculate_exchange()`` - J, D, A kernels
- ``calculate_exchange_twoindex()`` - CD/CC/SD/SC decomposition
- ``calculate_gilbert_damping()`` - Torque-correlation damping
- ``calculate_moment_of_inertia()`` - Torque-correlation inertia kernel
- ``calculate_jij_auxgreen()`` - LKAG formula with auxiliary GFs (no caller)
- ``calculate_jijk()`` - Spin-lattice coupling (no caller)

``calculate_dij()`` and ``calculate_aij()`` are empty.

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
