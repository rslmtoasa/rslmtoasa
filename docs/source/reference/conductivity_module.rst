.. _reference/conductivity_module:

=======================================
Conductivity Module
=======================================

Overview
========

The **conductivity module** (``conductivity.f90``) calculates the **electrical conductivity tensor** :math:`\sigma_{\alpha\beta}(\omega)` using the **Kubo-Bastin formalism** combined with **Chebyshev polynomial expansion** within the LMTO-ASA framework.

**Key features:**

- Optical conductivity :math:`\sigma(\omega)` as function of frequency
- DC conductivity :math:`\sigma(0)`
- Full :math:`3 \times 3` conductivity tensor (anisotropic transport)
- Orbital-resolved conductivity contributions
- Chebyshev recursion for efficient computation
- Stochastic trace evaluation for large systems
- Real and imaginary parts (dissipative + reactive)

**Applications:**

- Electronic transport properties
- Optical response
- Hall conductivity
- Magneto-transport
- Thermoelectric properties (with extension)

Module Structure
================

Type: conductivity
------------------

.. code-block:: fortran

   type conductivity
      ! Pointers to system objects
      class(control), pointer :: control
      class(lattice), pointer :: lattice
      class(self), pointer :: self
      class(energy), pointer :: en
      class(recursion), pointer :: recursion
      
      ! Conductivity parameters
      integer :: ll_cond  ! Recursion depth for conductivity
      
      ! Pre-computed coefficients
      complex(rp), dimension(:,:,:), allocatable :: gamma_nm
         ! Dimension: (channels_ldos+10, cond_ll, cond_ll)
         ! Chebyshev expansion coefficients
   end type

Theory: Kubo-Bastin Formalism
==============================

Linear Response and Kubo Formula
---------------------------------

The **Kubo formula** for electrical conductivity relates the current response to an applied electric field via linear response theory:

.. math::

   \sigma_{\alpha\beta}(\omega) = \frac{1}{\Omega} \lim_{\omega' \to \omega + i\eta} 
   \int_{-\infty}^{\infty} dt \, e^{i\omega' t} 
   \langle [\hat{j}_\alpha(t), \hat{j}_\beta(0)] \rangle

where:

- :math:`\Omega` is the system volume
- :math:`\hat{j}_\alpha` is the current operator in direction :math:`\alpha`
- :math:`\langle \ldots \rangle` denotes quantum/thermal average
- :math:`\eta \to 0^+` is a positive infinitesimal

Velocity Operator and Current
------------------------------

The **current operator** is related to the **velocity operator** :math:`\hat{v}`:

.. math::

   \hat{j}_\alpha = -e \sum_i \hat{v}_{i,\alpha} = -e \sum_i \frac{\partial \hat{H}}{\partial \mathbf{k}_\alpha}

In the tight-binding LMTO formalism:

.. math::

   \hat{v}_\alpha = \frac{1}{i\hbar}[\hat{H}, \hat{r}_\alpha] 
   = \frac{1}{\hbar} \sum_{nm} \langle n|\frac{\partial \hat{H}}{\partial k_\alpha}|m\rangle |n\rangle\langle m|

Kubo-Bastin Formula
--------------------

The **Kubo-Bastin formula** expresses conductivity in terms of Green's functions:

.. math::

   \sigma_{\alpha\beta}(\omega) = \frac{e^2\hbar}{\pi\Omega} \int_{-\infty}^{\mu} dE \, 
   \text{Tr}\left[ \hat{v}_\alpha \text{Im}[G^+(E)] \hat{v}_\beta \text{Im}[G^+(E + \hbar\omega)] \right]

where:

- :math:`G^+(E) = (E - \hat{H} + i\eta)^{-1}` is the retarded Green's function
- :math:`\mu` is the chemical potential (Fermi energy)
- :math:`\text{Im}[G^+] = A(E)/(2\pi)` is the spectral function

For **DC conductivity** (:math:`\omega = 0`):

.. math::

   \sigma_{\alpha\beta}(0) = \frac{e^2\hbar}{\pi\Omega} \int_{-\infty}^{\mu} dE \, 
   \text{Tr}\left[ \hat{v}_\alpha A(E) \hat{v}_\beta A(E) \right]

where :math:`A(E) = -2\pi \text{Im}[G^+(E)]` is the spectral function.

Chebyshev Polynomial Expansion
===============================

Kernel Polynomial Method (KPM)
-------------------------------

To avoid expensive direct evaluation, the **Chebyshev expansion** approximates the Green's function:

**Energy rescaling:** First, map the Hamiltonian spectrum :math:`[E_{\min}, E_{\max}]` to :math:`[-1, 1]`:

.. math::

   \tilde{H} = \frac{H - b}{a}, \quad a = \frac{E_{\max} - E_{\min}}{2 - \epsilon}, \quad b = \frac{E_{\max} + E_{\min}}{2}

where :math:`\epsilon \sim 0.3` provides margin.

**Chebyshev expansion** of Green's function:

.. math::

   G^+(\tilde{E}) = \sum_{n=0}^{N-1} \mu_n T_n(\tilde{E}) \frac{g_n}{\sqrt{1 - \tilde{E}^2}}

where:

- :math:`T_n(x)` are Chebyshev polynomials of the first kind
- :math:`\mu_n = 2 - \delta_{n,0}` are normalization factors
- :math:`g_n` are kernel damping factors (smoothing; the code uses the Lorentz kernel, see below)
- :math:`N` is the expansion order (recursion depth)

**Moments** :math:`\mu_n` are calculated via recursion:

.. math::

   \mu_n = \langle \psi | T_n(\tilde{H}) | \psi \rangle

using the **three-term recurrence**:

.. math::

   \begin{aligned}
   T_0(x) &= 1 \\
   T_1(x) &= x \\
   T_{n+1}(x) &= 2x T_n(x) - T_{n-1}(x)
   \end{aligned}

Gamma Function :math:`\Gamma_{nm}`

-------------------------------

The conductivity involves **products of two Green's functions**, leading to the :math:`\Gamma_{nm}` function:

.. math::

   \Gamma_{nm}(E) = \frac{1}{(1 - \tilde{E}^2)^2} 
   \left[ C_n(\tilde{E}) T_m(\tilde{E}) + \bar{C}_m(\tilde{E}) T_n(\tilde{E}) \right] g_n g_m w_n w_m

where, for orders :math:`n, m = 0, \ldots, N-1`:

- :math:`C_n(x) = \left(x - i n \sqrt{1 - x^2}\right) e^{i n \arccos x}`
- :math:`\bar{C}_m(x) = \left(x + i m \sqrt{1 - x^2}\right) e^{-i m \arccos x}`
- :math:`w_0 = 1/2`, :math:`w_n = 1` otherwise
- :math:`g_n` are Lorentz kernel factors

**Lorentz kernel** (``lorentz_kernel`` with :math:`\lambda = 6`):

.. math::

   g_n = \frac{\sinh\left[\lambda (1 - n/N)\right]}{\sinh \lambda}

The run stops if ``cond_ll`` < 2 or any scaled energy lies outside
:math:`(-1, 1)`.

**Implementation** in ``calculate_gamma_nm()``:

.. code-block:: text

   ! Energy rescaling
   a = (e_max - e_min) / (2 - 0.3)
   b = (e_max + e_min) / 2
   wscale(:) = (en%ene(:) - b) / a
   
   ! Chebyshev polynomials
   T_0 = 1.0
   T_1 = wscale
   do n = 2, cond_ll
      T_n = 2 * wscale * T_{n-1} - T_{n-2}
   end do
   
   ! Cn coefficients
   acos_x = acos(wscale)
   sqrt_term = sqrt(1.0 - wscale**2)
   cn(:,n) = (wscale - i*(n-1)*sqrt_term) * exp( i*(n-1)*acos_x)
   cm(:,n) = (wscale + i*(n-1)*sqrt_term) * exp(-i*(n-1)*acos_x)
   
   ! Gamma_nm matrix
   do n = 1, cond_ll
      do m = 1, cond_ll
         gamma_nm(:,n,m) = (cn(:,n) * T_m + cm(:,m) * T_n) / (1 - wscale**2)**2
         gamma_nm(:,n,m) = gamma_nm(:,n,m) * g_n * g_m * w_n * w_m
      end do
   end do

Stochastic Trace Evaluation
----------------------------

For large systems, exact trace evaluation :math:`\text{Tr}[\ldots]` is prohibitive. Instead, use **stochastic trace estimation**:

.. math::

   \text{Tr}[A] \approx \frac{1}{N_r} \sum_{r=1}^{N_r} \langle \xi_r | A | \xi_r \rangle

where :math:`|\xi_r\rangle` are **random vectors** with elements drawn from:

- Complex normal distribution: :math:`\xi_{i,r} \sim \mathcal{N}(0, 1) + i\mathcal{N}(0, 1)`
- Or: Random phase: :math:`\xi_{i,r} = e^{i\phi_r}`

**Variance** decreases as :math:`\mathcal{O}(1/\sqrt{N_r})`, so typically :math:`N_r = 10-50` random vectors suffice.

Conductivity Tensor Calculation
================================

Full Tensor Formula
-------------------

The **conductivity tensor** :math:`\sigma_{\alpha\beta}(E)` has 9 components (for 3D):

.. math::

   \sigma_{\alpha\beta}(E) = \frac{16}{\pi \Delta E^2} \sum_{n,m}^{N} 
   \Gamma_{nm}(E) \mu_{nm}^{\alpha\beta}

where:

- :math:`\mu_{nm}` are the Chebyshev moments of the velocity operator products
- :math:`\Delta E = E_{\max} - E_{\min}` is the energy span

The raw ``cond_*.out`` files carry only the factor :math:`16/(\pi \Delta E^2)`.
:math:`e^2/\hbar` and the cell volume or area enter only in
``write_physical_response`` (``response_*.out``, ``response_units.txt``).

**Velocity operator matrix elements:**

.. math::

   \langle n | \hat{v}_\alpha | m \rangle = \frac{1}{i\hbar} \langle n | [H, r_\alpha] | m \rangle

In LMTO, these are computed from the Hamiltonian structure function derivatives.

Orbital and Component Decomposition
------------------------------------

The trace can be decomposed:

1. **By orbital character:**

.. math::

   \sigma_{\alpha\beta}^{(l)}(E) = \sum_{m=-l}^{l} \langle l,m | \hat{v}_\alpha G \hat{v}_\beta G | l,m \rangle

giving :math:`s`, :math:`p`, :math:`d` contributions.

2. **By atomic type:**

.. math::

   \sigma_{\alpha\beta}^{(\text{type})}(E) = \sum_{i \in \text{type}} \langle i | \ldots | i \rangle

3. **By energy range:** Integrate from :math:`E_{\min}` to variable upper limit :math:`E`

Implementation Details
======================

Algorithm Overview
------------------

**Main procedure:** ``calculate_conductivity_tensor()``

1. **Precompute** :math:`\Gamma_{nm}`: Call ``calculate_gamma_nm()``
2. **Stochastic recursion:** Compute Chebyshev moments :math:`\mu_{nm}^{\alpha\beta}` via recursion module
3. **Construct integrand:**

.. code-block:: text

   integrand(α,β,E) = factor * sum_{n,m} gamma_nm(E,n,m) * mu_nm(α,β,n,m)

4. **Simpson integration:**

.. code-block:: text

   do E = E_min to E_max
      call simpson_f(σ_αβ(E), wscale, E, nv1, integrand_real(:), ...)
   end do

5. **Output results**

Per-Type Calculation
--------------------

``cond_calctype`` takes ``'per_type'`` (traces per atomic type) or
``'random_vec'`` (random starting vectors). Per type
(``cond_calctype = 'per_type'``):

.. code-block:: fortran

   do ntype = 1, lattice%ntype
      ! Compute stochastic trace restricted to atoms of type ntype
      integrand_type(:,:,:) = factor * sum gamma_nm * mu_nm_type
      ! Output to <element>_cond.out
   end do

Produces separate files for each element (e.g., ``Fe_cond.out``, ``Pt_cond.out``).

Numerical Parameters
====================

Key Control Parameters
----------------------

**In namelist &control:**

- ``cond_ll`` (default: 200): Chebyshev expansion order :math:`N`
  
  - Larger :math:`N` -> better energy resolution, smoother curves
  - Typical: 50-200
  - Convergence: check :math:`\sigma(E)` vs :math:`N`

- ``cond_calctype`` (default: 'per_type'): Calculation mode
  
  - ``'per_type'``: Per-atomic-species conductivity
  - ``'random_vec'``: With random vectors (stochastic)

- ``random_vec_num`` (default: 1): Number of stochastic vectors
  
  - Typical: 10-50
  - More vectors -> better statistics, slower

Convergence Criteria
--------------------

**Check the following for convergence:**

1. **Recursion depth** ``cond_ll``:

.. code-block:: bash

   # Run with different cond_ll values
   cond_ll = 50, 100, 150, 200
   # Check if σ(EF) stabilizes

2. **Energy mesh** ``channels_ldos``:

   - Default: 6000 points (``calctype`` B, S, L), 3000 (I)
   - Must resolve features near :math:`E_F`

3. **Random vectors** ``random_vec_num``:

.. code-block:: bash

   # Check variance in repeated runs
   σ_avg +- σ_std  # Should be < 5% for good statistics

Main Procedures
===============

calculate_gamma_nm
------------------

**Purpose:** Precompute the :math:`\Gamma_{nm}(E)` coefficients for Chebyshev expansion.

**Inputs:**

- Energy mesh ``en%ene(:)``
- Expansion order ``cond_ll``
- Energy bounds ``e_min, e_max``

**Algorithm:**

1. Rescale energies: :math:`\tilde{E} = (E - b)/a`
2. Compute Chebyshev polynomials :math:`T_n(\tilde{E})` recursively
3. Compute :math:`C_n(E)` coefficients
4. Apply Lorentz kernel :math:`g_n`
5. Construct :math:`\Gamma_{nm}(E)` matrix

**Output:**

- ``gamma_nm(energy, n, m)`` - Complex array of Chebyshev coefficients

calculate_conductivity_tensor
------------------------------

**Purpose:** Compute full :math:`3 \times 3` conductivity tensor :math:`\sigma_{\alpha\beta}(E)`.

**Algorithm:**

1. Loop over energy points :math:`E`
2. For each :math:`E`:

   a. Sum over Chebyshev indices: :math:`\sum_{nm} \Gamma_{nm}(E) \mu_{nm}^{\alpha\beta}`
   b. Multiply by prefactor: :math:`16/(\pi \Delta E^2)`
   c. Simpson integrate from :math:`E_{\min}` to :math:`E`: :math:`\sigma(E) = \int^E` 

3. Output:

   - Real part: Dissipative conductivity (absorption)
   - Imaginary part: Reactive conductivity (dispersion)

**Outputs:**

- Total: ``cond_total.out``
- Orbital-resolved: ``cond_total_orb_real.out``, ``cond_total_orb_im.out``
- Per-type: ``<element>_cond.out``

Output Files
============

cond_*.out (raw kernels)
------------------------

``cond_total.out``, ``cond_total_orb_real.out``, ``cond_total_orb_im.out``
and, per type, ``<element>_cond.out``, ``<element>_cond_orb_real.out``,
``<element>_cond_orb_im.out``. No header line. Column 1 is :math:`E - E_F`
in Ry. ``cond_total.out`` then has Re and Im of :math:`\sigma(E)`; the
``_orb_`` files have 18 spin-orbital columns. Values carry only the factor
:math:`16/(\pi \Delta E^2)` (no volume, no physical constants).

response_*.out (physical units)
-------------------------------

``response_total.out``, ``response_orbitals_real.out``,
``response_orbitals_imag.out`` and, per type, ``<element>_response.out``.
Column 1 is :math:`E - E_F` in eV; one header line names the units (S/cm
for the charge conductivity). ``response_units.txt`` records the
normalisation.

Physical Interpretation
=======================

DC Conductivity
---------------

**At Fermi level** :math:`E = E_F`, :math:`\sigma(E_F)` gives **DC conductivity**:

.. math::

   \sigma_{DC} = \sigma(E_F, \omega=0)

**Units:** ``response_total.out`` is in S/cm. ``write_physical_response``
multiplies by :math:`e^2/\hbar \times 10^{8}` and by :math:`1/\Omega`, with
:math:`\Omega` the cell volume in Å\ :sup:`3` (or area in Å\ :sup:`2` for
``response_geometry = 'sheet'``).

Optical Conductivity
--------------------

**For finite frequency** :math:`\omega \neq 0`, the **optical conductivity** describes:

- Absorption spectrum
- Dielectric function via Kramers-Kronig
- Optical transitions

**Real part** :math:`\text{Re}[\sigma(\omega)]`:

- Absorption/dissipation
- Peaks at interband transitions

**Imaginary part** :math:`\text{Im}[\sigma(\omega)]`:

- Dispersion/reactive response
- Related to plasma frequency

Drude Model Comparison
----------------------

At low frequencies, metallic systems follow **Drude model**:

.. math::

   \sigma(\omega) = \frac{\sigma_0}{1 - i\omega\tau}

where :math:`\tau` is the scattering time. From LMTO:

.. math::

   \sigma_0 = \frac{ne^2\tau}{m^*}

can be extracted by fitting near :math:`E_F`.

Integration and Usage
=====================

Input Namelist (&conductivity)
-------------------------------

The ``&conductivity`` namelist has three keys.

- ``ll_cond`` (integer): read but not used by the code. The expansion order
  is ``cond_ll`` in ``&control``.
- ``response_geometry`` (default ``'bulk'``): ``'bulk'`` or ``'sheet'``;
  any other value stops the run.
- ``response_measure`` (default 0): cell volume (Å\ :sup:`3`) or area
  (Å\ :sup:`2`) used for the ``response_*.out`` normalisation; 0 derives it
  from the lattice vectors and ``alat``.

**Example:**

.. code-block:: fortran

   &conductivity
     response_geometry = 'bulk'
   /

**See also:** ``cond_ll``, ``linear_in``, ``linear_out``, ``cond_calctype`` in
:ref:`keywords/control_parameters`.

Input Parameters
----------------

**Minimal setup:**

.. code-block:: fortran

   &calculation
      post_processing = 'conductivity'
   /
   &control
      cond_ll = 200              ! Chebyshev expansion order
      cond_calctype = 'per_type' ! or 'random_vec'
   /

**Advanced:**

.. code-block:: fortran

   &energy
      channels_ldos = 500       ! Energy mesh points
      energy_min = -1.0         ! absolute energy, Ry
      energy_max = 0.5          ! absolute energy, Ry
   /

Workflow
--------

1. **Complete SCF calculation** first
2. **Enable conductivity:** Set ``&calculation post_processing = 'conductivity'``
3. **Run calculation:** Conductivity computed in post-processing
4. **Check convergence:** Vary ``cond_ll``, check :math:`\sigma(E_F)` stability
5. **Analyze output:** Plot :math:`\sigma(E)`, identify features

Example: Conductivity of Cu
----------------------------

.. code-block:: fortran

   &calculation
      post_processing = 'conductivity'
   /
   &control
      calctype = 'S'
      cond_ll = 150
      cond_calctype = 'per_type'
   /
   &energy
      channels_ldos = 400
      energy_min = -8.0
      energy_max = 3.0
   /

**Expected result:** :math:`\sigma_{DC} \sim 10^7` S/m for Cu (literature: :math:`5.96 \times 10^7` S/m)

Performance Considerations
==========================

Computational Cost
------------------

**Scaling:**

.. math::

   \text{Cost} \sim \mathcal{O}(N_E \times N_{ll}^2 \times N_{orb}^2 \times N_r)

where:

- :math:`N_E`: Energy points (200-500)
- :math:`N_{ll}`: Recursion depth (50-200)
- :math:`N_{orb}`: Orbital basis (:math:`\sim 10-20`)
- :math:`N_r`: Random vectors (10-50)

**Typical timings:**

- Small system (10 atoms): ~10-30 minutes
- Large system (100 atoms): ~1-3 hours (with stochastic)

Memory Requirements
-------------------

**Main arrays:**

- ``gamma_nm(N_E, N_ll, N_ll)`` (complex, double precision, :math:`N_E` = ``channels_ldos`` + 10): ~330 MB for ``channels_ldos`` = 500, :math:`N_{ll}=200`
- ``mu_nm(18, 18, N_ll, N_ll, N_{type})``: ~100 MB

**Tips for large systems:**

1. Use stochastic method (``random_vec_num = 20-50``)
2. Reduce ``cond_ll`` to 100-150
3. Use ``cond_calctype = 'random_vec'``

Provenance
==========

**Implementation:**

- ``source/conductivity.f90`` - Main conductivity module
- ``source/recursion.f90`` - Chebyshev recursion and moments
- ``source/green.f90`` - Green's function and velocity operator

**Key procedures:**

- ``calculate_gamma_nm()`` - Chebyshev expansion coefficients
- ``calculate_conductivity_tensor()`` - Main conductivity calculation

**Dependencies:**

- ``energy_mod`` - Energy mesh
- ``recursion_mod`` - Stochastic Chebyshev moments
- ``hamiltonian_mod`` - Hamiltonian and velocity operator
- ``control_mod`` - Calculation parameters

**Authors:**

- Angela Klautau, Ramon Cardias, Lucas P. Campagna
- S. Frota-Pessôa, Pascoal R. Peduto
- Anders Bergman, Ivan P. Miranda
- S. B. Legoas, H. M. Petrilli

References
==========

1. **Kubo-Bastin formula:** A. Bastin, C. Lewiner, O. Betbeder-Matibet, and P. Nozieres,
   *J. Phys. Chem. Solids* **32**, 1811 (1971)

2. **Kernel Polynomial Method:** A. Weiße, G. Wellein, A. Alvermann, and H. Fehske,
   *Rev. Mod. Phys.* **78**, 275 (2006)

3. **Chebyshev expansion for transport:** J. H. García, L. Covaci, and T. G. Rappoport,
   *Phys. Rev. Lett.* **114**, 116602 (2015)

4. **Stochastic trace:** R. Iitaka and T. Ebisuzaki,
   *Phys. Rev. E* **69**, 057701 (2004)

5. **LMTO conductivity:** S. Frota-Pessôa and collaborators, implementation notes

See Also
========

- :ref:`theory/green_functions` - Green's function formalism
- :ref:`theory/recursion_method` - Chebyshev recursion details
- :ref:`user_guide/output_files` - Conductivity output files
- :ref:`keywords/output_options` - Conductivity parameters
