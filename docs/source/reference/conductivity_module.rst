.. _reference/conductivity_module:

=======================================
Conductivity Module
=======================================

Overview
========

The **conductivity module** (``conductivity.f90``) calculates the **electrical conductivity tensor** $\sigma_{\alpha\beta}(\omega)$ using the **Kubo-Bastin formalism** combined with **Chebyshev polynomial expansion** within the LMTO-ASA framework.

**Key features:**

- Optical conductivity $\sigma(\omega)$ as function of frequency
- DC conductivity $\sigma(0)$
- Full $3 \times 3$ conductivity tensor (anisotropic transport)
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
         ! Dimension: (energy_points, ll_cond, ll_cond)
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

- $\Omega$ is the system volume
- $\hat{j}_\alpha$ is the current operator in direction $\alpha$
- $\langle \ldots \rangle$ denotes quantum/thermal average
- $\eta \to 0^+$ is a positive infinitesimal

Velocity Operator and Current
------------------------------

The **current operator** is related to the **velocity operator** $\hat{v}$:

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

- $G^+(E) = (E - \hat{H} + i\eta)^{-1}$ is the retarded Green's function
- $\mu$ is the chemical potential (Fermi energy)
- $\text{Im}[G^+] = A(E)/(2\pi)$ is the spectral function

For **DC conductivity** ($\omega = 0$):

.. math::

   \sigma_{\alpha\beta}(0) = \frac{e^2\hbar}{\pi\Omega} \int_{-\infty}^{\mu} dE \, 
   \text{Tr}\left[ \hat{v}_\alpha A(E) \hat{v}_\beta A(E) \right]

where $A(E) = -2\pi \text{Im}[G^+(E)]$ is the spectral function.

Chebyshev Polynomial Expansion
===============================

Kernel Polynomial Method (KPM)
-------------------------------

To avoid expensive direct evaluation, the **Chebyshev expansion** approximates the Green's function:

**Energy rescaling:** First, map the Hamiltonian spectrum $[E_{\min}, E_{\max}]$ to $[-1, 1]$:

.. math::

   \tilde{H} = \frac{H - b}{a}, \quad a = \frac{E_{\max} - E_{\min}}{2 - \epsilon}, \quad b = \frac{E_{\max} + E_{\min}}{2}

where $\epsilon \sim 0.3$ provides margin.

**Chebyshev expansion** of Green's function:

.. math::

   G^+(\tilde{E}) = \sum_{n=0}^{N-1} \mu_n T_n(\tilde{E}) \frac{g_n}{\sqrt{1 - \tilde{E}^2}}

where:

- $T_n(x)$ are Chebyshev polynomials of the first kind
- $\mu_n = 2 - \delta_{n,0}$ are normalization factors
- $g_n$ are Jackson kernel damping factors (smoothing)
- $N$ is the expansion order (recursion depth)

**Moments** $\mu_n$ are calculated via recursion:

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

The conductivity involves **products of two Green's functions**, leading to the $\Gamma_{nm}$ function:

.. math::

   \Gamma_{nm}(E) = \frac{1}{(1 - \tilde{E}^2)^2} 
   \left[ C_n(E) T_m(\tilde{E}) + C_m(E) T_n(\tilde{E}) \right] g_n g_m w_n w_m

where:

- $C_n(E) = \frac{1}{\pi} \frac{(-1)^{n+1} \sin(n \arccos \tilde{E})}{\sqrt{1 - \tilde{E}^2}}$ is a Chebyshev derivative
- $w_n$ are integration weights
- $g_n$ are Jackson kernel factors

**Jackson kernel** (smoothing to reduce Gibbs oscillations):

.. math::

   g_n = \frac{1}{N+1} \left[ (N - n + 1)\cos\frac{\pi n}{N+1} + \sin\frac{\pi n}{N+1} \cot\frac{\pi}{N+1} \right]

**Implementation** in ``calculate_gamma_nm()``:

.. code-block:: text

   ! Energy rescaling
   a = (e_max - e_min) / (2 - 0.3)
   b = (e_max + e_min) / 2
   wscale(:) = (en%ene(:) - b) / a
   
   ! Chebyshev polynomials
   T_0 = 1.0
   T_1 = wscale
   do n = 2, ll_cond
      T_n = 2 * wscale * T_{n-1} - T_{n-2}
   end do
   
   ! Cn coefficients
   acos_x = acos(wscale)
   sqrt_term = sqrt(1.0 - wscale**2)
   cn(n) = ((-1)**(n+1)) * sin(n * acos_x) / (pi * sqrt_term)
   
   ! Gamma_nm matrix
   do n = 1, ll_cond
      do m = 1, ll_cond
         gamma_nm(:,n,m) = (cn(:,n) * T_m + cn(:,m) * T_n) / (1 - wscale**2)**2
         gamma_nm(:,n,m) = gamma_nm(:,n,m) * g_n * g_m * w_n * w_m
      end do
   end do

Stochastic Trace Evaluation
----------------------------

For large systems, exact trace evaluation $\text{Tr}[\ldots]$ is prohibitive. Instead, use **stochastic trace estimation**:

.. math::

   \text{Tr}[A] \approx \frac{1}{N_r} \sum_{r=1}^{N_r} \langle \xi_r | A | \xi_r \rangle

where $|\xi_r\rangle$ are **random vectors** with elements drawn from:

- Complex normal distribution: $\xi_{i,r} \sim \mathcal{N}(0, 1) + i\mathcal{N}(0, 1)$
- Or: Random phase: $\xi_{i,r} = e^{i\phi_r}$

**Variance** decreases as $\mathcal{O}(1/\sqrt{N_r})$, so typically $N_r = 10-50$ random vectors suffice.

Conductivity Tensor Calculation
================================

Full Tensor Formula
-------------------

The **conductivity tensor** $\sigma_{\alpha\beta}(E)$ has 9 components (for 3D):

.. math::

   \sigma_{\alpha\beta}(E) = \frac{16 e^2}{\pi \Omega \Delta E^2} \sum_{n,m}^{N} 
   \text{Tr}\left[ \hat{v}_\alpha \Gamma_{nm}(E) \hat{v}_\beta \mu_{nm} \right]

where:

- $\mu_{nm}$ are the Chebyshev moments of the velocity operator products
- $\Delta E = E_{\max} - E_{\min}$ is the energy span
- Factor $16/\pi$ comes from Chebyshev normalization

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

giving $s$, $p$, $d$ contributions.

2. **By atomic type:**

.. math::

   \sigma_{\alpha\beta}^{(\text{type})}(E) = \sum_{i \in \text{type}} \langle i | \ldots | i \rangle

3. **By energy range:** Integrate from $E_{\min}$ to variable upper limit $E$

Implementation Details
======================

Algorithm Overview
------------------

**Main procedure:** ``calculate_conductivity_tensor()``

1. **Precompute $\Gamma_{nm}$:** Call ``calculate_gamma_nm()``
2. **Stochastic recursion:** Compute Chebyshev moments $\mu_{nm}^{\alpha\beta}$ via recursion module
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

**Option 1:** Conductivity averaged over all atoms (``cond_calctype = 'average'``)

**Option 2:** Conductivity per atomic species (``cond_calctype = 'per_type'``):

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

- ``cond_ll`` (default: 100): Chebyshev expansion order $N$
  
  - Larger $N$ -> better energy resolution, smoother curves
  - Typical: 50-200
  - Convergence: check $\sigma(E)$ vs $N$

- ``cond_calctype`` (default: 'average'): Calculation mode
  
  - ``'average'``: Total conductivity
  - ``'per_type'``: Per-atomic-species conductivity
  - ``'random_vec'``: With random vectors (stochastic)

- ``random_vec_num`` (default: 10): Number of stochastic vectors
  
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

   - Default: 200-500 points
   - Must resolve features near $E_F$

3. **Random vectors** ``random_vec_num``:

.. code-block:: bash

   # Check variance in repeated runs
   σ_avg +- σ_std  # Should be < 5% for good statistics

Main Procedures
===============

calculate_gamma_nm
------------------

**Purpose:** Precompute the $\Gamma_{nm}(E)$ coefficients for Chebyshev expansion.

**Inputs:**

- Energy mesh ``en%ene(:)``
- Recursion depth ``ll_cond``
- Energy bounds ``e_min, e_max``

**Algorithm:**

1. Rescale energies: $\tilde{E} = (E - b)/a$
2. Compute Chebyshev polynomials $T_n(\tilde{E})$ recursively
3. Compute $C_n(E)$ coefficients
4. Apply Jackson kernel $g_n$
5. Construct $\Gamma_{nm}(E)$ matrix

**Output:**

- ``gamma_nm(energy, n, m)`` - Complex array of Chebyshev coefficients

calculate_conductivity_tensor
------------------------------

**Purpose:** Compute full $3 \times 3$ conductivity tensor $\sigma_{\alpha\beta}(E)$.

**Algorithm:**

1. Loop over energy points $E$
2. For each $E$:

   a. Sum over Chebyshev indices: $\sum_{nm} \Gamma_{nm}(E) \mu_{nm}^{\alpha\beta}$
   b. Multiply by prefactor: $(16 e^2)/(\pi \Omega \Delta E^2)$
   c. Simpson integrate from $E_{\min}$ to $E$: $\sigma(E) = \int^E$ 

3. Output:

   - Real part: Dissipative conductivity (absorption)
   - Imaginary part: Reactive conductivity (dispersion)

**Outputs:**

- Total: ``cond_total.out``
- Orbital-resolved: ``cond_total_orb_real.out``, ``cond_total_orb_im.out``
- Per-type: ``<element>_cond.out``

Output Files
============

cond_total.out
--------------

Total conductivity as function of energy:

.. code-block:: text

   # E-EF(eV)   Re[σ(E)] (arbitrary units)   Im[σ(E)]
     -5.000    125.43                         2.34
     -4.990    126.12                         2.45
     ...
      0.000    145.67                         5.21  ← Fermi level
     ...

**Columns:**

1. Energy relative to Fermi level $(E - E_F)$
2. Real part of conductivity (dissipation)
3. Imaginary part (dispersion)

cond_total_orb_real.out
-----------------------

Orbital-resolved real conductivity:

.. code-block:: text

   # E-EF  σ_s↑↑  σ_px↑↑  σ_py↑↑  σ_pz↑↑  σ_d_z2↑↑  ...  σ_d_xy↓↓
     -5.0  5.2    12.3    12.3    12.3    18.4      ...  15.2
     ...

**18 components:**

- Spin-up: $s, p_x, p_y, p_z, d_{z^2}, d_{xz}, d_{yz}, d_{x^2-y^2}, d_{xy}$ (9 orbitals)
- Spin-down: Same 9 orbitals

<element>_cond.out
------------------

Per-atomic-species conductivity (if ``cond_calctype = 'per_type'``):

**Files:** ``Fe_cond.out``, ``Pt_cond.out``, etc.

**Format:** Same as ``cond_total.out`` but for specific element

Physical Interpretation
=======================

DC Conductivity
---------------

**At Fermi level** $E = E_F$, $\sigma(E_F)$ gives **DC conductivity**:

.. math::

   \sigma_{DC} = \sigma(E_F, \omega=0)

**Units:** Typically in $(e^2/\hbar) \cdot a_0^{-1}$ or converted to SI.

**Conversion:** To S/m (Siemens per meter):

.. math::

   \sigma[\text{S/m}] = \sigma[\text{code units}] \times \frac{e^2}{\hbar a_0 \Omega}

where $\Omega$ is the unit cell volume.

Optical Conductivity
--------------------

**For finite frequency** $\omega \neq 0$, the **optical conductivity** describes:

- Absorption spectrum
- Dielectric function via Kramers-Kronig
- Optical transitions

**Real part $\text{Re}[\sigma(\omega)]$:**

- Absorption/dissipation
- Peaks at interband transitions

**Imaginary part $\text{Im}[\sigma(\omega)]$:**

- Dispersion/reactive response
- Related to plasma frequency

Drude Model Comparison
----------------------

At low frequencies, metallic systems follow **Drude model**:

.. math::

   \sigma(\omega) = \frac{\sigma_0}{1 - i\omega\tau}

where $\tau$ is the scattering time. From LMTO:

.. math::

   \sigma_0 = \frac{ne^2\tau}{m^*}

can be extracted by fitting near $E_F$.

Integration and Usage
=====================

Input Namelist (&conductivity)
-------------------------------

The ``&conductivity`` namelist accepts a single parameter that controls the
accuracy of the Chebyshev expansion used inside the transport calculation.

**ll_cond**

**Type:** Integer

**Purpose:** Chebyshev polynomial expansion order (recursion depth) for the
conductivity tensor calculation.

**Typical range:** 100–600

**Default:** set by ``cond_ll`` in ``&control`` if not specified here

**Example:**

.. code-block:: fortran

   &conductivity
     ll_cond = 400
   /

**Notes:**

- Higher values produce a smoother, more converged conductivity spectrum but
  increase computational cost as :math:`\mathcal{O}(N_{ll}^2)`.
- A convergence test is recommended: increase ``ll_cond`` until
  :math:`\sigma(E_F)` changes by less than ~1%.
- ``ll_cond`` in this namelist takes precedence over ``cond_ll`` in
  ``&control``; either can be used.

**See also:** ``cond_ll``, ``cond_type``, ``cond_calctype`` in
:ref:`keywords/control_parameters`.

Input Parameters
----------------

**Minimal setup:**

.. code-block:: fortran

   &control
      calc_conductivity = .true.
      cond_ll = 100             ! Chebyshev expansion order
      cond_calctype = 'average' ! or 'per_type'
      random_vec_num = 20       ! Stochastic vectors
   /

**Advanced:**

.. code-block:: fortran

   &energy
      channels_ldos = 500       ! Energy mesh points
      energy_min = -10.0        ! eV below EF
      energy_max = +5.0         ! eV above EF
   /

Workflow
--------

1. **Complete SCF calculation** first
2. **Enable conductivity:** Set ``calc_conductivity = .true.``
3. **Run calculation:** Conductivity computed in post-processing
4. **Check convergence:** Vary ``cond_ll``, check $\sigma(E_F)$ stability
5. **Analyze output:** Plot $\sigma(E)$, identify features

Example: Conductivity of Cu
----------------------------

.. code-block:: fortran

   &calculation
      calctype = 'S'
   /
   &control
      calc_conductivity = .true.
      cond_ll = 150
      cond_calctype = 'average'
   /
   &energy
      channels_ldos = 400
      energy_min = -8.0
      energy_max = 3.0
   /

**Expected result:** $\sigma_{DC} \sim 10^7$ S/m for Cu (literature: $5.96 \times 10^7$ S/m)

Performance Considerations
==========================

Computational Cost
------------------

**Scaling:**

.. math::

   \text{Cost} \sim \mathcal{O}(N_E \times N_{ll}^2 \times N_{orb}^2 \times N_r)

where:

- $N_E$: Energy points (200-500)
- $N_{ll}$: Recursion depth (50-200)
- $N_{orb}$: Orbital basis ($\sim 10-20$)
- $N_r$: Random vectors (10-50)

**Typical timings:**

- Small system (10 atoms): ~10-30 minutes
- Large system (100 atoms): ~1-3 hours (with stochastic)

Memory Requirements
-------------------

**Main arrays:**

- ``gamma_nm(N_E, N_ll, N_ll)``: ~500 MB for $N_E=500$, $N_{ll}=200$
- ``mu_nm(18, 18, N_ll, N_ll, N_{type})``: ~100 MB

**Tips for large systems:**

1. Use stochastic method (``random_vec_num = 20-50``)
2. Reduce ``cond_ll`` to 100-150
3. Use ``cond_calctype = 'average'`` (not 'per_type')

Provenance
==========

**Implementation:**

- ``source/conductivity.f90`` - Main conductivity module (379 lines)
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
