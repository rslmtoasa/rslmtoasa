.. _examples/conductivity_fcccu:

===================================================
Post-Processing: Conductivity of Cu
===================================================

**Location:** ``example/conductivity/fccCu/``

**System:** Electrical conductivity tensor for fcc Cu

**Physics:** Calculate conductivity using Kubo-Bastin formula with Chebyshev expansion (KPM).

Overview
========

This example calculates the electrical conductivity tensor σ_αβ(E) as a function of energy using:

- Kubo-Bastin formula for disordered systems
- Chebyshev polynomial expansion (Kernel Polynomial Method)
- Velocity operator matrix elements

Files
=====

- ``input.nml`` - Conductivity calculation setup
- ``Cu.nml`` - Cu tight-binding parameters

Input File: input.nml
=====================

.. code-block:: fortran

   &calculation
   pre_processing = 'bravais'
   verbose = T
   /
   
   &lattice
   rc = 60
   alat = 3.61411                 ! FCC Cu lattice constant
   crystal_sym = 'fcc'
   wav = 1.41237
   pbc = .true.
   b1 = .false.                   ! ← Open boundaries
   b2 = .false.                   ! (breaks translational symmetry)
   b3 = .false.
   n1 = 20
   n2 = 20
   n3 = 20
   ntype = 1
   ct(1) = 4.0
   r2 = 16.0
   /
   
   &atoms
   database = './'
   label(1) = 'Cu'
   /
   
   &self
   ws_all = .true.
   nstep = 100
   /
   
   &hamiltonian
   v_alpha = 0, 1, 0              ! ← Current direction (y)
   v_beta = 1, 0, 0               ! ← Field direction (x)
   js_alpha = 'z'                 ! ← Spin current (z)
   /
   
   &energy
   fermi = -0.088977
   energy_min = -1.0
   energy_max = 1.2
   channels_ldos = 2500
   /
   
   &control
   calctype = 'B'
   nsp = 2
   cond_ll = 500                  ! ← Chebyshev expansion order
   recur = 'chebyshev'            ! ← Use Chebyshev recursion
   lld = 150
   cond_type = 'spin'             ! ← Spin contribution
   cond_calctype = 'per_type'     ! ← Per-atom-type resolution
   /
   
   &mix
   mixtype = 'broyden'
   beta = 0.4
   /

Critical Conductivity Parameters
---------------------------------

- ``cond_ll = 500``: Chebyshev expansion order (higher → better resolution, but slower)
- ``cond_type = 'spin'`` or ``'orb'``: Calculate spin or orbital contribution
- ``cond_calctype = 'per_type'``: Atom-type resolved conductivity
- ``v_alpha, v_beta``: Velocity operator directions for tensor component σ_αβ
- ``b1 = b2 = b3 = .false.``: **Open boundaries required** for conductivity
- ``recur = 'chebyshev'``: Use Chebyshev recursion (KPM method)

Velocity Operator: Selecting Tensor Component
==============================================

The ``&hamiltonian`` namelist specifies which conductivity tensor element to calculate:

.. code-block:: text

   v_alpha = 0, 1, 0    ! y-direction
   v_beta = 1, 0, 0     ! x-direction
   => Calculates sigma_yx

**Common configurations:**

- **Longitudinal:** ``v_alpha = (1,0,0), v_beta = (1,0,0)`` → σ_xx
- **Hall:** ``v_alpha = (0,1,0), v_beta = (0,0,1)`` → σ_yz
- **Cross terms:** Any combination

For cubic Cu, by symmetry: σ_xx = σ_yy = σ_zz and off-diagonal terms are zero.

Running the Calculation
=======================

.. code-block:: bash

   cd example/conductivity/fccCu
   ../../../../build/bin/rslmto.x

Expected Output
===============

**Output files:**

- ``cond_total.out`` - Total conductivity vs energy
- ``cond_Cu.out`` - Cu contribution (if per_type)
- ``cond_spin.out`` or ``cond_orb.out`` - Spin/orbital resolved

**Expected results:**

For fcc Cu at Fermi level:

.. code-block:: text

   σ_xx(E_F) ≈ 6.0 × 10⁷ (Ω·m)⁻¹
   
   Experimental: 5.96 × 10⁷ (Ω·m)⁻¹
   Agreement: ~1%!

**View results:**

.. code-block:: bash

   grep -A 5 "Fermi" cond_total.out
   
   # Plot conductivity vs energy
   gnuplot -e "plot 'cond_total.out' u 1:2 w l; pause -1"

**Output format:**

.. code-block:: text

   # Energy (Ry)    σ (10⁷ Ω⁻¹m⁻¹)    Re[σ]    Im[σ]
   -0.100          2.15             2.15      0.01
   -0.090          4.82             4.82      0.02
   -0.089 (E_F)    6.03             6.03      0.00
   -0.088          5.71             5.71     -0.01

**Expected runtime:** 10-20 minutes.

Physical Interpretation
========================

**Conductivity peak at Fermi level:**

- High DOS at E_F (d-band contribution)
- Long scattering time (weak disorder)
- Cu is excellent conductor

**Energy dependence:**

- Below E_F: Decreasing σ (fewer available states)
- Above E_F: Also decreasing (fewer occupied states)
- Maximum at E_F for metals

Troubleshooting
===============

**Poor energy resolution:**

- Increase ``cond_ll`` to 1000 or higher
- Increase ``lld`` to 200-300

**Wrong conductivity value:**

- Check Fermi energy is correct
- Verify velocity operators are normalized
- Ensure b1=b2=b3=.false. (open boundaries)

See Also
========

- :doc:`conductivity_fccpt` - Spin Hall effect in Pt
- :doc:`conductivity_bccfe` - Anomalous Hall effect in Fe
- :doc:`../../reference/conductivity_module` - Conductivity calculation details
- :doc:`../../keywords/control_parameters` - cond_ll, cond_type parameters
