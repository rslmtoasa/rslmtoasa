.. _examples/exchange_bccfe:

=============================================
Post-Processing: Exchange Interactions (Fe)
=============================================

**Location:** ``example/exchange/bccFe/``

**System:** Magnetic exchange parameters J_ij for bcc Fe

**Physics:** Calculate Heisenberg exchange interactions using LKAG method.

Overview
========

After obtaining converged electronic structure, this example calculates magnetic exchange parameters J_ij that define the Heisenberg model:

.. math::

   H = -\sum_{i<j} J_{ij} \mathbf{S}_i \cdot \mathbf{S}_j

Using the Liechtenstein-Katsnelson-Antropov-Guslienko (LKAG) method based on magnetic force theorem.

Files
=====

- ``input.nml`` - Exchange calculation setup
- ``Fe.nml`` - Converged Fe parameters

Input File: input.nml
=====================

.. code-block:: fortran

   &calculation
   post_processing = 'exchange'   ! ← Post-processing mode
   verbose = T
   /
   
   &lattice
   ndim = 50000                   ! Large system for exchange
   rc = 80
   alat = 2.86120
   crystal_sym = 'bcc'
   wav = 1.40880
   
   njij = 2                       ! ← Number of exchange pairs
   ijpair(1, :) = 1, 2634         ! Pair 1: atoms 1 ↔ 2634
   ijpair(2, :) = 1, 2635         ! Pair 2: atoms 1 ↔ 2635
   /
   
   &atoms
   database = './'
   label(1) = 'Fe'
   /
   
   &self
   nstep = 100
   /
   
   &energy
   fermi = -0.069291              ! Use converged Fermi energy
   energy_min = -1.0
   energy_max = 1.2
   channels_ldos = 2500
   /
   
   &control
   calctype = 'B'
   nsp = 2
   lld = 21
   recur = 'block'
   /
   
   &mix
   beta = 0.5                     ! Can use larger mixing for post-processing
   mixtype = 'broyden'
   /

Critical Exchange Parameters
-----------------------------

- ``post_processing = 'exchange'``: Activates exchange calculation
- ``njij = 2``: Calculate 2 exchange pairs
- ``ijpair(i, :) = atom1, atom2``: Specify which atom pairs to calculate
- ``ndim = 50000``: Large system size to include many neighbors
- Use converged ``fermi`` from SCF calculation

**Atom indices:** The numbers 2634, 2635 correspond to specific neighbor shells. To find indices:

1. Check cluster generation output
2. Identify atoms at desired distances
3. Common shells for bcc: 1NN at a√3/2, 2NN at a, 3NN at a√11/2

No lattice.nml Needed
=====================

Structure is same as bulk bccFe SCF calculation.

Running the Calculation
=======================

**Prerequisites:** Must have converged ``Fe.nml`` from SCF calculation.

.. code-block:: bash

   cd example/exchange/bccFe
   ../../../../build/bin/rslmto.x

Expected Output
===============

**Output files:**

- ``exchange_jij.out`` - Exchange parameters in mRy and meV
- ``exchange_full.out`` - Full exchange tensor (isotropic J, DMI, anisotropy)
- ``Fe_out.nml`` - Updated parameters

**Expected results for bcc Fe:**

.. code-block:: text

   Neighbor Shell    Distance (Å)    J_ij (mRy)    J_ij (meV)
   1NN (2634)           2.48          +20 to +24    +0.45 to +0.55
   2NN (2635)           2.86          +2 to +5      +0.05 to +0.11
   3NN                  3.31          -1 to +1      ~0

**Sign convention:** Positive J_ij = ferromagnetic coupling.

View results:

.. code-block:: bash

   cat exchange_jij.out
   grep "J(" exchange_full.out | head -5

**Expected runtime:** 5-15 minutes.

Physical Interpretation
========================

**Exchange parameters determine:**

- **Curie temperature:** T_C ∝ Σ_j J_ij
- **Magnon spectrum:** ω(q) depends on J(R)
- **Domain wall properties:** Width ∝ √(A/K) where A ∝ J

**For bcc Fe:**

- Strong ferromagnetic 1NN interaction → high T_C (1043 K)
- Weaker 2NN coupling
- J_ij decays with distance

Using Exchange Parameters
=========================

Export J_ij for downstream simulations:

**Monte Carlo** (Curie temperature):

.. code-block:: text

   # Use J_ij in Metropolis or Wang-Landau
   H = -sum(J_ij * S_i dot S_j for all pairs)

**Atomistic spin dynamics:**

.. code-block:: text

   # Landau-Lifshitz-Gilbert equation
   dS/dt = -gamma S cross H_eff + alpha S cross dS/dt
   # where H_eff includes exchange field from J_ij

**Micromagnetic exchange stiffness:**

.. math::

   A = \\frac{J \\cdot S^2}{2a}

For Fe with J_1 ≈ 22 mRy, S = 1.1, a = 2.86 Å:

.. math::

   A \\approx 2 \\times 10^{-11} \\text{ J/m (agrees with experiment!)}

Troubleshooting
===============

**Exchange values seem wrong:**

- Verify converged SCF (check Fe.nml has reasonable magnetic moment)
- Ensure ``fermi`` energy is from converged calculation
- Check atom pair indices are correct

**Calculation very slow:**

- Reduce ``ndim`` (but ensure pairs are included)
- Reduce ``channels_ldos``

See Also
========

- :doc:`bulk_bccfe` - SCF calculation for this system
- :doc:`../../reference/exchange_module` - Exchange calculation details
- :doc:`../../theory/green_functions` - LKAG method theory
- :doc:`../../keywords/lattice_geometry` - njij, ijpair parameters
