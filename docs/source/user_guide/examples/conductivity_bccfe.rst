.. _examples/conductivity_bccfe:

===================================================
Post-Processing: Anomalous Hall Effect in Fe
===================================================

**Location:** ``example/conductivity/bccFe/``

**System:** Anomalous Hall conductivity of ferromagnetic bcc Fe

**Physics:** Calculate anomalous Hall effect arising from Berry curvature in ferromagnet.

Overview
========

In ferromagnets, a transverse (Hall) conductivity appears even without external magnetic field. This **anomalous Hall effect** (AHE) originates from:

- Berry curvature in electronic bands
- Spin-orbit coupling
- Broken time-reversal symmetry

Files
=====

- ``input.nml`` - AHE conductivity setup
- ``Fe.nml`` - Magnetic Fe parameters

Input File: input.nml
=====================

The SCF input in ``example/conductivity/bccFe/input.nml``:

.. code-block:: fortran

   &calculation
   pre_processing = 'bravais'
   verbose = T
   /
   &lattice
   rc = 080
   alat = 2.86120
   crystal_sym = 'bcc'
   wav = 1.40880
   pbc = .true.
   b1 = .false.
   b2 = .false.
   b3 = .false.
   n1 = 20
   n2 = 20
   n3 = 20
   ntype = 1
   ct(1) = 3.0d0
   r2 = 9.0d0
   /
   &atoms
   database = './'
   label(1) = 'Fe'
   /
   &self
   nstep = 100
   /
   &hamiltonian
   v_alpha = 0, 1, 0              ! Transverse current (y)
   v_beta = 1, 0, 0               ! Applied field (x)
   pol_alpha = 'z'                ! Magnetization (z)
   pol_beta = 'z'
   /
   &energy
   fermi = -0.069099
   energy_min = -1.0
   energy_max = 1.2
   channels_ldos = 2500
   /
   &control
   calctype = 'B'
   nsp = 2                        ! Collinear, with spin-orbit coupling
   cond_ll = 300
   lld = 150
   recur = 'chebyshev'
   linear_in = 'charge'
   linear_out = 'spin'            ! Spin current response
   cond_calctype = 'per_type'
   /
   &mix
   beta = 0.1
   mixtype = 'broyden'
   /

The conductivity step uses the folders ``sigma_xy_z-spin`` (``linear_out = 'spin'``) and
``sigma_xy_z-orb`` (``linear_out = 'orbital'``). Their ``input.nml`` equals the file above, except that
``pre_processing = 'bravais'`` is replaced by ``post_processing = 'conductivity'``.

Anomalous Hall Effect
=====================

The Hall resistivity in ferromagnets has two contributions:

.. math::

   \\rho_H = R_0 B + R_s M

where:

- R_0 B: Ordinary Hall effect (linear in magnetic field B)
- R_s M: **Anomalous Hall effect** (linear in magnetization M)

The anomalous Hall conductivity σ_xy^AHE is calculated from transverse velocity correlation.

Running the Calculation
=======================

The calculation has two steps. With ``post_processing = 'none'`` (the default) the file above runs
the SCF only:

.. code-block:: bash

   mkdir -p runs
   cp -r example/conductivity/bccFe runs/bccFe
   cd runs/bccFe
   ../../build/bin/rslmto.x

The conductivity step does not repeat the SCF. It reads the converged ``Fe_out.nml`` (the code reads
``<label>_out.nml`` before ``<label>.nml``), so copy it into the subfolder and run there:

.. code-block:: bash

   cp Fe_out.nml sigma_xy_z-spin/
   cd sigma_xy_z-spin
   ../../../build/bin/rslmto.x

Expected Output
===============

**Anomalous Hall conductivity:**

.. code-block:: text

   σ_xy^AHE ≈ 700-900 (Ω·cm)⁻¹
   
   Anomalous Hall angle:
   θ_AH = σ_xy / σ_xx ≈ 0.8-1.2%

**Spin-resolved conductivity:**

.. code-block:: bash

   cat cond_total.out       # E-E_F (Ry), Re, Im of the selected component (raw kernel)
   cat response_total.out   # same response in physical units, E-E_F in eV
   cat Fe_cond.out          # per-type contribution

Orbital-resolved parts are in ``cond_total_orb_real.out`` and ``cond_total_orb_im.out``.

**Expected runtime:** 15-25 minutes.

Physical Interpretation
========================

**Origin of AHE:**

1. **Intrinsic contribution:** Berry curvature in Bloch bands
   
   .. math::
   
      \Omega_n(\mathbf{k}) = \nabla_k \times \langle u_n | i\nabla_k | u_n \rangle

2. **Side-jump contribution:** Asymmetric scattering
3. **Skew-scattering contribution:** Spin-orbit coupling in scattering

The calculated σ_xy^AHE includes intrinsic + side-jump contributions.

**Applications:**

- **Magnetic sensors** (detect magnetization via Hall voltage)
- **Magnetic memory readout**
- **Topological spintronics**

Comparison with Ordinary Hall Effect
=====================================

**Ordinary Hall (non-magnetic Cu):**

- σ_xy ∝ B (external field)
- Small: ~10² (Ω·cm)⁻¹ at 1 Tesla

**Anomalous Hall (ferromagnetic Fe):**

- σ_xy ∝ M (magnetization)
- Large: ~10³ (Ω·cm)⁻¹ (no external field!)
- 10× larger signal

See Also
========

- :doc:`bulk_bccfe` - SCF calculation for Fe
- :doc:`exchange_bccfe` - Magnetic interactions
- :doc:`conductivity_fccpt` - Spin Hall effect
- :doc:`bulk_mn3sn` - Topological AHE in antiferromagnet
- :doc:`../../reference/conductivity_module` - Conductivity theory
