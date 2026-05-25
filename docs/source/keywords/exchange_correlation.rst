.. _keywords/exchange_correlation:

============================================
Exchange-Correlation Functional (&control)
============================================

Overview
========

The ``txc`` parameter in ``&control`` selects the exchange-correlation (XC) functional
used to describe electron-electron interactions.

txc Parameter
=============

**Type:** Integer

**Purpose:** Selection of XC functional

**Allowed values:**

.. list-table::
   :widths: 10 25 15 15
   :header-rows: 1

   * - txc
     - Functional
     - Type
     - Reference
   * - 1
     - Barth-Hedin (BH)
     - LDA
     - [1]
   * - 2
     - Slater X-α
     - LDA
     - [2]
   * - 3
     - Barth-Hedin-Janak
     - LDA
     - [3]
   * - 4
     - Vosko-Wilk-Nusair
     - LDA
     - [4]
   * - 5
     - Perdew-Zunger
     - LDA
     - [5]
   * - 8
     - Perdew-Burke-Enzerhof (PBE)
     - GGA
     - [6]

**Default:** 1 (Barth-Hedin)

**Example:**

.. code-block:: fortran

   txc = 1  ! Barth-Hedin LDA (standard)

Physical Meaning
================

**LDA (Local Density Approximation):**

XC energy depends only on charge density at each point:

.. math::

   E_{XC}^{\text{LDA}}[n] = \int d^3r \, \epsilon_{XC}(n(\mathbf{r})) n(\mathbf{r})

- Simple and fast
- Often excellent agreement with experiment
- Universal: works for many systems

**GGA (Generalized Gradient Approximation):**

XC energy also depends on density gradient:

.. math::

   E_{XC}^{\text{GGA}}[n] = \int d^3r \, \epsilon_{XC}(n(\mathbf{r}), \nabla n(\mathbf{r})) n(\mathbf{r})

- More sophisticated than LDA
- Better for geometries, surfaces, and magnetic moments
- Slightly slower than LDA

Functional Comparison
=====================

**Barth-Hedin (txc=1):**

- **Pros:** Standard for LMTO, well-tested with TB parameters
- **Cons:** Slight underestimate of band gaps
- **Best for:** Metals, general-purpose calculations
- **References:** J. Phys. C 5, 1629 (1972)

**Slater X-α (txc=2):**

- **Pros:** Simplest LDA functional
- **Cons:** Less accurate than Barth-Hedin
- **Best for:** Quick exploratory calculations
- **References:** Phys. Rev. B 12, 3060 (1975)

**Vosko-Wilk-Nusair (txc=4):**

- **Pros:** Accurate for materials science
- **Cons:** More complicated expression
- **Best for:** Semiconductors, transition metals
- **References:** Can. J. Phys. 58, 1200 (1980)

**PBE (txc=8):**

- **Pros:** Better band gaps, magnetization
- **Cons:** Requires revised TB parameters
- **Best for:** Modern high-accuracy calculations
- **References:** Phys. Rev. Lett. 77, 3865 (1996)

Implementation Details
======================

The XC functional is implemented in ``source/xc.f90`` with functions:

.. code-block:: fortran

   type :: xc
      procedure :: XCPOT       ! Main XC potential routine
      procedure :: EXCHPBE     ! Exchange part (PBE)
      procedure :: CORPBE      ! Correlation part (PBE)
      procedure :: LAGGGA      ! GGA functional (generalized)
      ! ... LDA functionals
   end type xc

**During SCF:**

1. Charge density integrated from Green's function
2. XC potential computed: $V_{XC} = \delta E_{XC} / \delta n$
3. Added to total Hamiltonian
4. New Green's functions computed
5. Iterate

Potential Parameters and txc
============================

**Important relationship:**

The tight-binding potential parameters (``center_band``, ``width_band``, etc.)
are typically **fitted to specific XC functional**.

**Example:**

Parameters in ``database/Si1.nml`` may be fitted assuming:

.. code-block:: fortran

   txc = 1  ! Barth-Hedin

**If you change txc:**

- Results may become inaccurate
- Band structure shapes may change
- Magnetic moments may shift
- Consider refitting parameters (advanced topic)

**Recommendation:**

- Use ``txc`` consistent with database potentials
- Document which txc value was used for parameter fitting
- Verify band structure against reference calculations

Temperature Dependence
======================

LDA and GGA XC potentials are **T=0 functionals**.

For **finite temperature** calculations:

- XC functional doesn't change formally
- Fermi-Dirac distribution broadens electronic levels
- Effective electronic temperature matters for stability

Parameter: ``temperature`` in ``&energy`` namelist.

Magnetic Moment Predictions
===========================

For **magnetic systems** (nsp > 1):

- LDA typically **underestimates** magnetic moments by 5-15%
- GGA (PBE) often **closer** to experimental values
- Spin-orbit coupling (nsp=2,4) reduces moments

**Example:** Fe

- Experimental: m ≈ 2.22 μ_B
- LDA (Barth-Hedin): m ≈ 2.0 μ_B
- GGA (PBE): m ≈ 2.1 μ_B
- LDA + SOC: m ≈ 1.9 μ_B

Meta-GGA and Beyond
====================

**Not currently in RS-LMTO:**

- MGGA functionals (TPSS, SCAN, etc.)
- Hybrid functionals (PBE0, HSE, etc.)
- Range-separated functionals

**If needed:**

- Fit TB parameters to ab-initio GGA/MGGA calculations
- Use those parameters with RS-LMTO LDA (approximate)
- Or contribute new XC functional to code

Choosing the Right Functional
=============================

**For most calculations:** txc = 1 (Barth-Hedin)

- Standard choice
- Well-tested parameters available
- Good balance of accuracy and speed

**For semiconductors/insulators:** txc = 8 (PBE)

- Better band gaps
- More accurate structural properties

**For magnetic systems:** 

- LDA (txc=1): works reasonably
- GGA (txc=8): slightly better moments
- Experiment determines best choice for your system

**For exploratory work:** txc = 2 (Slater X-α)

- Simplest and fastest
- Adequate accuracy for survey calculations

**For publication-quality:** 

- Match functional to database parameters
- Verify against experimental data
- Consider both LDA and GGA results

Provenance
==========

XC functionals implemented in:

- **Type definition:** ``source/xc.f90::type xc``
- **Initialization:** ``source/xc.f90::constructor()``
- **Potential calculation:** ``source/xc.f90::XCPOT()``
- **LDA functionals:** ``source/xc.f90`` (Barth-Hedin, Slater, etc.)
- **GGA functionals:** ``source/xc.f90`` (PBE implementation)

Parameter fitting typically involves:

- ``source/potential.f90`` - Holds fitted parameters
- Database files: ``element_name.nml`` - Pre-fitted for specific txc

See Also
========

- :ref:`keywords/control_parameters` - Related: nsp, verbose
- :ref:`theory/lmto_asa_overview` - XC functional theory
- :doc:`../user_guide/examples` - Example txc values
- ``source/xc.f90`` - Implementation source code

References
==========

[1] Barth, U.v., & Hedin, L. (1972). J. Phys. C: Solid State Phys., 5(13), 1629.
[2] Slater, J.C. (1951). Physical Review, 81(3), 385.
[3] Janak, J.F. (1978). Physical Review B, 18(12), 7165.
[4] Vosko, S.H., Wilk, L., & Nusair, M. (1980). Can. J. Phys., 58(8), 1200.
[5] Perdew, J.P., & Zunger, A. (1981). Physical Review B, 23(10), 5048.
[6] Perdew, J.P., Burke, K., & Ernzerhof, M. (1996). Phys. Rev. Lett., 77(18), 3865.
