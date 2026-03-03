.. _theory/constraining_fields:

================================================
Constraining Fields for Non-Collinear Magnetism
================================================

Introduction
============

In non-collinear spin density functional theory, the local magnetization on each
atomic site is a three-dimensional vector :math:`\mathbf{m}_i` whose direction and
magnitude are both determined self-consistently.  For many applications — exchange
interaction calculations, spin-spiral total energies, or constrained spin-dynamics
simulations — it is necessary to **fix (or steer) the magnetization direction**
on one or more sites toward a prescribed reference direction
:math:`\hat{\mathbf{e}}_i^{\,\mathrm{ref}}`, while allowing the magnitude to
relax freely.

This is achieved by augmenting the Kohn-Sham Hamiltonian with a site-dependent
**constraining field** :math:`\mathbf{B}_i^{\,\mathrm{con}}`:

.. math::

   \hat{H} = \hat{H}_{\mathrm{KS}} + \sum_i \mathbf{B}_i^{\,\mathrm{con}} \cdot \hat{\boldsymbol{\sigma}}_i,

where :math:`\hat{\boldsymbol{\sigma}}_i` is the vector of Pauli matrices projected
onto site :math:`i`.  The field :math:`\mathbf{B}_i^{\,\mathrm{con}}` is updated
iteratively so that the self-consistent magnetization direction converges to
:math:`\hat{\mathbf{e}}_i^{\,\mathrm{ref}}`.

RS-LMTO-ASA implements four distinct algorithms for determining
:math:`\mathbf{B}_i^{\,\mathrm{con}}`, selected by the input parameter
``constraints_i_cons``.  They are described in detail in the sections below.

.. note::

   The constraining field enters the Hamiltonian in Rydberg units internally.
   A unit-conversion factor :math:`b_{2T} = 235298.924` converts the internal
   field to Tesla when needed for output.

----

Method 2 — Lagrange Multiplier (no orthogonalization)
======================================================

**Physical idea**

The simplest approach adds a penalty energy proportional to the squared deviation
of the unit magnetization vector from the reference direction:

.. math::

   E_{\mathrm{con}} = \lambda_t \sum_i \left| \hat{\mathbf{m}}_i - \hat{\mathbf{e}}_i^{\,\mathrm{ref}} \right|^2,

where :math:`\hat{\mathbf{m}}_i = \mathbf{m}_i / |\mathbf{m}_i|` is the
normalised magnetization on site :math:`i`.

The constraining field on site :math:`i` is the functional derivative of
:math:`E_{\mathrm{con}}` with respect to :math:`\mathbf{m}_i`, pointing
opposite to the energy gradient to drive the system toward the minimum:

.. math::

   \mathbf{B}_i^{\,\mathrm{con}} = -\frac{\partial E_{\mathrm{con}}}{\partial \mathbf{m}_i}.

**Implementation**

In each SCF iteration the error vector and penalty energy for site :math:`i` are:

.. math::

   \boldsymbol{\delta}_i &= \hat{\mathbf{m}}_i - \hat{\mathbf{e}}_i^{\,\mathrm{ref}}, \\[4pt]
   E_{\mathrm{con}} &\mathrel{+}= \lambda_t \, |\boldsymbol{\delta}_i|^2.

Adaptive :math:`\lambda_t` update (end of iteration):

.. math::

   \lambda_t \leftarrow \min\!\left(1.2\,\lambda_t,\; 10^4\right)
   \quad \text{if } E_{\mathrm{con}} < 10^{-2},

ensuring the penalty strength grows whenever the constraint is well satisfied,
until it reaches a ceiling of :math:`10^4`.

.. note::

   For method 2 the constraining field vector is derived implicitly from
   :math:`E_{\mathrm{con}}` and fed into the Hamiltonian; the ``bfield`` array
   carries the accumulated field that enters the exchange-correlation potential
   on subsequent iterations.

**Code reference:** ``source/include_codes/abspinlib/constrain.f90``,
``i_cons == 2`` branch (lines 84–89).

----

Method 3 — Lagrange Multiplier with Orthogonalization (:math:`\mathbf{B} \perp \mathbf{m}`)
==============================================================================================

**Physical idea**

A known issue with the plain Lagrange approach is that the constraining field
may develop a component *parallel* to the local magnetization.  Such a
longitudinal component merely shifts the exchange splitting without changing the
moment direction, wasting iterations.

Method 3 applies a Gram-Schmidt projection to remove the component of
:math:`\boldsymbol{\delta}_i` parallel to the reference direction, enforcing
:math:`\mathbf{B}_i^{\,\mathrm{con}} \perp \hat{\mathbf{e}}_i^{\,\mathrm{ref}}` at each step.

**Implementation**

.. math::

   \boldsymbol{\delta}_i &= \hat{\mathbf{m}}_i - \hat{\mathbf{e}}_i^{\,\mathrm{ref}}, \\[4pt]
   \boldsymbol{\delta}_i^{\perp} &= \boldsymbol{\delta}_i
     - \left(\boldsymbol{\delta}_i \cdot \hat{\mathbf{e}}_i^{\,\mathrm{ref}}\right)
       \hat{\mathbf{e}}_i^{\,\mathrm{ref}},\\[4pt]
   E_{\mathrm{con}} &\mathrel{+}= \lambda_t \, |\boldsymbol{\delta}_i^{\perp}|^2.

The Gram-Schmidt step (second line) ensures that the penalty — and hence the
constraining field — lies in the plane perpendicular to :math:`\hat{\mathbf{e}}_i^{\,\mathrm{ref}}`.

Adaptive :math:`\lambda_t` update after the loop:

.. math::

   \lambda_t \leftarrow \min\!\left(
     \lambda_t \left(2 - \frac{\lambda_t}{\lambda_{\max}}\right),\;
     \lambda_{\max}
   \right).

This is a logistic-type update that accelerates growth when
:math:`\lambda_t \ll \lambda_{\max}` and saturates smoothly at the maximum
value :math:`\lambda_{\max} = 10`.

**Code reference:** ``source/include_codes/abspinlib/constrain.f90``,
``i_cons == 3`` branch (lines 90–96).

----

Method 4 — PID Controller (AMN Non-Collinear Constraints)
==========================================================

**Physical idea**

The Lagrange approaches update the field based only on the instantaneous error
(proportional control).  A **PID (Proportional-Integral-Derivative)** controller
uses three components to improve convergence speed and stability, especially for
slowly varying or oscillating moments:

- **P** (proportional): current error, provides restoring force.
- **I** (integral):   accumulated error history, eliminates steady-state offset.
- **D** (derivative): rate of change of error, provides damping.

**Implementation**

The error signal for site :math:`i` is the component of :math:`\mathbf{m}_i`
perpendicular to the reference direction :math:`\hat{\mathbf{e}}_i`:

.. math::

   \hat{\mathbf{e}}_i &= \frac{\hat{\mathbf{e}}_i^{\,\mathrm{ref}}}
                              {|\hat{\mathbf{e}}_i^{\,\mathrm{ref}}|}, \\[4pt]
   \mathbf{P}_i^{(n)} &= -\left[\mathbf{m}_i^{(n)}
     - \left(\mathbf{m}_i^{(n)} \cdot \hat{\mathbf{e}}_i\right)
       \hat{\mathbf{e}}_i\right].

The integral and derivative terms are accumulated/differenced across iterations:

.. math::

   \mathbf{I}_i^{(n)} &= \mathbf{I}_i^{(n-1)} + \mathbf{P}_i^{(n)}, \\[4pt]
   \mathbf{D}_i^{(n)} &= \mathbf{P}_i^{(n)} - \mathbf{P}_i^{(n-1)}.

The constraining field is then:

.. math::

   \mathbf{B}_i^{\,\mathrm{con}} = \lambda \left(
     1.30\, \mathbf{P}_i^{(n)}
     + 0.35\, \mathbf{I}_i^{(n)}
     - 0.10\, \mathbf{D}_i^{(n)}
   \right).

The PID gains :math:`(k_P, k_I, k_D) = (1.30, 0.35, -0.10)` were tuned
empirically for bcc Fe-type magnetic systems.

On the **first iteration** (when :math:`|\mathbf{P}_i^{(0)}| < 10^{-15}`) the
error is reduced by a factor of 10 to avoid an overshoot from the uninitialized
derivative term.

The constraining energy is accumulated Zeeman-like:

.. math::

   E_{\mathrm{con}} \mathrel{+}= \mathbf{B}_i^{\,\mathrm{con}} \cdot \mathbf{m}_i.

**Code reference:** ``source/include_codes/abspinlib/constrain.f90``,
``i_cons == 4`` branch (lines 97–155).

.. tip::

   Method 4 prints detailed per-atom PID diagnostics to standard output at every
   iteration, which is useful for understanding convergence behaviour.

----

Method 5 — Ma-Dudarev Approach
================================

**Physical idea**

The Ma-Dudarev scheme [#madudarev]_ provides a direct, physically motivated
update to the constraining field based on the instantaneous deviation of the
computed moments from the target, followed by Gram-Schmidt orthogonalization and
field mixing to control the step size.

**Implementation**

For each atom :math:`i` (skipped if :math:`|\mathbf{m}_i| < m_{\mathrm{thresh}} = 0.5`):

1. **Scale reference to match current magnitude:**

   .. math::

      \hat{\mathbf{e}}_i = \frac{\hat{\mathbf{e}}_i^{\,\mathrm{ref}}}
                                 {|\hat{\mathbf{e}}_i^{\,\mathrm{ref}}|}
                           \cdot |\mathbf{m}_i|.

2. **Compute raw constraining field:**

   .. math::

      \mathbf{B}_i^{\,\mathrm{new}} = \lambda_t \left(\mathbf{m}_i - \hat{\mathbf{e}}_i\right).

3. **Orthogonalize** (:math:`\mathbf{B} \perp \hat{\mathbf{e}}_i`):

   .. math::

      \mathbf{B}_i^{\,\mathrm{new}} \leftarrow \mathbf{B}_i^{\,\mathrm{new}}
        - \frac{\hat{\mathbf{e}}_i \cdot \mathbf{B}_i^{\,\mathrm{new}}}
               {|\mathbf{m}_i|^2}\, \hat{\mathbf{e}}_i.

4. **Mix with old field** (controls step size via :math:`\beta`):

   .. math::

      \mathbf{B}_i^{\,\mathrm{con}} = (1 - \beta)\,\mathbf{B}_i^{\,\mathrm{old}}
        - \beta\, \mathbf{B}_i^{\,\mathrm{new}},

   with default :math:`\beta = 1` (full replacement each step).

5. **Constraining energy** (Zeeman-like):

   .. math::

      E_{\mathrm{con}} \mathrel{+}= \lambda_t \left(
        |\mathbf{m}_i|
        - \frac{\mathbf{m}_i \cdot \hat{\mathbf{e}}_i}{|\mathbf{m}_i|}
      \right).

**Adaptive** :math:`\lambda_t` — two update rules are applied each iteration:

.. math::

   \lambda_t^{(n+1/2)} &= \min\!\left(
     \lambda_t^{(n)} + 1 - \frac{\lambda_t^{(n)}}{\lambda_{\max}},\;
     \lambda_{\max}
   \right), \\[4pt]
   \lambda_t^{(n+1)} &= \min\!\left(\lambda_t^{(n+1/2)} + 2,\; \lambda_{\max}\right).

The first rule is a logistic saturation; the second adds a constant increment,
providing rapid initial growth of :math:`\lambda_t` toward :math:`\lambda_{\max} = 10`.

**Code reference:** ``source/include_codes/abspinlib/constrain.f90``,
``i_cons == 5`` branch (lines 158–200).

.. [#madudarev] P.-W. Ma and S. L. Dudarev, *Phys. Rev. B* **91**, 054420 (2015).
   DOI: `10.1103/PhysRevB.91.054420 <https://doi.org/10.1103/PhysRevB.91.054420>`_

----

Comparison of Methods
=====================

.. list-table::
   :widths: 8 22 22 22 10
   :header-rows: 1

   * - ``i_cons``
     - Method
     - :math:`\mathbf{B}` derivation
     - Orthogonalization
     - Best for
   * - 2
     - Lagrange (plain)
     - Gradient of :math:`\lambda_t|\boldsymbol{\delta}|^2`
     - None
     - Simple collinear-like cases
   * - 3
     - Lagrange (:math:`\mathbf{B}\perp\mathbf{m}`)
     - Same + Gram-Schmidt
     - Yes (:math:`\mathbf{B}\perp\hat{\mathbf{e}}^{\,\mathrm{ref}}`)
     - General non-collinear
   * - 4
     - PID controller
     - Proportional + integral + derivative
     - Implicit (perpendicular error signal)
     - Oscillating or slowly converging moments
   * - 5
     - Ma-Dudarev
     - Direct moment deviation + orthogonalization
     - Yes (:math:`\mathbf{B}\perp\hat{\mathbf{e}}_i`)
     - Robust general-purpose non-collinear

----

Adaptive Lambda Schedule
=========================

All methods share a common infrastructure for adaptively growing the Lagrange
multiplier (penalty strength) :math:`\lambda_t` across SCF iterations, starting
from :math:`\lambda_t^{(0)} = \lambda_{\max}` at initialisation.  The update
rules applied at the end of each call to ``constrain()`` are:

- **Methods 2 & 3:** if :math:`E_{\mathrm{con}} < 10^{-2}`, multiply
  :math:`\lambda_t` by 1.2 (capped at :math:`10^4`); additionally for method 3
  apply the logistic saturation.
- **Method 4:** :math:`\lambda_t` is held at the fixed maximum :math:`\lambda`.
- **Method 5:** additive growth by 2 per iteration (capped at
  :math:`\lambda_{\max}`).

The motivation is to start with a moderate penalty (avoids over-steering at the
beginning of SCF) and tighten it as the calculation progresses.

----

Integration into the SCF and Exchange Loops
============================================

The constraining field is active whenever ``constraints_enable = .true.`` in the
input.  The call chain is:

.. code-block:: text

   constructor (self / exchange / conductivity)
       └─ initialize_cfd(nrec, 1, i_cons, code_prefac)
              Allocates PID arrays; sets lambda_t = lambda_max

   each SCF / exchange iteration
       └─ if (constraints_enable)
              mom_in  ← current computed moments
              mom_ref ← reference moments from input or initial guess
              bfield  ← 0 (re-initialised each call in self.f90)
              call constrain(mom_in, mom_ref, bfield, nrec)
              bfield  → added to effective exchange-correlation field
                        in the Hamiltonian

The ``bfield`` array (shape ``(3, natoms)``) is added to the exchange-correlation
potential before the next diagonalisation step, steering the moments toward
:math:`\hat{\mathbf{e}}^{\,\mathrm{ref}}`.

----

Provenance
==========

- **Core implementation:** ``source/include_codes/abspinlib/constrain.f90``
  — module ``cfd``, subroutines ``initialize_cfd`` and ``constrain``
- **Namelist definition:** ``source/include_codes/namelists/constraints.f90``
- **Control integration:** ``source/control.f90`` (lines 192–196, 340–360, 418–420)
- **SCF application:** ``source/self.f90`` (lines 445–449, 862–881)
- **Exchange application:** ``source/exchange.f90`` (lines 120–124, 1471–1492)
- **Transport application:** ``source/conductivity.f90`` (lines 102–106)

See Also
========

- :ref:`keywords/constraints` — Input parameters for the ``&constraints`` namelist
- :ref:`theory/scf_cycle` — SCF loop in which the constraint is applied
- :ref:`theory/spin_dynamics` — Broader context of non-collinear magnetism
