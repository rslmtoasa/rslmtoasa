.. _examples/interface_fcccu111:

=========================================================
Layered interfaces: fcc Cu(111), ``calctype = 'L'``
=========================================================

**Location:** ``example/interface/fccCu111_AA/``,
``example/interface/fccCu111_AB_misaligned/`` and
``example/interface/fccCu111_Avac/``

**System:** one fcc Cu(111) atomic layer embedded between two frozen
reference regions — two metallic ones, or metal and vacuum.

**Physics:** two-sided embedding — the layered counterpart of the one-sided
surface path. Demonstrates the region registry, deviation-variable
electrostatics, and the alignment solver that reconciles two independently
converged parameter sets.

.. warning::

   This page documents what the layered path **does** do, and is deliberately
   explicit about what it **does not**. Read :ref:`interface-limitations`
   before using ``calctype = 'L'`` for production work. One of the four scoped
   geometries is **not implemented**, one namelist knob does not mean what its
   name suggests, and the layered dipole barrier — while now physically
   meaningful — is not yet quantitatively pinned against the one-sided route.

Overview
========

RS-LMTO-ASA does everything by *embedding*: a finite real-space cluster is
built, an interaction zone is allowed to relax self-consistently, and
everything outside it is **frozen** at potential parameters imported from a
previous converged calculation. Bulk, impurity and surface modes are all
special cases of that one idea.

The layered mode adds one thing: **more than one frozen reference region may
coexist in a single cluster**. A "lead" is nothing but a set of frozen
potential parameters assigned to a subset of cluster sites. There is no new
Green-function machinery, no boundary self-energy :math:`\Sigma(E)`, and no
decimation — the recursion is entirely unaware that regions exist.

.. code-block:: text

   region A (frozen)  |  active zone  |  region B (frozen)
   ------------------- + ------------- + -------------------
   imported parameters   relaxes SCF     imported parameters
   low z                                 high z

With ``region_b_kind = 'metal'`` (the default) both boundaries are metallic
and this is *not* a surface calculation; see :ref:`interface-not-a-surface`.
With ``region_b_kind = 'vacuum'`` region B is semi-infinite vacuum and the
geometry *is* a surface — see ``fccCu111_Avac`` below.

The three examples
==================

``fccCu111_AA`` — the identity case
------------------------------------

All three regions are the same material, loaded from the same converged bulk
Cu parameter set (``vmad ~ 0``), so the active atom is physically
indistinguishable from bulk fcc Cu. Every iteration must report exactly:

.. code-block:: text

   interfacepot: Q=   0.0000E+00  P=   0.0000E+00
   interfacepot: deep-probe potentials v_lo= 0.000000 v_hi= 0.000000 step= 0.000000
   align_regions: V(B)= 0.000000 Ry
   align_regions: max alignment residual=   0.0000E+00

Verified: all of these hold to printed precision for 30 consecutive
iterations. This is the cheapest sanity check on the layered path — any
nonzero value means something is wrong.

``Q = 0`` is a **kernel precondition, not merely physics**. The existing 2D
Madelung ``AM`` matrix is built in a half-space gauge whose difference from
the symmetric kernel is a global constant *if and only if* :math:`Q = 0`; at
:math:`Q \neq 0` the gauge injects a spurious uniform field of
:math:`(2\pi/A)Q`. The code warns if this is violated.

``fccCu111_AB_misaligned`` — a negative control
------------------------------------------------

Same geometry, but the input *declares* the two regions to sit on different
absolute energy scales:

.. code-block:: fortran

   &charge
   fermi_a = -0.089874    ! the true converged Cu value
   fermi_b = -0.109874    ! a lie: 20 mRy below it
   /

while both are in fact the same parameter set. **This run is supposed to
warn**, every iteration:

.. code-block:: text

   <warning> align_regions: converged alignment shift for region B disagrees
   with the analytic contact potential E_F(anchor) - E_F(region) by 0.020000
   Ry, which exceeds the 0.0050 Ry threshold. ... treat the reported dipole
   barrier as unvalidated.

``0.020000`` Ry is exactly the injected discrepancy, recovered. The example
exists to show the check has teeth.

``fccCu111_Avac`` — a surface, with the generated vacuum lead
--------------------------------------------------------------

One active fcc Cu(111) layer with frozen Cu below and **semi-infinite vacuum
above**, selected with a single namelist entry:

.. code-block:: fortran

   &lattice
   region_b_kind = 'vacuum'   ! region B is vacuum, not a second metal
   /

   &atoms
   label(1) = 'A'      ! frozen region A reference (bulk Cu)
   label(2) = 'Vac'    ! vacuum slot — parameters GENERATED, not read
   label(3) = 'Ac1'    ! active layer
   /

The run announces the mode and dumps the generated parameter set before the
first recursion:

.. code-block:: text

   <info> buildinterface: region B is VACUUM. Its frozen parameters are
   generated per run by vacuum_lead and regenerated each iteration at the
   solved vacuum level (B7 §1.6).
   --- vacuum_lead: empty-lattice frozen parameter set ---
     ws_r  =     2.668992 bohr
     v0    =     0.000000 Ry
     band onset (C_s) =     0.346865 Ry

Both spin channels come out identical, as they must — vacuum is
non-magnetic — and the alignment diagnostics name the region ``vacuum``
rather than ``B``:

.. code-block:: bash

   grep "region B is VACUUM"     out.log   # mode confirmation
   grep "band onset"             out.log   # must sit well above E_F
   grep "align_regions: V(vacuum)" out.log

.. warning::

   This example produces a surface dipole of ``step = -0.0949`` Ry
   (``V(vacuum) = +0.0949`` Ry, about 1.29 eV), converging cleanly. The
   independent ``buildsurf`` route gives ``-0.1236`` Ry on the same system —
   same sign and magnitude, ~23 % apart. See
   :ref:`interface-electrostatics-zero` for why the two differ and what is
   still open.

Why alignment needs a solver at all
===================================

Each bulk calculation carries its own potential zero. Freezing two of them
side by side with no relative shift silently imposes

.. math::

   \text{contact potential} \equiv 0

whether or not that is true. The run converges, looks entirely plausible, and
reports a wrong dipole barrier and wrong charge transfer. This is identified
in the design as *the worst failure mode available to us*, because every
other failure mode in the feature is loud.

The physical requirement is that deep inside region :math:`r`, the site must
*be* neutral bulk :math:`r`:

.. math::

   E_F - V_r = E_F^{(r)}, \qquad
   V_B - V_A = E_F^{(A)} - E_F^{(B)}

so only differences are physical, and the shifts are fixed up to one gauge.
Region A is the gauge anchor (:math:`V_A \equiv 0`).

Crucially, :math:`V_r` is **not** computed from that closed form, because
that would require a trustworthy cross-run absolute energy scale. It is
solved as an SCF fixed point on the deep-probe residual, mixed with ``vmix``
like any other SCF quantity.

``fermi_a`` / ``fermi_b`` therefore do **not** set the alignment. They only

1. seed the initial guess, and
2. drive the consistency check.

This is why ``fccCu111_AB_misaligned`` converges to :math:`V_B = 0` despite
declaring 20 mRy: the two regions really are the same material, and the fixed
point correctly reports the physics rather than the declaration. The fixed
point is the answer; the declared Fermi levels are the check.

.. _interface-limitations:

Limitations — read this before production use
==============================================

.. _interface-electrostatics-zero:

Electrostatics: fixed, and what the numbers mean now
------------------------------------------------------

Earlier builds of the layered path reported ``Q``, ``P``, ``step`` and every
deep probe as exactly zero for *every* ``calctype = 'L'`` run, so no layered
calculation produced a physical dipole barrier. Three index-space bugs caused
it, and all three are fixed:

1. The active site's deviation charge was written to a **frozen** Madelung row
   (the active zone starts at row ``nlay_a+1``, not row 1), and the
   compensation step then subtracted the whole residual from that same row —
   exact cancellation.
2. The N(E_F) compensation weighting read the same wrong index space, returned
   zero on both sides, and fell back to a 50/50 split — putting **half the
   compensation charge into vacuum**, which does not perturb the work function
   so much as *set* it.
3. Frozen boundary rows carried no valid reference type, which is what made
   (2) silent.

The ``A | A`` identity masked all of this perfectly: its correct answer *is*
zero for all five reported quantities, so a spuriously zero result was
indistinguishable from a right one.

**After the fix**, the identity still holds exactly while the surface produces
a real barrier:

.. list-table::
   :header-rows: 1
   :widths: 30 12 20 20 18

   * - Case
     - Q
     - P
     - step
     - V(B / vacuum)
   * - ``A | A`` identity
     - 0
     - 0
     - 0
     - 0
   * - ``A | vacuum``
     - 0
     - -1.12e-2
     - -0.0949 Ry
     - +0.0949 Ry

Cross-checked against the independent one-sided ``buildsurf`` route on the same
system (fcc Cu(111), one active layer), which reports
``vmad1 = vm1 - vbulk = -0.1236`` Ry against the interface route's ``-0.0949``
Ry — same sign and magnitude, about 23 % apart. The two routes probe different
points of the same profile (``buildsurf`` uses rows 1 and ``nbas``; the
interface route uses the frozen-region extremes), so exact agreement is not
expected. Closing that gap is a B7.7 item; until then, treat the layered
barrier as physically meaningful but not yet quantitatively pinned.

.. note::

   **Leave ``&charge nlay_a`` / ``nlay_b`` unset.** They position the active
   zone within the Madelung stack, and it must be *centred*: otherwise one deep
   probe sits next to the active charge and the other ``nbas`` rows away, and
   the solver correctly reports a nonzero :math:`V_B` for a physically
   symmetric cell. Measured on ``fccCu111_AA``: ``1/1`` gives
   :math:`V_B = -0.0109` Ry, centred ``24/24`` gives exactly 0. Unset, the code
   centres it automatically and logs the choice. To thicken the frozen buffer,
   use the ``&lattice`` pair instead.

Vacuum regions: ``A | vacuum`` works, the vacuum gap does not
--------------------------------------------------------------

The feature was scoped for four geometries:

.. list-table::
   :header-rows: 1
   :widths: 32 14 54

   * - Geometry
     - Status
     - Notes
   * - ``A | A`` identity
     - **works**
     - ``fccCu111_AA``
   * - ``A | B`` metallic
     - **works**
     - Alignment machinery exercised by ``fccCu111_AB_misaligned``
   * - ``A | vacuum`` (surface)
     - **works**
     - ``fccCu111_Avac``; set ``&lattice region_b_kind = 'vacuum'``
   * - ``A | vacuum-gap | B``
     - **not implemented**
     - Needs a four-region layout; the cluster builder hardcodes three

Setting ``region_b_kind = 'vacuum'`` makes region B semi-infinite vacuum.
Its frozen parameters are **not** read from an ``&atoms`` label: they are
generated per run by ``source/vacuum_lead.f90``, which sets up an empty
sphere of the run's own Wigner–Seitz radius as an ordinary :math:`Z = 0`
atom and pushes it through the code's own radial solver. The representation,
the :math:`\gamma/o` convention, the :math:`E_\nu` convention and the
normalization therefore come out automatically consistent with every other
parameter set in the run. The analytic spherical-Bessel result is the test
*oracle*, never the implementation.

The vacuum level :math:`V_0` is **self-consistent, not an input**. It is the
alignment shift the solver already maintains for every non-anchor region,
and the vacuum parameters are regenerated at the updated level after every
alignment step. Iteration 0 runs at :math:`V_0 = 0`, so the older
"generate once from a given vacuum level" behaviour arrives as the first
step of the self-consistent scheme rather than as a separate code path.

.. note::

   The ``Vac.nml`` file in the example is a *slot*, not a data source. Its
   potential parameters are overwritten by the generator before the first
   recursion and after every alignment update; only its element block
   (:math:`Z = 0`, no valence) is meaningful. It exists because the
   ``&atoms`` reader requires one file per declared label.

What the vacuum lead will **not** do is change converged good-metal surface
numbers much. For :math:`\Phi \approx 4` eV the evanescent decay constant is
:math:`\kappa = \sqrt{\Phi} \approx 0.54\ a_0^{-1}`, so vacuum tails die
within about 2 Å — which is exactly why one or two hand-set empty-sphere
layers already worked. The acceptance criterion is *reduced sensitivity to
the vacuum layer count*, not a changed answer. The regimes where it genuinely
earns its keep are vacuum-gap interlayer exchange, tunnelling barriers and
low-DOS spacers, where the tail is the physics rather than a boundary detail.

One honest limitation: LMTO-ASA reproduces empty-lattice free-electron
dispersion only over a limited window above the band bottom, so the frozen
vacuum parameters are meaningful only while :math:`E_F` sits well below the
vacuum band onset. That holds comfortably for ordinary metals; it is
uncomfortable for very low work functions. The run checks this margin every
time it regenerates and warns when it is thin.

Two different ``nlay_a`` / ``nlay_b``: widen the buffer with ``&lattice``
--------------------------------------------------------------------------

The same two names appear in two namelists and mean different things. This is
deliberate, and it is the sharpest trap on this path.

- **``&lattice nlay_a`` / ``nlay_b`` count atomic layers.** This is the
  frozen-buffer width, and the knob to reach for. ``build_interface_full``
  bins cluster sites into layers on exactly these counts.
- **``&charge nlay_a`` / ``nlay_b`` count rows of the synthetic 2D Madelung
  stack**, whose size is a fixed constant (``nbas = 49`` for this geometry)
  deliberately decoupled from the physical layer count so the ``set2d``
  NLAMA/NLAMB split stays balanced.

For a one-atom-per-layer cluster the two happen to coincide at 1; in general
they do not. Raising the ``&charge`` values relabels *interior* Madelung rows
— rows that carry real deviation charge — as frozen boundary, which moves the
deep probe onto a charged row and breaks the fixed point.

Measured on the ``fccCu111_AA`` geometry, varying one at a time:

.. list-table::
   :header-rows: 1
   :widths: 26 24 26 24

   * - ``&lattice`` a/b
     - ``&charge`` a/b
     - converged ``V(B)``
     - identity
   * - ``1 / 1``
     - ``1 / 1``
     - ``0.000000`` Ry
     - holds ✔
   * - ``4 / 4``
     - ``1 / 1``
     - ``0.000000`` Ry
     - holds ✔
   * - ``1 / 1``
     - ``4 / 4``
     - ``-0.4498`` Ry
     - broken ✘

So a thicker frozen buffer **is** available: raise ``&lattice nlay_a`` /
``nlay_b`` and leave the ``&charge`` pair at 1.

The identity is exact only for a single active layer, though. Widening the
*active* zone to ``nlay = 3`` (with ``&lattice nlay_a/nlay_b = 4``) stays
stable — ``Q = 0`` exactly, alignment residual ~1.8e-5 — but ``V(B)`` settles
near **7.6 mRy** rather than 0, and ``P`` at ~-1.1e-3. That is the imperfect
screening of a three-layer active zone between frozen boundaries, not a
divergence; it is the size of the effect §1.5 predicts when
:math:`\sum_{\text{active}} \delta q` fails to vanish. Quantifying how it
shrinks with active-zone thickness is exactly the compensation-insensitivity
study deferred to B7.7 / gate G-B7-3, and it has not been done.

The drift warning is noisier than its wording suggests
--------------------------------------------------------

With a free Fermi level, "deep-A is neutral bulk" is not imposed — it is
inherited from a sufficiently thick frozen buffer. The code reports when that
inheritance looks doubtful:

.. code-block:: text

   <warning> align_regions: deviation charge ... at frozen anchor-region site 1
   exceeds the 1.00E-03 threshold. ... this means the active zone is too thin
   AND that the alignment gauge is inconsistent. Widen the active zone.

Follow it via ``&lattice``, per the table above. Note however that the warning
still fires in runs where the buffer has been widened and the identity holds
*exactly* (``V(B) = 0``, residual ``0``), so on this path it is better read as
"check your buffer" than as proof that anything is wrong. Recalibration of the
``deep_drift_tol`` threshold is deferred to B7.7 / gate G-B7-3, where the
tolerances accepted "for now" in G-B7-2 are revisited.

Orientation-dependent residual on (111)
-----------------------------------------

Cu treated as an interface between Cu regions must reproduce bulk fcc Cu.
``surftype = '0 0 1'`` does so **exactly**, at every sampled DOS row.
``surftype = '1 1 1'`` deviates by **2.05e-3** near the d-band peak — roughly
200× the deviation of the impurity route, against a peak DOS of 48.3.

The Hamiltonian is *not* the cause: the on-site block and all 19 fcc
neighbour hoppings are bit-identical across the bulk, impurity and layered
routes for both orientations, and ``etot`` agrees to ~1e-8 Ry. The residual
is downstream of the Hamiltonian and specific to the (111) normal. It is
unexplained, deliberately pinned by the committed test reference rather than
tolerated away, and tracked in ``tests/KNOWN_ISSUES.md``.

The examples on this page use (111) — the orientation with the known
residual — because that is what the surface example they mirror uses. For the
cleanest identity check, ``surftype = '0 0 1'``.

What the vacuum lead is, and is not, for
=========================================

Even once wired, the vacuum lead is **not expected to change converged
good-metal surface numbers**, and it would be a mistake to adopt it hoping
for that. For a typical metal the evanescent decay constant is
:math:`\kappa = \sqrt{\Phi} \approx 0.54\ a_0^{-1}` at :math:`\Phi \approx
4` eV, so vacuum tails die within about 2 Å. That is precisely why one or two
hand-set empty-sphere layers already works.

The acceptance criterion for the vacuum lead is therefore **reduced
sensitivity to the empty-sphere layer count, not a changed answer**. Its real
payoff is a well-defined vacuum level, removal of an ad hoc convergence knob,
and geometries where the tail genuinely matters: vacuum-gap interlayer
exchange, tunneling barriers, and low-DOS spacers.

One further honest caveat: LMTO-ASA reproduces empty-lattice free-electron
dispersion only over a limited window above the band bottom. Frozen vacuum
parameters are meaningful only when :math:`E_F` sits well below the vacuum
band onset — true for ordinary metals, uncomfortable for very low work
functions.

Expected physics limitations
=============================

Long-period interlayer-exchange oscillations arise from fine Fermi-surface
features, and real-space recursion resolves those only with long recursion
chains. Fe/Cr and Co/Cu spacer-thickness benchmarks are the right validation
targets, but they require a convergence study on the *oscillation period* in
the recursion-level count ``lld``, not only on the layer moments. A period
discrepancy should not be treated as a bug before that study exists. Neither
the benchmarks nor the study have been run.

The layered path also delivers interface **electronics and magnetics** only —
charge transfer, dipole barrier, layer-resolved moments and DOS, and
interlayer exchange coupling. It does **not** deliver Landauer conductance;
the word "lead" here carries no transport machinery.

.. _interface-not-a-surface:

This is not a surface calculation
==================================

For the one-sided ``vacuum | active | bulk`` geometry, use
``pre_processing = 'buildsurf'`` (:doc:`surface_fcccu001`). ``buildsurf`` is
untouched by the interface work and remains the permanent regression oracle
for the surface case. ``example/surface/fccCu111/Single`` is the direct
one-sided counterpart of the examples on this page.

Running the examples
====================

.. code-block:: bash

   cd example/interface/fccCu111_AA
   ../../../build/bin/rslmto.x

Each SCF iteration takes ~10 s on one core; the 30 iterations configured take
roughly five minutes.

**What to check:**

.. code-block:: bash

   grep "interfacepot: Q"       out.log   # Q and P must be 0
   grep "align_regions: V(B)"   out.log   # must be 0 for A|A
   grep "max alignment residual" out.log  # must be 0 for A|A

See Also
========

- :doc:`surface_fcccu001` — the one-sided surface path
- :doc:`impurity_b2feco` — embedded-cluster calculation
- ``docs/dev/plans/B7_interfaces_and_vacuum_leads.md`` — design and rationale
- ``docs/dev/CONTRACT_FROZEN_REGION.md`` — what a parameter set must persist
  to be usable as a frozen region
- ``tests/KNOWN_ISSUES.md`` — the (111) residual
