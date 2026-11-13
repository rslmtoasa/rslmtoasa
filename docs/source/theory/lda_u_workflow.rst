.. _theory_lda_u_workflow:

=======================================
LDA+U and LDA+U+V Workflow
=======================================

Purpose
=======

This page documents the current LDA+U and LDA+U+V implementation in RS-LMTO-ASA:

1. Reference equations and conventions used in development.
2. Developer-oriented mapping from equations to Fortran routines and variables.
3. Practical ``input.nml`` examples for the supported modes.

Primary authorship credit
=========================

The LDA+U/LDA+U+V development documented here is based on work carried out by
**Viktor Frilen** and **Emil Beiersdorf** (project documentation, June 2024).
Please keep this attribution when reusing or extending this workflow description.

Credits and references
======================

This implementation and documentation build on the internal project document by
Viktor Frilen and Emil Beiersdorf (June 2024), and on the standard literature:

- Liechtenstein, Anisimov, Zaanen (rotationally invariant DFT+U).
- Dudarev et al. (``U_eff = U - J`` form).
- Campo and Cococcioni (extended DFT+U+V form).
- Agapito, Curtarolo, Nardelli (ACBN0 framework).

Reference equations
===================

Rotationally invariant on-site interaction (Liechtenstein form):

.. math::

   E_U = \frac{1}{2}\sum_{m m' m'' m''' \sigma}
   \left[
   \langle m m''|V_{ee}|m' m'''\rangle n^{\sigma}_{m m'} n^{-\sigma}_{m'' m'''}
   +
   \left(
   \langle m m''|V_{ee}|m' m'''\rangle
   -
   \langle m m''|V_{ee}|m''' m'\rangle
   \right)
   n^{\sigma}_{m m'} n^{\sigma}_{m'' m'''}
   \right].

Double counting (FLL):

.. math::

   E_{dc} =
   \frac{U}{2}n(n-1) -
   \frac{J}{2}\sum_{\sigma} n^{\sigma}(n^{\sigma}-1).

Resulting potential matrix (on-site channel):

.. math::

   V^{\sigma}_{m m'} =
   \frac{\partial (E_U - E_{dc})}{\partial n^{\sigma}_{m m'}}.

Simplified Dudarev form:

.. math::

   E_U - E_{dc} = \frac{U}{2}\sum_{\sigma}
   \mathrm{Tr}\left[n^{\sigma}(1-n^{\sigma})\right].

The validated implementation limit is ``J=0`` with nonzero ``U``; the full
Liechtenstein U/J implementation is described above.

Extended +V correction (nearest-neighbor inter-site contribution):

.. math::

   E_{U+V} - E_{dc} =
   \sum_{I,\sigma}\frac{U^I}{2}\mathrm{Tr}[n^{I\sigma}(1-n^{I\sigma})]
   -
   \sum_{I,J,\sigma}^{*}\frac{V^{IJ}}{2}\mathrm{Tr}[n^{IJ\sigma}n^{JI\sigma}].

Implementation status
=====================

- ``liechtenstein`` potential path:
  the ordinary on-site U/J potential is wired and is the validated scope of
  [VAL-04]. Validated scope: KS and RS bulk energies for fixed d-shell U/J;
  RS impurity density only. Surfaces are experimental. The Dudarev validation
  limit is ``J=0`` with ``U``; this does not validate arbitrary J.
  Chebyshev recursion with fixed U is validated for collinear bulk and impurity runs: the
  on-site Green function and occupation matrices are deterministic; a bcc Fe block-versus-Chebyshev
  comparison at U=2 eV gave an E_F difference of 9.75e-4 Ry (limit 1.84e-3) and an E_Hub relative
  difference of 4.1e-3 (limit 1e-2); and in local mode Tr(M1)=q0*(c'+t/q0) holds to 1e-15 in a
  Chebyshev bulk and a magnetic Ni impurity run (see ``tests/KNOWN_ISSUES.md``).
- ``acbn0`` potential path:
  optional path is available, but should currently be treated as an
  implementation-in-progress path for validation studies.
- ``+V``:
  **Development**.  The current nearest-neighbor orbital-channel correction
  uses local occupations and a practical diagonal-in-m approximation.  This
  is a **proxy** implementation: the full inter-site matrix form
  :math:`-V^{IJ} n^{JI,\sigma}_{m' m}` is not yet active and must not be
  described as a validated inter-site DFT+U+V implementation.

C5 impurity radial-moment mode
-----------------------------

.. note::

   ``hubbard_radial_moments='auto'`` now always resolves to ``spectral``, the standard LDA+U construction
   (moments of the Green function of the full H = H_LDA + V_U), which also works with ``hoh``, ``ccor_2c``,
   impurities and surfaces. ``local`` (moments of the LDA part only, V_U's own contribution removed; not the
   standard construction) is an explicit opt-in: set ``hubbard_radial_moments = 'local'`` in ``&hamiltonian``.
   It is the only route with the routed total energy, at first order (no ``hoh``) for the bulk and for impurity and surface clusters (cluster sites only, host types excluded).
   The sections below that say ``auto`` selects ``local`` describe the rule before 2026-11-07; the same
   scope and host-provenance conditions apply when ``local`` is requested explicitly.

For fixed Liechtenstein U/J in the supported scalar-collinear d-shell scope,
``hubbard_radial_moments='auto'`` selects ``local`` for an RS impurity only
when every U-carrying frozen host potential records ``local`` mode. This
applies to both an RS-local host feeding an RS impurity and a KS-local bulk
host feeding an RS impurity. Impurity U/J may be supplied through the general
``hubbard_u_general``/``hubbard_j_general`` fields or through the
``hubbard_u_impurity``/``hubbard_j_impurity`` override fields. The selected
radial mode is independent of the host's electronic-structure backend.

If any U-carrying frozen host records ``spectral`` mode, or lacks mode
provenance, ``auto`` selects ``spectral`` and logs the reason. A legacy host
without ``hubbard_radial_moments`` has unknown mode; it is treated as
``spectral`` and emits a warning asking the user to restart the host once with
the current version so its mode is recorded. Updated impurity potentials are
written with the run's selected mode and backend. A forced explicit
``spectral`` choice remains spectral.

For impurity decks, copy a working deck, then set ``fermi`` and the energy
window to the host endpoint values before launch. The C5.4 preflight compares
the deck Fermi value with the host run output's endpoint ``E_F`` and requires
a difference no greater than ``1e-6 Ry``. It records the host output path and
value, and rejects a mismatch before the solver starts.

With ``hubbard_radial_moments='local'`` at first order, ``E_Hub`` of the cluster sites' P_out
(host types excluded) is added to ``physical_total_energy`` as in the bulk
(``tests/ldau/test_cluster_energy_u.py``: ``nsp=3`` against ``nsp=1``); with ``spectral`` or ``hoh`` impurity energies
exclude the Hubbard term. The energy formula is unchanged. C5's numerical
KS-host-to-RS-impurity equivalence for Ni remains unvalidated, and RD4 has no
energy measurement: its optimized Release diagnostic has only been checked
through one complete iteration. See the C5 records in
``docs/validation/fixed_ldau/VALIDATION_SUMMARY.md`` and
``C5_10_RESULTS.md``.

C6 surface density-only mode
-----------------------------

For ``calctype='S'`` in the same supported fixed Liechtenstein, scalar-collinear
d-shell scope, ``hubbard_radial_moments='auto'`` selects ``local`` when every
U-carrying frozen host potential records ``local`` mode. This applies to an
RS-local or KS-local bulk host. The surface recursion layers remain the active
surface density sites; the host mode is resolved from the frozen-host ownership
map. A legacy U-carrying host with missing mode provenance is treated as
``spectral`` and does not qualify for automatic local mode. Explicit
``spectral`` remains available.

To prepare a surface deck, copy an existing surface fixture. Replace only the
frozen host potential with the accepted local host potential, set U/J, and set
``fermi`` and the energy window to the host endpoint values. Run the host-E_F
preflight before launch. Use the same deck for RS and KS surface calculations;
the local density mode is independent of the host's backend.

Surface Hubbard support is density-only. With ``hubbard_radial_moments='local'`` at first order the surface energy includes ``E_Hub`` of the
surface-layer P_out (host types excluded; ``tests/ldau/test_cluster_energy_u.py``); with ``spectral`` or ``hoh`` it excludes the Hubbard term. Surface density support is experimental. A Ni(001) probe showed a narrow-window
non-finite/NaN projected-DOS observation; widening the energy window did not
remove it. See the C6 records in
``docs/validation/fixed_ldau/VALIDATION_SUMMARY.md`` and
``C6_AMENDMENT.md``.

SOC+U for bulk (nsp=2)
======================

Fixed Liechtenstein U/J with spin-orbit coupling (``nsp=2``) runs by default for ``calctype='B'`` clean bulk,
``lmax=2``, in both the RS and the k-space (KS) solvers. The d shell is treated as a 10x10 spinor block: the
density matrix P and the first energy moment M1 keep the spin-off-diagonal blocks, the potential V and the
energies come from the spinor functional, and the radial moments use the spinor local correction (``auto`` resolves to
``local``). The Hubbard energy of the spinor P_out is added to ``physical_total_energy`` (candidate b). The former
switch ``hubbard_soc_experimental`` is accepted and ignored (info line "SOC+U bulk is now default").

Restrictions (each stops with a fatal): ``calctype`` ``S`` with U, self-consistent U, intersite V, HOH with
k-space ``kspace_ham_order='first'``, onsite U with ``kspace_ham_order='orthogonal'``, with HOH for ``calctype`` other than ``B`` outside ``nsp=1`` without SOC and the ``nsp=2``, ``3`` and ``4`` spinor routes,
at ``nsp=2`` correlated sites whose spin axes are not parallel or antiparallel (message "noncollinear axes: use nsp=3 or 4"),
and RS block recursion with a spin axis not parallel to global z (see "Spin axis" below).

``nsp=3`` and ``4`` with U are accepted (spinor U as at ``nsp=2``; with HOH the V is added to ``enim``). They are tested
for collinear setups only: bcc Fe ferromagnet, a two-atom Fe antiferromagnet and NiO AFM-II, where results agree with the
``nsp=1`` (``nsp=3``) and ``nsp=2`` (``nsp=4``) routes. Tilted and noncollinear states are checked by symmetry only
(global spin rotation, exchange of the two sublattices, cubic equivalence of the spin axis with SOC, restart round trip; no external
reference), with Chebyshev recursion and k-space ``second``. RS block recursion at ``nsp>=3`` with noncollinear axes is not
rotation invariant with the default terminator and nearly so (1e-4 muB) with ``sym_term = .true.`` (``tests/KNOWN_ISSUES.md``). ``local_axis``, ``gbt_single_q``, constraints
and runs with U at ``nsp>=3`` on surfaces other than the Ni(001) test are untested; ``local`` at ``nsp=3`` on impurity clusters, with and without hoh (spinor correction with V by type as at ``nsp=2``; with hoh V_U is in ``enim``), is checked against ``nsp=1`` by ``tests/ldau/test_impurity_nsp3_u.py`` (B2FeCo, nsp=3 k-space host); the nsp=3 surface route is checked on the Ni(001) surface below. ``hubbard_radial_moments='local'`` at ``nsp=3`` (bulk;
real space and k-space) and ``nsp=4`` (bulk real space and k-space, impurity clusters) takes the spinor correction along each site axis, with any site axes (``tests/ldau/test_local_nsp4_u.py``: ``nsp=4`` with SOC held at zero against ``nsp=3``, a ferromagnet with SOC against ``nsp=2``, and z against x); ``tests/ldau/test_local_nsp4_workflows_u.py`` checks ``nsp=4`` in k-space (antiferromagnet with SOC off against ``nsp=3``, with and without hoh; ferromagnet with SOC against ``nsp=2``) and on the B2FeCo impurity (SOC off against ``nsp=3``, with and without hoh; surfaces: ``tests/ldau/test_surface_nsp3_u.py`` checks ``nsp=3`` against ``nsp=1`` on the Ni(001) surface, with and without hoh, with ``density_policy='constrained_spiral'`` (the default policy differs once the axis of a small-moment site flips, ``tests/KNOWN_ISSUES.md``); ``tests/ldau/test_surface_nsp4_u.py`` compares ``nsp=4`` with SOC held at zero against ``nsp=3`` on the same surface, with and without hoh (largest difference 0 over 4 iterations and 112 quantities; the SOC matrix is zero and the correction code is shared for ``nsp>=3``, so this does not test the correction itself; ``nsp=4`` surfaces with SOC on are not tested)); ``tests/ldau/test_local_noncollinear_u.py`` checks that rotating or exchanging a noncollinear pair leaves E_F and the site moments unchanged. The two-atom Fe antiferromagnet reproduces the ``nsp=1`` ``local``
route (``tests/ldau/test_nsp34_u.py``, ``af_fe_local``, real space; ``tests/ldau/test_kspace_local_nsp3.py``, k-space, antiferromagnet and ferromagnet). The Hubbard energy of the spinor P_out is added to ``physical_total_energy`` for bulk ``local`` at ``nsp=3`` and ``4`` as at ``nsp=2`` (``tests/ldau/test_local_energy_u.py``, real space: ``nsp=3`` antiferromagnet against ``nsp=1``, ``nsp=4`` with SOC held at zero against ``nsp=3``, a ferromagnet with SOC ``nsp=4`` against ``nsp=2``; k-space energies are not tested); impurity and surface energies are routed at first order with ``local`` (``tests/ldau/test_cluster_energy_u.py``: ``nsp=3`` against ``nsp=1``; ``hoh`` is not routed; ``nsp=4`` clusters take the same route and are not tested). ``local`` with HOH removes V_U from the moments of the orthonormal-basis Green
function in the same way as at first order (H = H0 + V_U with V_U additive in ``enim``); the nsp=3 against nsp=1 check
``af_fe_local_hoh`` tests the spin frame, not the correction. On impurity (``calctype I``) and surface (``S``) clusters with HOH (``nsp=1``, and the ``nsp=2`` spinor route) V_U is captured at the ``enim`` insertion, as at the ``ee`` or ``hall`` insertion without HOH; ``tests/ldau/test_impurity_hoh_u.py`` compares the ``nsp=2`` run with SOC held at zero against ``nsp=1`` for the B2FeCo impurity (the surface route with HOH is compared the same way on the Ni(001) surface by ``tests/ldau/test_surface_hoh_u.py``, ``nsp=2`` with SOC held at zero against ``nsp=1``; ``tests/ldau/test_impurity_equals_bulk.py`` checks that an Fe impurity in an Fe host with the bulk Fermi energy reproduces the bulk moment, see ``tests/KNOWN_ISSUES.md``). To set up a collinear antiferromagnet use the noncollinear
convention: the same potential file (majority along the site axis) for both sublattices with ``mom = 0, 0, 1`` and
``0, 0, -1``. Potential files in the ``nsp=1`` convention (the second sublattice with up and down swapped, ``mom = +z``) run at
``nsp>=3`` with the default ``density_policy='relaxed_reference'`` silently turn ferromagnet-like after the first step
(measured, also at U=0); with ``density_policy='constrained_spiral'`` they reproduce ``nsp=1``.
Explicit ``hubbard_radial_moments='spectral'`` is accepted and its energy is not routed (as at nsp=1).

Spin axis (nsp=2)
-----------------

The magnetisation of an ``nsp=2`` run is collinear along one axis n, which need not be z. P, M1 and V stay in the
global spin frame; the local-moment correction projects them with the spin projectors Pi_pm(n) = (1 +- n.sigma)/2 (x 1 on
the orbital index). For n = z the kernel is bitwise equal to the earlier z-only loop; for n = -z it fixes a channel swap
the earlier loop made (it swapped the up and down channels).

The axis of a site is, in order: its ``constraint_target`` if non-zero, else the potential-file ``mom``, else z. With
``density_policy='constrained_spiral'`` the axis is that reference; with ``'relaxed_reference'`` it is the current moment divided by its length. Both policies, and a ``constraint_target``, are accepted. All correlated sites must share one axis (n and
-n count as one), checked once where V is built. The axis of each site is logged as a ``HUBBARD_AXIS`` info line the
first time and whenever it moves by more than 1e-8.
The potential file of a U-carrying host records the axis as ``hubbard_axis`` (written only when non-zero); a run along a
different axis stops with a fatal, and a legacy file without the field is read as z with an info line. At
``hubbard_log_level='detailed'`` the ``HUBBARD_SPINOR_DM`` line adds ``axis``, ``n_up``, ``n_dn`` (the d-shell trace of
P along n, up and down) and ``m_L``, the d-shell orbital moment along n. The line is printed with U, or at U=0 with
``hubbard_spinor_dm_diagnostic``.

Measured rotation invariance, bcc Fe, U=2/J=0, KS, ``soc_scale=0`` with ``fix_soc``, axes z, x, [111], -z, full 12x12x12
mesh: the F_b difference against z is 6.2e-11 (x), 5.6e-10 ([111]) and 8.3e-10 Ry (-z); the E_Hub difference is 1.2e-11 or smaller; the d
occupations along n agree to the printed 8 digits. With SOC on the x and z runs agree in F_b to 2.2e-10 Ry (the cubic
axes are equivalent, so this is the noise floor of the comparison) and the d occupations along n to the printed 8 digits.

Constraint targets set the axis. The constraint controller itself did not converge in the one KS smoke deck tried
(``tests/KNOWN_ISSUES.md``); use a constraint target for the axis only until that is triaged.

RS recursion and tilted axes. RS ``recur='chebyshev'`` accepts any shared axis. RS block recursion is not invariant under
a rotation of the spin axis, even at U=0 and SOC to zero (iteration-1 dE_F 3e-4 Ry, dmoment 1.8e-4 muB, x against z;
with ``sym_term`` the converged differences are still 2.4e-5 Ry and 6.1e-5 muB); it is therefore not supported with U for
a tilted axis (fatal), and RS block recursion along z is unchanged. See ``tests/KNOWN_ISSUES.md``.

RS Chebyshev is not exactly spin-rotation invariant: d occupations and moments differ by up to ~5e-6 between equivalent
axes; total energies by ~1e-7 Ry (bulk) and ~2e-6 Ry (impurity). Anisotropy energies below this floor are not
resolvable in RS; use KS (invariant to ~1e-10) for strict anisotropy comparisons. Block recursion with a tilted axis is not
supported with U (see KNOWN_ISSUES).

Measured floor, z against x at ``soc_scale=0`` with ``fix_soc``, RS Chebyshev, converged to ``conv_thr`` 1e-8 (``etot`` is
the ``etot`` of the output potential file, the only total energy printed at U=0; candidate b is printed only with U):

===================  ====  ===========  ===========  ===========  ============
system               U     d_up         d_dn         moment       etot (Ry)
===================  ====  ===========  ===========  ===========  ============
Fe bulk              0     3.9e-8       3.3e-7       3.1e-7       3.3e-7
Fe bulk              2     1.34e-6      1.08e-6      3.3e-6       2.2e-7
Ni host              0     1.30e-6      1.22e-6      2.3e-6       4.2e-7
Ni host              3     1.74e-6      2.55e-6      4.7e-6       5.6e-8
Ni impurity          0     3.0e-8       5.4e-7       1.1e-6       2.0e-6
Ni impurity          3     1.48e-6      4.6e-7       1.2e-6       1.6e-6
===================  ====  ===========  ===========  ===========  ============

Candidate b, z against x: Fe U=2 1.2e-7 Ry (3.7e-8 after a restart at beta 0.005 and conv_thr 1e-10, where etot gives
7.1e-8), Ni host U=3 2.2e-8 Ry, Ni impurity U=3 1.7e-6 Ry. The densities are the d shell along n from the output ``ql``;
the moment is the sum of ``ql`` up minus down. E_F agrees to the printed 6 digits in every pair.

Magnetic anisotropy. KS [111] against [001], U=2/J=0, SOC on, 12x12x12 mesh, with no further k-point refinement:
F_b([111]) - F_b([001]) = -4.8e-7 Ry, against 2.2e-10 Ry for x against z. The d occupations along n differ by 5.5e-5 (up)
and 8.4e-5 (down), the moment by 1.05e-4 muB. This is a record on one mesh, not a converged anisotropy energy: anisotropy
energies of this size need dense k-meshes, and the force theorem (anisotropy from a non-self-consistent rotation) is not
tested.

Validation numbers (bcc Fe, U=2 eV; details in ``docs/validation/fixed_ldau/SPINOR_DM_CONTRACT.md``). With
``soc_scale=0`` and ``fix_soc``, nsp=2 reproduces nsp=1: KS candidate b differs by 5.8e-10 Ry and E_Hub by 2.5e-13 Ry.
Hellmann-Feynman checks at nsp=2 with SOC on, dF_b/dU and dF_b/dJ against dE_Hub/dX at fixed spinor P: the difference
is 1.4e-4 (KS, U), 3.9e-4 (KS, J) and 2.4e-4 (RS, U) Ry/eV against allowances of 2.7e-3, 3.1e-3 and 4.8e-3 Ry/eV; the
spectral candidate misses by 53, 42 (KS) and 30 (RS) allowances. The RS J derivative was not tested. The base ASA
functional's W1 derivative residual (-0.14% KS, -0.16% RS) bounds the claim, as at nsp=1.

SOC+U for impurities and surfaces (nsp=2, calctype I and S)
----------------------------------------------------------

Impurity (``calctype='I'``) and surface (``'S'``) runs with U accept ``nsp=2``. The density is spinor as in the bulk;
with ``local`` at first order the impurity and surface energies include ``E_Hub`` of the cluster sites' P_out as at nsp=1; otherwise they exclude the Hubbard term.

Host requirements. Every U-carrying frozen host must be local mode, written at ``nsp=2``, and carry the spinor
density matrix (``spinor_d = 1`` in its ``_dm`` file). The potential file records ``hubbard_nsp``; a host written at another
nsp stops with "Frozen U-carrying host ... was written at nsp=X, run nsp=Y. Regenerate the host at nsp=Y." (legacy files
without the field: nsp is inferred from ``spinor_d`` with an info line). Both host workflows work: an RS bulk run
(``calctype='B'``) or a KS run (``use_kspace``) at ``nsp=2`` with U, local mode. ``fermi`` and the energy window must equal the host
endpoint, as for collinear impurities. ``auto`` resolves to ``local`` for these impurities.

Spin axis: as for the bulk, with the host axis matching the impurity axis; RS impurities along a tilted axis need
``recur='chebyshev'``.

Mixing: ``ldm_beta`` mixes the spinor P that builds V (full 10x10, same parameter and history as the collinear ldm).
Before this change it had no effect at nsp>=2. Ni impurities at U=3/J=0.7 oscillated without it.

Measured (RS Chebyshev Ni, U=3/J=0.7, SOC on; host and impurity runs read from the potential-file ``ql``; records in
``tests/KNOWN_ISSUES.md``): impurity minus bulk, d up 6.9e-4, d down -1.9e-3, spin moment 2.1e-3 muB, orbital moment
-1.8e-4 muB, against 1.3e-3, -2.8e-3, 3.7e-3 and -1.1e-4 at U=0. With ``soc_scale=0`` and ``fix_soc`` the nsp=2 impurity
agrees with the nsp=1 impurity to 9e-9 (d up), 1.9e-8 (d down) and 2.6e-7 muB (moment).

Surfaces (RS, ``recur='chebyshev'`` or ``'block'`` along z; tilted axes need Chebyshev) follow the impurity rules: every
U-carrying frozen host written at ``nsp=2`` in local mode with the spinor ``_dm`` file and a matching axis, ``fermi`` equal to the host E_F,
``ldm_beta=0.3``, and the energy window checked against the U-shifted d band. The empty sphere carries no U. ``auto`` resolves to ``local``
under the impurity host rule; ``nsp=3/4`` with U was not run for surfaces, and KS surfaces are not supported. Seed a surface from the converged
``nsp=1`` surface (zero spin coherences) or, for SOC on, expect the spin-orbit transient: a host converged without SOC moves when SOC is
switched on (``ni_host_soc`` in ``decks` (archive tag `archive/fable_v4b_u-2026-11-05`)`).

Surface energies exclude the Hubbard term, so U enters them only through the density.

Measured (Ni(001), U=3/J=0.7, short runs; per-step tables in the ``ni001_soc_u3``, ``ni001_soc_u0`` and ``ni_host_soc`` READMEs): with
``soc_scale=0`` and ``fix_soc`` the nsp=2 surface reproduces the nsp=1 surface at iteration 1 (``ql`` 6.9e-18, Diff 3.1737e-6 in both;
a host converged independently at nsp=1 differs by 2.5e-10 in ``ql``). A surface along x against the exactly rotated z state differs at
iteration 1 by 7.0e-9 in ``ql`` and 1.5e-8 Ry in ``etot``; no converged rotation floor or surface anisotropy was computed. With SOC on, 5
steps from the nsp=1 seeds: U=0, orbital moments flat from step 2 (Ni-S 0.0647, Ni-S-1 0.0402 muB); U=3, the Ni-S orbital moment rises
0.066 -> 0.109 muB with decreasing increments and Diff stays at 1.2e-4 to 1.9e-4, as the (unconverged) bulk host does. Converged
surface-versus-bulk values are not available.
A 20-step continuation at U=3 (SOC on) from the step-5 state, with the converged bulk host (59 iterations, orbital moment 0.0751 muB),
shows Diff rising from 1.2e-4 to 2.1e-4, a jump in the layer occupation at step 17, and a surface orbital moment still rising
(Ni-S 0.118 -> 0.175 muB). The same happens with J=0 on the layers, so it is not attributed to J; it is recorded as an undiagnosed
surface SCF instability (``tests/KNOWN_ISSUES.md``). Frozen hosts must be converged; the testing policy is
``docs/validation/fixed_ldau/TESTING_POLICY.md``.

Validated ordinary on-site scope
================================

Base ASA stationarity measurement
---------------------------------

The fixed-shift U=0 control has a base-energy derivative residual of −0.14% in
KS and −0.16% to −0.17% in RS. This is a property of the ASA energy assembly
and bounds the accuracy claimed for the matching fixed-U derivative tests.
It does not imply that explicit spectral mode adopts the local Hubbard-energy
correction; that mode remains spectral/non-adopted. See
``docs/validation/fixed_ldau/ENERGY_ACCOUNTING.md`` for scope and receipts.

Fixed-U energy-accounting status (A4.4)
---------------------------------------

For native reciprocal KS with fixed Liechtenstein U/J and scalar-collinear
``nsp=1`` on clean bulk structures (``nbulk=0`` and ``nrec=ntype``), the
current numerical M2.2/M2.3 result supports candidate
``b = F_local + E_Hub(P_OUT)`` at the stated accuracy: “valid to the ASA
base-functional accuracy measured in A5.1.” Here ``F_local`` is the local-mode
base free energy, including the KS entropy term, and ``E_Hub`` is evaluated by
the production Hubbard functional on the actual output density matrix
``P_OUT``. The U=2,J=0 and J=0.5/0.9 eV derivatives lie within the
pre-registered charge-only ASA allowances; the spectral candidate ``c`` lies
outside them. Candidate ``a`` is reported without an a-versus-b claim.

The tested derivative separates into the fixed-density Hubbard response and
the change of the captured density:

.. math::

   \frac{d b}{d x} =
   \left.\frac{\partial E_{\rm Hub}(P;x)}{\partial x}\right|_{P=P_{\rm OUT}}
   + \mathrm{Tr}\!\left[V_{\rm Hub}(P_{\rm OUT};x)
     \frac{dP_{\rm OUT}}{dx}\right]
   + \frac{dF_{\rm local}}{dx}, \qquad x\in\{U,J\}.

The first term is measured independently by applying the production
``evaluate_hubbard_energy`` function to the direct captured ``P_OUT`` at the
center parameter and nearby parameter values. The remaining response is
checked against the charge-only ASA allowance derived from finite-W1 spin
response and measured charge changes. Agreement is numerical and
tolerance-limited for these KS states; it does not assert that the last two
terms vanish or that the relation holds for arbitrary radial states.

This evidence does not establish an exact general ASA energy identity. The
selected formula is an empirically tested local base-functional correction
within the measured numerical accuracy, and does not prove a universal
``SUMEV`` double-counting coefficient. The A4.4 source routing, restart
continuity, and matched serial/MPI checks pass for native reciprocal KS. The
reviewed B3 scope removes the backend selector from the shared local-Hubbard
eligibility predicate, so qualifying clean-bulk RS and KS runs can use the
same local correction and include ``E_Hub(P_OUT)`` in the physical total
energy. The retained eligibility conditions are ``calctype='B'``, clean bulk
(``nbulk=0`` and ``nrec=ntype``), fixed Liechtenstein U/J, scalar-collinear
``nsp=1``, ``lmax=2``, and disabled constraints. Accepted-endpoint RS route reproduction and unchanged-mode/no-U controls are
recorded in ``docs/validation/fixed_ldau/VALIDATION_SUMMARY.md``; this is
route evidence, not a new ASA derivative-gate result. Explicit spectral selection remains
spectral/non-adopted. Impurity/surface hosts, noncollinear states, nsp=2 outside the SOC+U bulk scope (next section), HOH,
self-consistent U, intersite V, impurity U, and nonzero test fields remain
excluded. Constraints-enabled input selects spectral mode before a zero-target
setup fatal; successful constrained SCF behavior remains untested. See
``docs/validation/fixed_ldau/ENERGY_ACCOUNTING.md`` for scope and status.

For the ordinary fixed-U/J route, ``calculate_hubbard_u_potential_general``
builds the rotationally invariant Liechtenstein orbital potential from the
real, spin-resolved local occupation matrices ``ldm(l,ispin,m,m')``.  The
direct and exchange factors are generated from the supplied U and J.  For the
default ``liechtenstein`` form, the diagonal terms additionally receive the
implemented FLL-like shift

.. math::

   -U(n-\tfrac{1}{2}) + J(n_\sigma-\tfrac{1}{2}),

where :math:`n=n_\uparrow+n_\downarrow`.  This is the double-counting
convention actually used by the Hamiltonian route.  Namelist U/J values are
in eV and are converted to Ry on input.

The current route applies this potential to the on-site block of both the
real-space Hamiltonian and the shared reciprocal Fourier assembly.  Without
``hoh`` it is added to the on-site block of ``ee`` (the first-order ``h``).  With ``hoh`` it is added to
``enim``, the additive on-site term of :math:`H = h - h\,o\,h + \mathrm{enim} + \mathrm{lsham}`, so that
:math:`H = H_2 + V_U`; adding it to ``ee`` would put :math:`V_U` inside :math:`h\,o\,h` as well.
Real-space ``gbt_single_q`` follows the same rule: with ``hoh`` V_U goes into ``enim`` and the onsite term of ``h`` is ``cex``, as in the collinear builder
(``LdaUGbtHoh`` compares it with ``periodic_nc`` at q = 0; antiferromagnets are not covered, see ``tests/KNOWN_ISSUES.md``, F1).
K-space ``gbt_single_q`` takes V_U in the same way (first order: added to the onsite block of ``h``; ``hoh``: added to the onsite ``enu`` after the spin swap of a site with ``mom(3) < 0``, since V_U is built from the rotating-frame P); ``LdaUKspaceGbt`` compares it with ``periodic_nc`` at q = 0, ferromagnet and antiferromagnet.
``hoh`` with ``kspace_ham_order='first'`` stops with a fatal.  The
Hamiltonian source currently does not accumulate a separate scalar
:math:`E_U-E_{dc}` correction.  Therefore ``etot`` and ``sumev`` in a
validation run are reported as finite response observables under the active
Hamiltonian convention; they are not evidence of a fully variational,
double-counting-corrected total-energy functional.

Code mapping (developer view)
=============================

Main entry points
-----------------

- ``source/hamiltonian_hubbard.f90``: ``calculate_hubbard_u_potential_general`` builds
  the on-site Hubbard potential matrix.
- ``source/hamiltonian_hubbard.f90``: ``calculate_hubbard_v_potential`` builds the
  current proxy nearest-neighbor +V potential contribution (diagonal/local-occupation based).
- ``source/hamiltonian_build.f90``: ``build_bulkham`` and ``build_locham`` inject
  ``hubbard_u_pot`` and ``hubbard_v_pot`` into Hamiltonian blocks.
- ``source/bands.f90``: ``calculate_hubbard_u_sc`` performs the self-consistent
  ``U_eff`` update loop for channels selected by mask.
- ``source/math.f90``: ``Coulomb_mat`` and ``a_k`` provide Coulomb matrix
  elements in real spherical-harmonic form.
- ``source/math.f90``: ``tabulated_slater_integrals`` provides tabulated radial
  integral values used by SC-U flow.

Core arrays and variables
-------------------------

- ``ldm(na,l,ispin,m,m')``: local density matrix per atom type and channel.
- ``hubbard_u(l)``, ``hubbard_j(l)``: on-site U/J per ``l`` channel (internal Ry).
- ``hubbard_u_pot(i,j,na)``: assembled on-site correction in basis space.
- ``hubbard_v(i,j,li,lj)`` and ``hubbard_v_pot``: inter-site correction data (internal Ry).
  In the current implementation, ``hubbard_v_pot`` is populated by a diagonal proxy.
- ``hubbard_u_sc(itype,l)``: mask enabling SC-U update for a channel.
- ``hubbard_u_potential_form``:
  ``'liechtenstein'`` (default) or ``'acbn0'``.

Namelist inputs
===============

The relevant keys are in ``&hamiltonian``.

Units
-----

- ``hubbard_u*``, ``hubbard_j*`` and ``hubbard_v`` in ``input.nml`` are interpreted in **eV**.
- Internally, RS-LMTO-ASA stores and applies these values in **Ry**.
- Conversion is done during namelist load in ``source/hamiltonian.f90`` via division by ``ry2ev``.

1. Fixed on-site U/J (recommended baseline)
-------------------------------------------

.. code-block:: fortran

   &hamiltonian
     hubbard_u_general(1,:) = 0.0, 0.0, 2.5, 0.0
     hubbard_j_general(1,:) = 0.0, 0.0, 0.9, 0.0
     hubbard_u_potential_form = 'liechtenstein'
   /

2. Self-consistent U update (SC-U)
----------------------------------

.. code-block:: fortran

   &hamiltonian
     ! Enable SC-U on d channel for atom type 1.
     hubbard_u_sc(1,3) = 1
     ! Do not combine with explicit hubbard_u_general/hubbard_j_general.
     hubbard_u_potential_form = 'liechtenstein'
   /

3. Optional ACBN0-like potential path
-------------------------------------

.. code-block:: fortran

   &hamiltonian
     hubbard_u_sc(1,3) = 1
     hubbard_u_potential_form = 'acbn0'
   /

4. Add nearest-neighbor +V
--------------------------

.. code-block:: fortran

   &hamiltonian
     hubbard_u_general(1,:) = 0.0, 0.0, 2.5, 0.0
     hubbard_j_general(1,:) = 0.0, 0.0, 0.9, 0.0
     ! Example: V between type 1 d and type 1 d (l-index: s=1, p=2, d=3, f=4)
     hubbard_v(1,1,3,3) = 0.6
   /

Notes and constraints
=====================

- Internal energy unit is Ry; Hubbard inputs in namelists are interpreted in eV and converted on read.
- ``hubbard_u_sc`` and explicit ``hubbard_u_general/hubbard_j_general`` are mutually exclusive in the current flow.
- For ``lmax=3`` (spdf), orbital blocks are sized from the active basis dimensions (no hardwired spd-only matrix extents in the current implementation path).
- ``+V`` currently lacks the inter-site density-matrix machinery required for the
  full matrix expression. The implemented path uses nearest-neighbor, orbital-diagonal
  shifts derived from local occupations. Treat ``+V`` as a proxy/experimental correction
  until inter-site :math:`n^{JI,\sigma}_{m' m}` integration is wired from inter-site Green functions.

See `VAL-04 developer record <../../dev/VAL-04_LDA_U.md>`_ for the material response envelope,
occupation-matrix checks, convergence result, and the explicit RS/KS,
noncollinear, and SOC support classification.


Restart density-matrix files
=============================

For each atom type carrying U or J, the solver writes ``<label>_dm.nml`` beside
``<label>_out.nml``. The density file records the label, radial mode, backend,
per-channel U/J in eV, run id, iteration, and square real/imaginary blocks
indexed by angular momentum and row/column spin. Collinear output contains only
spin-diagonal blocks. The potential restart stores ``dm_file``, the name of that file next to the
potential file, and the density matrix is read from it only through that field. The
file is not checked against the potential file: keep the two from the same run. A
``dm_sha256`` field written by older versions is read and ignored; legacy inline
``ldm_flatten`` fields remain readable.

Hubbard screen output
=====================

``hubbard_log_level`` controls Hubbard diagnostics: ``summary`` (default)
prints run setup and compact per-iteration site summaries, ``detailed`` adds
per-site density, physicality, and energy diagnostics, and ``debug`` prints
the full diagnostic stream.

``hubbard_spinor_dm_diagnostic`` (default ``.false.``) extracts the full 10x10
spinor d-shell occupation ``P`` and energy-weighted ``M1`` at ``nsp>=2`` (RS and
KS) and logs ``HUBBARD_SPINOR_DM`` per site. It is a diagnostic: the Hamiltonian
does not use these matrices and fixed U with ``nsp>=2`` still stops. The spin-off-diagonal
blocks of ``P`` are written to the dm file when the array is present
(see ``docs/validation/fixed_ldau/SPINOR_DM_CONTRACT.md``).
