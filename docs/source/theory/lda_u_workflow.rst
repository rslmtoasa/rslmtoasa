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

   E_U - E_{dc} = \frac{U_{\mathrm{eff}}}{2}\sum_{\sigma}
   \mathrm{Tr}\left[n^{\sigma}(1-n^{\sigma})\right], \quad
   U_{\mathrm{eff}} = U - J.

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
  [VAL-04].  The physical validation is scoped to collinear bcc-Fe with a
  fixed d-shell U/J input; it is not a blanket validation of every magnetic,
  SOC, or self-consistent-U combination.
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

Selecting local impurity density moments does not add ``E_Hub(P_imp)`` to the
reported impurity energy; impurity energies exclude the Hubbard term. Impurity
Hubbard-energy routing remains off. The
clean-bulk eligibility and energy formula are unchanged. C5's numerical
KS-host-to-RS-impurity equivalence for Ni remains unvalidated, and RD4 has no
energy measurement: its optimized Release diagnostic has only been checked
through one complete iteration. See the C5 records in
``docs/validation/fixed_ldau/rs/C5_AMENDMENT.md`` and
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

Surface Hubbard support is density-only. Surface energies exclude the Hubbard
term; surface Hubbard-energy routing remains off. The same-material RS/KS
surface smoke and layer-moment comparison are still OPEN because this checkout
has no accepted local U-carrying Cu bulk endpoint for its only surface fixture,
fcc Cu(001). See the C6 records in
``docs/validation/fixed_ldau/rs/C6_T0_RESULTS.md`` and
``C6_AMENDMENT.md``.

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
recorded in ``docs/validation/fixed_ldau/rs/C4_RC3_RESULTS.md``; this is
route evidence, not a new ASA derivative-gate result. Explicit spectral selection remains
spectral/non-adopted. Impurity/surface hosts, SOC/noncollinear states, HOH,
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
real-space Hamiltonian and the shared reciprocal Fourier assembly.  The
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
