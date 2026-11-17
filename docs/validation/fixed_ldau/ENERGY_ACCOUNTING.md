# Fixed-U energy accounting evidence

> **Note (2026-11-07).** `hubbard_radial_moments='auto'` now always resolves to `spectral`. Statements below that `auto`
> selects or resolves to `local` describe the earlier rule; set `hubbard_radial_moments = 'local'` to get that behaviour.

## Current A4.4 status (KS scalar-collinear scope)

The KS M2.2/M2.3 numerical comparison passes at the explicitly limited level
“valid to the ASA base-functional accuracy measured in A5.1.” The selected
candidate is **b = F_local + E_Hub(P_OUT)**: the local-mode base free energy
plus production `evaluate_hubbard_energy` evaluated on the captured direct
`P_OUT` matrix. At U=2 eV, J=0 and at J=0.5/0.9 eV, this derivative agrees
with the fixed-P production Hubbard derivative within the preregistered
charge-only ASA allowance. Candidate c lies outside that allowance; candidate
a is report-only. No a-versus-b preference is claimed. The measured ranges
and endpoint provenance are summarized in `VALIDATION_SUMMARY.md` (raw data: tag archive/fable_v4b_u-2026-11-05).

This is a scoped numerical identity test, not a derivation of an exact general
ASA energy functional. In particular, it does not establish the radial
`SUMEV` double-count coefficient/sign for every backend or representation.
The implemented route is limited to ordinary native reciprocal KS
(`use_kspace=.true.`, `calctype='B'`) and qualifying clean-bulk RS
(`use_kspace=.false.`, `calctype='B'`), with `nbulk=0`, `nrec=ntype`, explicit
fixed general Liechtenstein U/J, scalar-collinear `nsp=1`, `lmax=2`, and
disabled constraints. SOC/noncollinear, self-consistent U, intersite V,
impurity U, impurity/surface hosts, and nonzero test fields were excluded at A4.4 (HOH is routed since 2026-11-13, see below). For
clean-bulk RS, the C4 RC2 candidate gate passes at U=2 eV and J=0, 0.5, and
0.9 eV under the amended occupation-response allowance; the independent
arithmetic is summarized in `VALIDATION_SUMMARY.md` (raw data: tag archive/fable_v4b_u-2026-11-05). The reviewed B3 extension enables the
same local correction for qualifying clean-bulk RS runs; explicit spectral
selection remains spectral/non-adopted. The measured base ASA non-stationarity
remains −0.14% in KS and −0.16% to −0.17% in RS; the RS gate does not claim
stricter stationarity or a reduction in that residual.
The scoped energy route and restart/serial-MPI checks were run; the serial route
receipts, the matched serial/MPI pair and the RS spectral scope control are
summarized in `VALIDATION_SUMMARY.md` (raw data: tag archive/fable_v4b_u-2026-11-05).
Constraints-enabled input reaches the spectral/non-adopted guard before its
zero-target setup fatal; a successful constrained SCF route remains untested.

### B3 reviewed scope extension

B3 removes the backend selector from the shared local-Hubbard eligibility
predicate. The existing guards remain: fixed Liechtenstein U/J, scalar
collinear `nsp=1`, `lmax=2`, `calctype='B'`, clean bulk (`nbulk=0`,
`nrec=ntype`), and no constraints, SOC/noncollinearity,
self-consistent U/V, impurity U, or nonzero test-operator field. Thus the
automatic local route is available to both RS and native KS when those
conditions hold. In either backend, local `E_Hub(P_OUT)` is included in the
physical total energy. Explicit spectral selection remains spectral and does
not adopt the local correction; impurity/surface and constrained cases remain
outside this eligibility scope. This describes the reviewed implementation scope; accepted-endpoint route
reproduction and unchanged-mode/no-U controls are summarized in
`VALIDATION_SUMMARY.md` (raw data: tag archive/fable_v4b_u-2026-11-05). This is not a new RS ASA
derivative-gate verdict.

### hoh (2026-11-13)

With `hoh = .true.` V_U is added to `enim` and the radial moments come from the orthonormal-basis Green function; `bands.f90` and `self.f90`
have no `hoh` branch, so the correction `c - t/q0`, `q2'` and the radial step are the ones of the first-order route. `E_Hub(P_OUT)` is routed
into `physical_total_energy` by the same predicate (`is_eligible_local_hubbard`), now with `allow_hoh=.true.`.
Finite differences, bulk bcc Fe (`tests/scf/cases/bulk/bccFe_hubbard_nsp1`, `nsp=1`, RS block, rc=80, J=0, `local`, beta 0.1, `conv_thr` 1e-9,
all runs converged, 62 to 98 steps), central differences of the converged `physical_total_energy` against `E_Hub/U` at U=2 eV
(E_Hub is linear in U at J=0, so `dE_Hub/dU` at fixed P is `E_Hub/U`); criterion fixed before the hoh runs: relative difference within 5e-3:

| run | FD, Ry/eV | E_Hub/U, Ry/eV | relative |
| --- | ---: | ---: | ---: |
| hoh, U=1.8/2.2 | 0.05342474 | 0.05355919 | -2.5e-3 |
| hoh, U=1.9/2.1 | 0.05350306 | 0.05355919 | -1.0e-3 |
| hoh, five-point U=1.8 to 2.2 | 0.05352916 | 0.05355919 | -5.6e-4 |
| first order, U=1.8/2.2 | 0.05424312 | 0.05442266 | -3.3e-3 |

`tests/ldau/test_hoh_energy_u.py` (`LdaUHohEnergy`) checks `nsp=3` against `nsp=1` with hoh (a consistency check, not an energy reference):
k-space and real-space antiferromagnet, B2FeCo impurity, Ni(001) surface; max |dE_Hub| 9.2e-13, 1.7e-13, 7.4e-12, 3.5e-10 Ry and max |dE_total|
7.4e-10, 3.0e-10, 4.5e-10, 1.4e-8 Ry. Not measured: a hoh finite difference on a cluster (the first-order impurity values are in `tests/KNOWN_ISSUES.md`),
`ccor_2c`, nsp=4.

## Gate result (historical pre-A5 status)

The paragraph below records the earlier pre-A5 assessment and is superseded
for the scoped KS path by the implementation and reproducibility results
above. It remains relevant to the broader unvalidated functional and RS scope.

The production energy-closure gate is **not established** for fixed Liechtenstein U/J in this implementation. No total-energy correction was added. The available source-connected captures validate the finite-H operator trace and the existing radial replay, but they do not derive the identity needed to place the fixed-U contribution in the current ASA radial total-energy functional.

This is an evidence limitation, not a proof that an appropriate functional cannot be constructed. The next proof must derive the energy derivative for this implementation’s radial and LMTO representations, including the actual projector metric, Q2 terms, and finite-sphere boundary contribution, or propose a direct fixed-basis kinetic-energy functional with a validated no-U limit. Any proposed change to the total-energy formula requires independent review before implementation.

## Source-derived energy path

The production `atomsc` assembly in [`source/self.f90`](../../../source/self.f90) is

```text
EKIN = SUMEV + SUMEC - RHOVH - RHOMU
UTOT = 0.5 * (RHOVH + ZVNUCL)
ETOT = EKIN + UTOT + RHOEPS
```

Thus `ETOT = SUMEV + SUMEC - 0.5*RHOVH - RHOMU + 0.5*ZVNUCL + RHOEPS`. In `NEWRHO`, `SUMEV` is accumulated as `EVAL*Q0 + Q1`; `EVAL` comes from the radial `RSEQSR` solve using `PL` and the radial potential. The scalar-relativistic implementation records Q0, Q1, Q2 moments, maps the spectral center through `gravity_center` into `PL`, then reconstructs the radial density and energy. This is the production path, not a direct sum of periodic Hamiltonian eigenvalues.

The usual DFT+U correction `E_U - Tr(V_U P)` follows when the base one-electron energy term contains the expectation value of the same Hamiltonian and density used to construct `V_U`. We have not established that identity between this code’s radial `SUMEV` and the full anisotropic, complex d-shell operator trace. In particular, this source path carries spherical per-l radial moments and radial partial waves, while the fixed-U potential acts as an orbital matrix. The observed nonzero Q2 and finite-sphere radial boundary terms are part of the radial state and cannot be omitted by assuming a direct `M1 = SUMEV` identity.

The evidence does **not** establish that an exact identity is impossible. It establishes that the simple band-trace formula cannot be inserted into `ETOT` without a derivation matching this implementation’s radial reconstruction, metric, and boundary conventions.

## Accepted-state production captures

Each capture used one SCF update from its own backend-specific accepted U=2 eV, J=0, nsp=1 converged Release seed. They are finite one-step replays, not converged material energies and not a common-H RS/KS comparison. Both runs used `nstep=1` and ended with `Not converged`, as expected for a one-step diagnostic replay. Input, output and checksum records: (raw data: tag archive/fable_v4b_u-2026-11-05).

The capture-only production `ATOMSC_ETOT` components (Ry) were:

| Component | RS | KS |
| --- | ---: | ---: |
| SUMEV | -2.5086424632182527 | -2.5059765641774083 |
| SUMEC | -1487.1142883530001 | -1487.1138967502879 |
| RHOVH | -3920.7610896942851 | -3920.7546365558183 |
| RHOMU | -142.85075975826041 | -142.84828369955358 |
| VNUCL | 234.43657027622018 | 234.43646675895775 |
| ZVNUCL | -6095.3508271817245 | -6095.3481357329019 |
| EKIN | 2573.9889186363271 | 2573.9830469409062 |
| UTOT | -5008.0559584380044 | -5008.0513861443596 |
| RHOEPS | -108.15715744611761 | -108.15536326271487 |
| ETOT | -2542.2241972477950 | -2542.2237024661681 |

These values are diagnostic outputs from separate accepted backend states. They are not evidence of a backend energy discrepancy or of the Hubbard energy derivative.

## Independent checks

The independent `energy_replay_check.py` (archive tag `archive/fable_v4b_u-2026-11-05`) reconstructs the captured production `NEWRHO` valence-density and `SUMEV` equations. On the fresh accepted-state captures, maximum valence-density residuals were `5.16e-15` (RS) and `5.83e-15` (KS); `SUMEV` reconstruction error was zero for both. Re-running the captured radial calculation in the same process gave zero residual for rho, SUMEC, SUMEV, EC, and EV. Raw data: tag archive/fable_v4b_u-2026-11-05. This is a replay of the source-connected radial path, not the Hubbard correction to `ETOT`.

The independent `radial_energy_oracle.py` (archive tag `archive/fable_v4b_u-2026-11-05`) checks the scalar-relativistic radial eigen-equation including its finite-sphere surface term. Maximum channel identity residuals were `1.4529108e-6 Ry` (RS) and `1.4525208e-6 Ry` (KS). The production charge normalization uses `GFAC*u^2 + f^2`; the eigen-equation identity uses `u^2 - f^2`, differing by roughly `2e-4` in these channels. The checker also reconstructs the captured valence density integral from Q moments within `3e-14 Ry`. Raw data: tag archive/fable_v4b_u-2026-11-05. Its explicit claim limit is that it does not identify `SUMEV` minus the density-potential integral with the complete kinetic functional, and it does not settle the Hubbard double-counting term.

The `energy_operator_oracle.py` (archive tag `archive/fable_v4b_u-2026-11-05`) was run on the fresh KS finite-H capture. It reports `Tr(V_U P) = -0.2543992705723978 Ry`, occupation count `7.999999999992392`, maximum eigensystem residual `1.14e-15`, and maximum metric-orthogonality residual `3.22e-15`. Finite-temperature free-energy slopes at lambda steps `1e-3`, `1e-4`, and `1e-5` agree with `Tr(V_U P)` to at most `1.90e-9 Ry`. Raw data: tag archive/fable_v4b_u-2026-11-05. This checks the captured finite-H operator and occupation derivative for that KS state. It does not establish that the radial `SUMEV` carries this same derivative or prove a total-energy correction.

An earlier one-step RS trace had spectral `M0=8`, `M1=-2.20667225497 Ry`, and radial `SUMEV=-2.05402849160 Ry`. Their difference was `-0.15264376337 Ry`; therefore `M1` must not be substituted for `SUMEV`. The fresh accepted-state captures also retain the actual center, PL, Q0/Q1/Q2, radial wavefunctions, and full `ATOMSC_ETOT` terms for future derivation.

Matched trace-on/off controls with identical pre-run inputs and Fe seeds produced byte-identical `Fe_out.nml` files for each backend. The controlled outputs matched exactly in ETOT/SUMEV/SUMEC. The historical one-off observer delta remains unexplained; the controlled evidence does not reproduce it. Details: `VALIDATION_SUMMARY.md` (raw data: tag archive/fable_v4b_u-2026-11-05).

## Convergence gate evidence

The fixed-U P/V residual gate is active only for fixed Liechtenstein U/J. A deliberately loose charge threshold (`conv_thr=1e6`) did not permit a one-step U=1.8 eV run to converge: its P and V residuals remained `1.1028e-1` and `1.4590e-2 Ry`, with stability streak zero. With U/J explicitly zero in both input and restart files, the no-U control retained its previous first-update convergence behavior. The first no-U attempt left U in the restart file and is excluded.

Release U=2 eV, J=0, nsp=1 production convergence and restart runs used the fixed P/V absolute tolerances `1e-7` and three consecutive stable updates. RS converged after 59 updates; its restart passed in 3. KS converged after 32 updates; its restart took 4 updates because the first update had P residual `1.0423e-7` and correctly reset the streak before the next three stable updates. These prove the stopping gate and restart stability, not the total-energy identity. Residual records: (raw data: tag archive/fable_v4b_u-2026-11-05).

Measured: the Debug unit tests (113) all passed on the source used for the captures (raw data: tag archive/fable_v4b_u-2026-11-05).

## What is still needed

To clear the production energy gate, derive and independently verify how a fixed-H operator perturbation changes the current ASA total-energy functional at a fixed production density. The derivation must connect `P_in`, `P_out`, and `P_mix`, the actual orbital projector and LMTO metric, spectral Q0/Q1/Q2 moments, radial `PL`/`EVAL`, the radial eigen-equation and surface term, and all terms in `ATOMSC_ETOT`. A radial-center/PL carrier test alone is insufficient because changing PL changes radial functions and density. A normal SCF lambda sweep is also insufficient because occupations, mixing, charge, and radial potential change together.

Primary references checked for scope and conventions:

- O. K. Andersen, “Linear methods in band theory,” *Phys. Rev. B* **12**, 3060 (1975), [doi:10.1103/PhysRevB.12.3060](https://doi.org/10.1103/PhysRevB.12.3060). This develops the linear variational basis and Hamiltonian/overlap formulation; it does not by itself establish the energy identity for this repository’s fixed-U radial implementation.
- O. K. Andersen, O. Jepsen, and D. Glötzel, “Canonical description of the band structures of metals,” in *Highlights of Condensed-Matter Theory* (Varenna, 1985), [author-hosted notes](https://www2.fkf.mpg.de/andersen/users/jepsen/Varenna_1985.pdf). The notes separate solving the one-electron problem, constructing charge density, and evaluating total energy after self-consistency; they motivate checking the actual energy path rather than assuming an eigenvalue-sum identity.
- R. W. Tank et al., “Improved LMTO-ASA methods part II: Total energy,” (1997), [author-hosted abstract](https://www2.fkf.mpg.de/andersen/docs/pub/abstract/and1997-25.html). This work addresses ASA total-energy accuracy by retaining nonspherical charge in electrostatic evaluation; it supports treating the spherical radial approximation as a relevant limitation, not a proof of a particular missing Hubbard term.
- The [Questaal ASA documentation](https://www.questaal.org/docs/code/asaoverview/) describes the ASA band/atomic cycle as exchanging Q0/Q1/Q2 moments and potential parameters. This is supporting documentation, while the findings about this implementation above come from direct source tracing and the recorded production captures.

## C5 impurity energy status

C5 Stage RD has not adopted impurity energy routing. RD2 is a measured gate
failure; its scalar-shift failure remains unresolved, and the separately
proposed state-selection change is not adopted. RD4 is NOT MEASURED: the
authorized W1 Release run exited 139 before a complete SCF iteration, so no
``dE/dλ``, ``Tr(WP)``, residual, or energy-gate verdict exists. The bounded
same-deck Debug run reached 36 loops without a bounds error; the clean8b
binary rejected the deck before SCF as input-incompatible. Neither result
establishes a fix or a pre-existing defect. Impurity Hubbard-energy routing
remains off; no correction formula or tolerance has been added or changed.
The no-U RS/KS identity checks passed; explicit-spectral identity remains
unrun, and RD5 automatic adoption and merge remain blocked by RD2. See
`VALIDATION_SUMMARY.md` (raw data: tag archive/fable_v4b_u-2026-11-05).
