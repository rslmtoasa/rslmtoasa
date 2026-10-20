# Fixed-U energy accounting evidence

## Current A4.4 status (KS scalar-collinear scope)

The KS M2.2/M2.3 numerical comparison passes at the explicitly limited level
“valid to the ASA base-functional accuracy measured in A5.1.” The selected
candidate is **b = F_local + E_Hub(P_OUT)**: the local-mode base free energy
plus production `evaluate_hubbard_energy` evaluated on the captured direct
`P_OUT` matrix. At U=2 eV, J=0 and at J=0.5/0.9 eV, this derivative agrees
with the fixed-P production Hubbard derivative within the preregistered
charge-only ASA allowance. Candidate c lies outside that allowance; candidate
a is report-only. No a-versus-b preference is claimed. The measured ranges
and endpoint provenance are in [the U result](evidence/a5_m22_u2_derivative_result.md),
[the J result](evidence/a5_m23_j_derivative_result.md), and the
[A4.3 like-for-like table](evidence/a4_3_like_for_like_observables.md).

This is a scoped numerical identity test, not a derivation of an exact general
ASA energy functional. In particular, it does not establish the radial
`SUMEV` double-count coefficient/sign for every backend or representation.
The implemented route is limited to ordinary native reciprocal KS
(`use_kspace=.true.`, `calctype='B'`) and qualifying clean-bulk RS
(`use_kspace=.false.`, `calctype='B'`), with `nbulk=0`, `nrec=ntype`, explicit
fixed general Liechtenstein U/J, scalar-collinear `nsp=1`, `lmax=2`, and
disabled constraints. SOC/noncollinear, HOH, self-consistent U, intersite V,
impurity U, impurity/surface hosts, and nonzero test fields are excluded. For
clean-bulk RS, the C4 RC2 candidate gate passes at U=2 eV and J=0, 0.5, and
0.9 eV under the amended occupation-response allowance; the independent
arithmetic is in `rs/C4_RC2_RESULTS.md`. The reviewed B3 extension enables the
same local correction for qualifying clean-bulk RS runs; explicit spectral
selection remains spectral/non-adopted. The measured base ASA non-stationarity
remains −0.14% in KS and −0.16% to −0.17% in RS; the RS gate does not claim
stricter stationarity or a reduction in that residual.
The scoped energy route and restart/serial-MPI checks pass. Individual serial route
receipts, the distinct matched serial/MPI pair, and the RS spectral scope
control are indexed in the [A4.4 control evidence](evidence/a4_4_control_evidence_index.json);
the matched pair itself is [archived here](evidence/a4_4_matched_serial_mpi_repro.json).
Constraints-enabled input reaches the spectral/non-adopted guard before its
zero-target setup fatal; a successful constrained SCF route remains untested
([paired guard receipts](evidence/a4_4_constraints_guard_paired_enabled.json),
[positive control](evidence/a4_4_constraints_guard_paired_positive.json)).

### B3 reviewed scope extension

B3 removes the backend selector from the shared local-Hubbard eligibility
predicate. The existing guards remain: fixed Liechtenstein U/J, scalar
collinear `nsp=1`, `lmax=2`, `calctype='B'`, clean bulk (`nbulk=0`,
`nrec=ntype`), and no constraints, SOC/noncollinearity, HOH,
self-consistent U/V, impurity U, or nonzero test-operator field. Thus the
automatic local route is available to both RS and native KS when those
conditions hold. In either backend, local `E_Hub(P_OUT)` is included in the
physical total energy. Explicit spectral selection remains spectral and does
not adopt the local correction; impurity/surface and constrained cases remain
outside this eligibility scope. This describes the reviewed implementation scope; accepted-endpoint route
reproduction and unchanged-mode/no-U controls are recorded in
[C4 RC3 results](rs/C4_RC3_RESULTS.md). This is not a new RS ASA
derivative-gate verdict.

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

The fresh source-connected captures were run with Debug executable SHA256 `080d803f339461237fc4da944138567b7478f2e031632f5badf2cac91056131e`, built from the remediation source with `source/self.f90` SHA256 `379d9c2e787a9e9038220efabf68f9e35bd23023b57d76ca4c754a389ec46499`. Each run used one SCF update from its own backend-specific accepted U=2 eV, J=0, nsp=1 converged Release seed. They are finite one-step replays, not converged material energies and not a common-H RS/KS comparison.

| Backend | Input SHA256 | Starting Fe_out SHA256 | Trace CSV SHA256 | Run log SHA256 | Final Fe_out SHA256 |
| --- | --- | --- | --- | --- | --- |
| RS | `2f60ebb9e14ae8047b6910ce5d173ced426decb99eb456c5e84e71f5e81ca518` | `f4482119ecb33e1bd443ef9d9481f87c0b99c18a3fea4ee7acde61d49d572e3a` | `fca17f4db149fa14f6834a7ba4452a267c303667f641bfb44915bb4d05204d27` | `b003e9290a9801c0d2921a4fd66a54e989382db99003c582aa418b34006a35c4` | `179abe64d41fe175760340c12f82ff610fc58fb895d7c06fa75f60c4ab8b8010` |
| KS | `d967bdbc17f0c5f3ea93bf35502b401ad7196145ab7af24ca3510bd55a6ec696` | `6b9ddf4d0c552302d734aa9b5c38e55e8e7d87b49c54b412d2944255cbff9f25` | `ea565449cd264ac9dcf6793adf772765367ed767a32b40d0d207f54d2bac1324` | `68278257aa35f849abd2cc7496f814eb4541b31755b249ab7fcea8f2b90269f0` | `8fff0510faec7441ba9e0eee7841f1f22d63d9a85936329efffb5ba6554df2c1` |

The corresponding captures and exact manifests are under `/tmp/rslmto-energy-trace/accepted_final_debug_rs/` and `/tmp/rslmto-energy-trace/accepted_final_debug_ks/`. Each directory contains `input.before.nml`, `Fe_out.before.nml`, `records.csv`, `manifest.json`, `run_provenance.json`, `testrun.log`, and independent checker JSON outputs. Both runs exited 0 with `nstep=1`, and each log ends with `Not converged`, as expected for a one-step diagnostic replay.

The `run_provenance.json` SHA256 values are `f040b1e6f75a16ce2c85c489c1e1c94b5333a0526ecdde95e01233c1b2dce63b` (RS) and `7ecf13b341ec1dcacb99f7154dda9aa85d3f34f85f253c0b025677a44881e5b8` (KS). Manifest SHA256 values are `2a14766acf03445670ec215daeba30ab2db15159bb0b8ae05689065028b8a78c` (RS) and `f26e26726b59fac2074faf8e1a472210f638e1cebf571faa1f2a1bf1cc5bae10` (KS). The KS manifest uses the exact captured `fermi_level` and `temperature_ry` rows; the oracle’s maximum occupation residual is `6.34e-23`.

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

The independent [`energy_replay_check.py`](../../../tests/validation/fixed_ldau/energy_replay_check.py) reconstructs the captured production `NEWRHO` valence-density and `SUMEV` equations. On the fresh accepted-state captures, maximum valence-density residuals were `5.16e-15` (RS) and `5.83e-15` (KS); `SUMEV` reconstruction error was zero for both. Re-running the captured radial calculation in the same process gave zero residual for rho, SUMEC, SUMEV, EC, and EV. Checker JSON SHA256 values are `147c5dad24a1cb00c8cc7db7756f62496ef295a3b8bc292622811654217734b4` (RS) and `bacec9ababa5929272010aae2bc17d7eff7fe79816114614f9749478dad506f3` (KS). The compact results are preserved in [`evidence/accepted_rs_replay.json`](evidence/accepted_rs_replay.json) and [`evidence/accepted_ks_replay.json`](evidence/accepted_ks_replay.json). This validates the source-connected radial replay, not the Hubbard correction to `ETOT`.

The independent [`radial_energy_oracle.py`](../../../tests/validation/fixed_ldau/radial_energy_oracle.py) checks the scalar-relativistic radial eigen-equation including its finite-sphere surface term. Maximum channel identity residuals were `1.4529108e-6 Ry` (RS) and `1.4525208e-6 Ry` (KS). The production charge normalization uses `GFAC*u^2 + f^2`; the eigen-equation identity uses `u^2 - f^2`, differing by roughly `2e-4` in these channels. The checker also reconstructs the captured valence density integral from Q moments within `3e-14 Ry`. Full outputs are preserved in [`evidence/accepted_rs_radial_oracle.json`](evidence/accepted_rs_radial_oracle.json) and [`evidence/accepted_ks_radial_oracle.json`](evidence/accepted_ks_radial_oracle.json). Its explicit claim limit is that it does not identify `SUMEV` minus the density-potential integral with the complete kinetic functional, and it does not settle the Hubbard double-counting term.

The [`energy_operator_oracle.py`](../../../tests/validation/fixed_ldau/energy_operator_oracle.py) passes the fresh KS finite-H capture. It reports `Tr(V_U P) = -0.2543992705723978 Ry`, occupation count `7.999999999992392`, maximum eigensystem residual `1.14e-15`, and maximum metric-orthogonality residual `3.22e-15`. Finite-temperature free-energy slopes at lambda steps `1e-3`, `1e-4`, and `1e-5` agree with `Tr(V_U P)` to at most `1.90e-9 Ry`. Oracle JSON SHA256 is `ebf1bd5cb6b2a6826cce733e6363968f206a9bb3b472b2373fb371d79c3f3e04`. This certifies the captured finite-H operator and occupation derivative for that KS state. It does not establish that the radial `SUMEV` carries this same derivative or prove a total-energy correction.

An earlier one-step RS trace had spectral `M0=8`, `M1=-2.20667225497 Ry`, and radial `SUMEV=-2.05402849160 Ry`. Their difference was `-0.15264376337 Ry`; therefore `M1` must not be substituted for `SUMEV`. The fresh accepted-state captures also retain the actual center, PL, Q0/Q1/Q2, radial wavefunctions, and full `ATOMSC_ETOT` terms for future derivation.

Matched trace-on/off controls with identical pre-run inputs and Fe seeds produced byte-identical `Fe_out.nml` files for each backend. The controlled outputs matched exactly in ETOT/SUMEV/SUMEC. The historical one-off observer delta remains unexplained; the controlled evidence does not reproduce it. Details are in [`TRACE_NEUTRALITY.md`](TRACE_NEUTRALITY.md).

## Convergence gate evidence

The fixed-U P/V residual gate is active only for fixed Liechtenstein U/J. A deliberately loose charge threshold (`conv_thr=1e6`) did not permit a one-step U=1.8 eV run to converge: its P and V residuals remained `1.1028e-1` and `1.4590e-2 Ry`, with stability streak zero. With U/J explicitly zero in both input and restart files, the no-U control retained its previous first-update convergence behavior. The first no-U attempt left U in the restart file and is excluded.

Release U=2 eV, J=0, nsp=1 production convergence and restart runs used the fixed P/V absolute tolerances `1e-7` and three consecutive stable updates. RS converged after 59 updates; its restart passed in 3. KS converged after 32 updates; its restart took 4 updates because the first update had P residual `1.0423e-7` and correctly reset the streak before the next three stable updates. These prove the stopping gate and restart stability, not the total-energy identity. Provenance and residual records are in the corresponding `run_provenance.json` and logs under `/tmp/rslmto-energy-trace/scf_*_release/`.

The exact-source Debug build used for the fresh captures passed 113/113 unit tests. Command: `ctest --test-dir /tmp/rslmto-ldau-remediation-debug -L unit --output-on-failure`; log `/tmp/rslmto-remediate-debug-unit-379d.log`, SHA256 `ced452343db5c7b5b1c60c7dac67d8786fbdc9ae5c1f990712710fa3fb1191a9`.

## What is still needed

To clear the production energy gate, derive and independently verify how a fixed-H operator perturbation changes the current ASA total-energy functional at a fixed production density. The derivation must connect `P_in`, `P_out`, and `P_mix`, the actual orbital projector and LMTO metric, spectral Q0/Q1/Q2 moments, radial `PL`/`EVAL`, the radial eigen-equation and surface term, and all terms in `ATOMSC_ETOT`. A radial-center/PL carrier test alone is insufficient because changing PL changes radial functions and density. A normal SCF lambda sweep is also insufficient because occupations, mixing, charge, and radial potential change together.

Primary references checked for scope and conventions:

- O. K. Andersen, “Linear methods in band theory,” *Phys. Rev. B* **12**, 3060 (1975), [doi:10.1103/PhysRevB.12.3060](https://doi.org/10.1103/PhysRevB.12.3060). This develops the linear variational basis and Hamiltonian/overlap formulation; it does not by itself establish the energy identity for this repository’s fixed-U radial implementation.
- O. K. Andersen, O. Jepsen, and D. Glötzel, “Canonical description of the band structures of metals,” in *Highlights of Condensed-Matter Theory* (Varenna, 1985), [author-hosted notes](https://www2.fkf.mpg.de/andersen/users/jepsen/Varenna_1985.pdf). The notes separate solving the one-electron problem, constructing charge density, and evaluating total energy after self-consistency; they motivate checking the actual energy path rather than assuming an eigenvalue-sum identity.
- R. W. Tank et al., “Improved LMTO-ASA methods part II: Total energy,” (1997), [author-hosted abstract](https://www2.fkf.mpg.de/andersen/docs/pub/abstract/and1997-25.html). This work addresses ASA total-energy accuracy by retaining nonspherical charge in electrostatic evaluation; it supports treating the spherical radial approximation as a relevant limitation, not a proof of a particular missing Hubbard term.
- The [Questaal ASA documentation](https://www.questaal.org/docs/code/asaoverview/) describes the ASA band/atomic cycle as exchanging Q0/Q1/Q2 moments and potential parameters. This is supporting documentation, while the findings about this implementation above come from direct source tracing and the recorded production captures.
