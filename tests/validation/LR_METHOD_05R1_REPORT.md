# LR-METHOD-05R1 — accepted Fe Ward decomposition

**BLOCKED — LOCAL BXC TO LMTO TRANSVERSE-VERTEX MAPPING NOT CLOSED**

**CORE-CONTRIBUTES-BUT-NOT-DOMINANT**

The accepted coefficient Hamiltonian satisfies the zero-broadening Ward identity to 5.57e-14. Its rigid-rotation torque differs from the compact local Bxc source by 23.45%. Core removal leaves a 23.97% compact residual. These are substantial finite discrepancies, not a proposed tolerance or a repair. A full-Hamiltonian source also exposes an additional response representation loss; the source mismatch is not claimed to explain the entire residual.

## Provenance and changes

Starting committed HEAD: `5a6ba54641a0ae22b65d16a953602ee0912c636f` on `fable_v4`. No commits were made.

Inherited 05R changes: `source/radial_ground_state.f90`, `source/linear_response.f90`, `source/linear_response_kernel_dyson.f90`, `source/linear_response_run.f90`; `docs/DECISIONS.md`, `docs/linear_response/{CONVENTIONS,FORMULATION,PUBLIC_API}.md`; `tests/unit/linear_response/kernel_dyson/{test_alsda_kernel,test_compact_interaction}.f90`; untracked `tests/validation/{lr_method_05r.py,LR_METHOD_05R_REPORT.md}`. Five inherited input modifications are `example/linear_response/{bccFe_rotation,bccFe_tddft}/input.nml` and `tests/regression/{bccFe_block,bccFe_chebyshev,bccFe_lanczos}/input.nml`. The full inherited generated-artifact inventory and patch are saved in `build-lr-method-03/Testing/05r1-evidence/{inherited-status.txt,inherited-05r.patch}`.

05R1 adds `source/lr_ward_diagnostics.f90`, `tests/validation/lr_method_05r1.py`, and this report. `source/CMakeLists.txt` registers the diagnostic module; `source/linear_response_run.f90` passes the read-only Hamiltonian to an export hook enabled only by `RSLMTO_LR_WARD_DUMP`. The exporter has input-only production objects and uses diagnostic copies for radial adapter and no-HOH comparisons.

Production Kxc remains bxc_SR/m_SR. No changes to susceptibility equations, Dyson, Hamiltonian, accepted state, physical radial bases, Jülich, Mills, rotation, RSGF or exchange.f90. No scale fit, Goldstone correction, denominator shift, reference regeneration or tolerance weakening.

The exported accepted k-space state and radial snapshot are byte-identical to the earlier 05R Fe evidence. The production compact output has exactly zero numerical difference from the existing 05R `LinearResponseCompactAlsda` output (2 rows × 15 columns).

## Accepted Fe evidence

| Quantity | Value |
|---|---:|
| EF (Ry) | -0.0622448467555 |
| Temperature (K) | 300 |
| k mesh / count | 4×4×4 / 64 |
| SCF moment (muB) | 2.00000130935 |
| Coefficient valence moment | 2.00000130935 |
| M_SR total | 2.00000720602 |
| M_P total | 2.05454931688 |
| M_P valence | 2.05455046068 |
| M_P core | -1.14380609001e-06 |
| Point norm m_P total | 0.698494252459 |
| Point norm m_P valence | 0.690749657024 |
| Point norm m_P core | 0.0107398737028 |
| Compact norm m_P total | 0.698494251887 |
| Compact norm m_P valence | 0.690749657017 |
| Compact norm m_P core | 0.0107398430056 |
| Core / valence metric norm | 0.0155481419261 |
| Rotation torque norm | 0.255666628305 |
| Compact Bxc source norm | 0.222816200405 |
| Relative torque/source mismatch | 0.234464457425 |
| Radial Pauli Bxc source norm | 0.235712309652 |
| Radial Pauli torque mismatch | 0.120243731432 |
| HOH derivative norm | 0.0328511600654 |

| eta (Ry) | Coefficient r_H | Compact total | Compact valence | Radial Lmax=0 valence | Full-H source point Pauli | Full-H source compact |
|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 0.07879167794 | 0.2513705833 | 0.2481657528 | 0.1013107208 | 0.07634352275 | 0.2354483569 |
| 0.005 | 0.03993862466 | 0.2451270195 | 0.2418311787 | 0.0844292141 | 0.05201799831 | 0.2288925558 |
| 0.0025 | 0.02004337054 | 0.2435264473 | 0.2402065274 | 0.07961058762 | 0.04380863706 | 0.2272094816 |
| 0.00125 | 0.01003115753 | 0.2431236774 | 0.2397976484 | 0.07835686915 | 0.04149852044 | 0.2267857902 |
| 0.000625 | 0.005016769697 | 0.2430228182 | 0.2396952562 | 0.0780401155 | 0.04090029841 | 0.2266796818 |
| 0.0003125 | 0.002508533932 | 0.242997593 | 0.2396696474 | 0.07796071487 | 0.04074935113 | 0.2266531431 |
| 0 | 5.56932734e-14 | 0.2429891836 | 0.2396611101 | 0.07793422903 | 0.04069890922 | 0.2266442959 |

All finite eta rows are diagnostic evaluations of the unchanged band-loop equation on one accepted state. Eta=0 is its discrete spectral limit, not a production broadening change. Exact coincident-energy terms have zero occupation difference and are omitted.

## Convention and coefficient oracle

Repository sigma+ has up/down entry 1; the circular bubble prefactor is `2 (fu-fd)/(Eu-Ed+i eta)`. With D=Hd-Hu, `<u|D|d>=(Ed-Eu)<u|d>`. The existing y-axis rotation derivative has `T_up,down=(Hu-Hd)/2=-D/2`. Thus the repository bubble acting on T gives rho_up-rho_down. The D-source coefficient oracle uses `-(fu-fd)/(Eu-Ed+i eta)` (minus one-half the repository bubble); this sign and factor follow algebraically, without fitting.

The target density matrix is constructed independently from occupied accepted eigenvectors. Residual and target norms are the square root of the k-weighted sum of squared Frobenius norms, rather than norms after averaging away k-dependent errors. All accepted bands and occupations are used.

| Deliberately incomplete/wrong source | eta-zero relative residual |
|---|---:|
| full_D | 5.56932733976e-14 |
| half_D | 0.5 |
| minus_D | 2 |
| diagonal_D | 0.117509745242 |
| d_sector_D | 0.144186856361 |
| onsite_diagonal_D | 0.135495321386 |

`diagonal_D` keeps only the orbital diagonal at each k. The additional stricter `onsite_diagonal_D` keeps the diagonal of the weighted k average, separating local diagonal from diagonal hopping. The d-sector control keeps orbital indices 4–8 only. All negative controls fail independently. A separate noncommuting Hermitian synthetic oracle also closes for full D and fails these controls.

Fixture versus accepted Hamiltonian max error: 4.443e-16; rotation commutator max error: 3.261e-16; accepted eigenpair max error: 8.607e-13. The fixture is derived by the existing validated rotation adapter; the oracle itself uses accepted `hk_bulk`.

## Local source and torque decomposition

Compact source coordinates are projected by the existing compact projection; transition coordinates are exported by the existing product transition machinery. The independent L=0 radial integration of the same six branch coefficients agrees with the compact projected source to 8.17e-16 relative error. All source/torque norms include both Hermitian spin-flip blocks: sqrt(2 sum_k wk ||T_up,down(k)||F²). No occupation weighting enters this source comparison.

| Piece | Rotation norm | Local source norm | Difference norm |
|---|---:|---:|---:|
| d_sector | 0.234076024112 | 0.210376789285 | 0.0275067563456 |
| non_d_or_cross | 0.102829177593 | 0.0734075315823 | 0.0532611479766 |
| orbital_diagonal | 0.246879970159 | 0.220302821924 | 0.0294785678191 |
| orbital_offdiagonal | 0.0664507724791 | 0.0333725308513 | 0.0521956470194 |
| R0_onsite | 0.246350083542 | 0.219986569235 | 0.0280154764814 |
| R_nonzero_hopping | 0.0683890427438 | 0.0353972953657 | 0.0529953262652 |

The d/non-d-or-cross partition and orbital diagonal/offdiagonal partition are each exhaustive. R0 onsite is the weighted k average; the remaining finite-mesh Fourier translations are grouped as hopping. This is a 4×4×4 mesh decomposition subject to translation aliasing, not a claim to resolve infinite-range bonds. The primitive cell has one site, so distinct-site offdiagonal norm is zero; translated hopping is nonzero.

The HOH contribution is the existing full derivative minus a derivative computed from a diagnostic fixture copy with HOH disabled and Enu retained. Its norm is 0.0328511600654; the no-HOH derivative norm is 0.254844761748. These derivative pieces are not orthogonal and their norms must not be added. This isolates the represented -QB term without changing the accepted Hamiltonian.

## Core and residual localization

Core signed integrated moment is -1.14380609001e-06 muB, despite a metric norm ratio of 0.0155481. Ninety percent of the absolute core moment lies within 1.85731271447 bohr. The largest positive-measure density magnitude is 0.528436975591 at 2.75304688731e-06 bohr. Core additivity max error is 3.268e-12.

The core difference is eta-independent: point metric absolute norm 0.0107398737028; relative to the valence target 0.0155481419261; Y00 infinity norm and location are recorded below. A signed integrated moment alone therefore understates the core contribution.

| eta-zero residual | Absolute metric norm | Relative norm | Point Ylm infinity | Radius of maximum (bohr) |
|---|---:|---:|---:|---:|
| Compact total | 0.169726548028 | 0.242989183618 | 1.83464145576 | 2.75304688731e-06 |
| Compact valence | 0.165545829621 | 0.239661110128 | 0.338083779854 | 0.367450459292 |
| Radial total | 0.0627081485859 | 0.0897761840776 | 1.5358681795 | 2.75304688731e-06 |
| Radial valence | 0.0538330419753 | 0.0779342290335 | 0.337392125195 | 2.75304688731e-06 |
| Core difference / valence target | 0.0107398737028 | 0.0155481419261 | 1.87326030469 | 2.75304688731e-06 |

Every eta row, absolute/relative/infinity norm, radial index and radius is retained in `diagnostics.json`, including SR-point and full-H-source actions. Compact coefficient infinity norms at eta=0 are 0.109011614223 for both targets; their reconstructed point maxima differ. Core removal changes the limiting compact relative residual from 0.242989 to 0.239661 and moves its point maximum from the innermost positive radius to 0.367450 bohr. It contributes, but is not dominant.

## Maps and response-space loss ladder

The radial coefficient of a spherical scalar is sqrt(4 pi) times its physical density/field. Point metric is sum_r W(r)|f_LM(r)|² with W=radial_weights (r²dr); all angular channels are summed when present. The origin has zero measure and is excluded from infinity-location diagnostics, without a radius floor.

The L=0 point transition is formed exactly from the existing six second-order radial branches, Gaunt factor 1/sqrt(4 pi), and up/down eigenvector overlaps within each l. Radial production uses the unchanged radial_points Pauli adapter. Compact uses the accepted complete-SR adapter, all L through 4 and dimension 348, with angular cutoff 4. Its metric is ordinary Euclidean norm in orthonormal weighted product coordinates. Magnetization is projected with the existing compact_project_magnetization service. Point reconstruction uses weighted product modes divided by sqrt(W) at positive-measure radii. Hence radial and compact each use their own mapped valence target, not an undocumented direct vector comparison.

| Full-H-source loss ladder, eta=0 | Relative residual |
|---|---:|
| Accepted coefficient density matrix | 5.56932733976e-14 |
| Pauli-adapter point, L=0 | 0.0406989092234 |
| Complete-SR-adapter point, L=0 | 0.0540925141592 |
| Complete compact, all L<=4 | 0.226644295851 |

The first loss already occurs on mapping the full represented source to the Pauli point density. The larger full compact residual includes higher angular channels and the complete-SR response adapter. The two point rows explicitly separate adapter effects; this ladder does not prove that compact truncation alone accounts for their difference. Full-H source substitution occurs only in diagnostic actions, never in production susceptibility or interaction. For the actual local field, radial valence residual is 0.0779342 while compact is 0.239661; the complete-SR L=0 local-field row is 0.115041.

## Available versus used radial ingredients

| Ingredient | Accepted snapshot available | Radial-point adapter | Used radial-point vertex | Compact adapter | Used compact vertex |
|---|---|---|---|---|---|
| phi_large | yes, pauli_large | accepted Pauli | yes | accepted Pauli | yes |
| phidot_large | yes, Pauli and SR | Pauli derivative | yes | SR derivative | yes |
| phi_small | yes, pauli_small | initialized zero | no | accepted Pauli small | no |
| phidot_small | yes, SR | initialized zero | no | SR derivative | no |
| phiddot_large | yes, SR | initialized zero | no nonzero term | SR derivative | yes |
| phiddot_small | yes, SR | initialized zero | no | SR derivative | no |
| GFAC | yes, SR | default one | no | accepted SR | no |

Both response routes use `lmto_product_second_order_radial_branch`: only large phi/dot/ddot and the adapter energy parameter enter its six polynomial branches (00,10,01,11,20,02). Radial ddot branches vanish because that adapter initializes ddot to zero. Compact ddot is nonzero. The Pauli valence magnetization service uses the squared large component linear in energy about Pauli Enu and adds frozen core separately; it does not use small, GFAC or ddot terms. Availability does not establish that unused terms should be added; no such formula is implemented here.

## Validation and reproducibility

Built `rslmto.x` and `UnitLrKernelDyson` successfully. Existing 05R focused suite: **8/8 passed** (`UnitLrAlsdaKernel`, both provenance/zero-magnetization rejection tests, `UnitLrProjectedLcmm`, `UnitLrCompactInteraction`, `UnitLrCompactDysonOracle`, `UnitLrDyson`, `UnitLrMills1U`). Synthetic independent coefficient oracle plus five negative controls passed. Accepted material diagnostic assertions passed, including full-D closure, all negative controls, accepted eigenpairs, Hamiltonian fixture, commutator and source projection checks. `git diff --check` passed.

```bash
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /Users/andersb/envs/p311/bin/python3 \
  tests/validation/lr_method_05r1.py --binary build-lr-method-03/bin/rslmto.x \
  --scratch build-lr-method-03/Testing/05r1

# Reanalyze the immutable accepted export without running SCF again:
/Users/andersb/envs/p311/bin/python3 tests/validation/lr_method_05r1.py \
  --dump build-lr-method-03/Testing/05r1/ward_inputs.bin \
  --scratch build-lr-method-03/Testing/05r1
```

Evidence: `build-lr-method-03/Testing/05r1/{ward_inputs.bin,diagnostics.json,profiles.npz,run.log,compact_alsda.dat}`. The diagnostic is intentionally restricted to this one-site collinear chi_plus Fe fixture and native little-endian binary export; it is not a general production interface. The algebra tolerance 1e-10 checks double-precision identities, not acceptable material Ward closure. No complete regression campaign was rerun or old failures relabeled.

## Decision

**Case B: BLOCKED — LOCAL BXC TO LMTO TRANSVERSE-VERTEX MAPPING NOT CLOSED.**

**CORE-CONTRIBUTES-BUT-NOT-DOMINANT.** The finite source mismatch is demonstrated on the accepted Hamiltonian; additional point/compact representation losses remain visible. Stop at this classification. No resulting repair is implemented.
