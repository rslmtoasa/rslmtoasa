# Log parser inventory for S0

This inventory covers Python scripts under `tests/`, `tools/`, and `docs/validation/` that refer to captured solver output and contain literal solver-log markers. The token column records markers each script matches or searches for; numeric regexes, return codes, namelist values, and data-file columns are not log tokens. Scripts that only capture a log without interpreting its contents are omitted.

| Script | Solver log tokens read |
|---|---|
| `tests/benchmarks/accp1b_physical_scf.py` | `Converged!`, `Not converged!` |
| `tests/magnetic_scf/run_magnetic_scf_regression.py` | `Converged!` |
| `tests/regression/wp9_validation/commensurate_supercell/run_wp9_commensurate_supercell.py` | `Band energy`, `fatal` |
| `tests/regression/wp9_validation/gammaH_sweep/run_wp9_gammaH_sweep.py` | `fatal` |
| `tests/run_test.py` | `fatal` |
| `tests/scf/test_si_dos_equivalence.py` | `Fermi energy` |
| `tests/validation/fixed_ldau/a3_m03_iteration_control_replay.py` | `Charge transfer at atom`, `Converged!`, `Excess charge`, `FIXED_U_RESIDUAL`, `Not converged!` |
| `tests/validation/fixed_ldau/a4_4_energy_route_repro.py` | `HUBBARD_MOMENT_ENERGY_SITE`, `HUBBARD_MOMENT_ENERGY_TOTAL`, `HUBBARD_RUN` |
| `tests/validation/fixed_ldau/a4_m21_assemble_lambda0.py` | `Converged!`, `HUBBARD_MOMENT_ENERGY_TOTAL` |
| `tests/validation/fixed_ldau/a4_m21_assembler_test.py` | `Converged!`, `HUBBARD_MOMENT_ENERGY_TOTAL` |
| `tests/validation/fixed_ldau/a4_m21_response_cache_test.py` | `Converged!`, `HUBBARD_MOMENT_ENERGY_TOTAL` |
| `tests/validation/fixed_ldau/a5_guard_boundary_cli.py` | `HUBBARD_MOMENT_ENERGY_TOTAL`, `fatal` |
| `tests/validation/fixed_ldau/a5_ks_m22_m23_runner.py` | `Converged!`, `FIXED_U_RESIDUAL`, `Not converged!` |
| `tests/validation/fixed_ldau/a5_local_w1_controls.py` | `Converged!`, `Not converged!` |
| `tests/validation/fixed_ldau/a5_spin_w1_probe.py` | `Converged!` |
| `tests/validation/fixed_ldau/a5_w2_endpoint_capture.py` | `Converged!`, `Not converged!` |
| `tests/validation/fixed_ldau/analyze_rb4_prefix.py` | `Converged!`, `Not converged!` |
| `tests/validation/fixed_ldau/c2_hubbard_operator_audit.py` | `HUBBARD_AUDIT`, `HUBBARD_AUDIT_FILE` |
| `tests/validation/fixed_ldau/c2_operator_mapping_regression.py` | `HUBBARD_AUDIT`, `HUBBARD_AUDIT_FILE` |
| `tests/validation/fixed_ldau/c2_residual_mapping_bounds_smoke.py` | `FIXED_U_RESIDUAL`, `HUBBARD_AUDIT`, `HUBBARD_MOMENT_ENERGY_SITE` |
| `tests/validation/fixed_ldau/energy_operator_oracle.py` | `Band energy` |
| `tests/validation/fixed_ldau/m02_floor_oracle.py` | `Converged!`, `Not converged!` |
| `tests/validation/fixed_ldau/rs_embedding_ra1.py` | `HUBBARD_MOMENT_ENERGY_SITE` |
| `tests/validation/fixed_ldau/rs_embedding_ra2.py` | `fatal` |
| `tests/validation/fixed_ldau/rs_embedding_ra3.py` | `fatal` |
| `tests/validation/fixed_ldau/rs_rc2_analyze.py` | `HUBBARD_MOMENT_ENERGY_TOTAL`, `Not converged!` |
| `tests/validation/fixed_ldau/rs_rc2_campaign.py` | `Not converged!` |
| `tests/validation/fixed_ldau/rs_rc2_monitored_runner.py` | `HUBBARD_MOMENT_ENERGY`, `Not converged!` |
| `tests/validation/fixed_ldau/test_a5_convergence_markers.py` | `Converged!`, `Not converged!` |
| `tests/validation/val02_reciprocal_scf.py` | `Band energy`, `Charge transfer at atom`, `Converged!`, `Not converged!`, `Spin moment of`, `Total RMS Diff`, `Total energy`, `fatal` |
| `tests/validation/val03_libxc_scf.py` | `Band energy`, `Charge transfer at atom`, `Converged!`, `Fermi energy`, `Not converged!`, `Spin moment of`, `Total RMS Diff`, `Total energy` |
| `tests/validation/val04_lda_u.py` | `Converged!`, `HUBBARD +V`, `HUBBARD_LDM_CHECK`, `HUBBARD_VDIAG`, `Not converged!`, `fatal` |
| `tests/validation/val04_lda_u_inputs.py` | `HUBBARD_RUN`, `fatal` |
| `tests/validation/val07_exchange_formulations.py` | `fatal` |
| `tests/validation/val12_lmto_fields_torques.py` | `Total energy` |
| `tests/validation/val13_abinitio_spin_dynamics.py` | `Total energy` |
| `tests/validation/val16_gbt_supercell.py` | `Band energy`, `Converged!`, `Fermi energy` |
| `tests/validation/val18_gbt_wp01_q0.py` | `fatal` |
| `tests/validation/val21_gbt_wp12_scf.py` | `fatal` |
| `tests/validation/val22_lr_radial_ground_state.py` | `Converged!`, `fatal` |
| `docs/validation/fixed_ldau/rs/C5_2_PART_D_assets/fe_rd3_b2_record/fe_local_candidate/scan_fe_i_beta005.py` | `Converged!`, `Not converged!` |
| `docs/validation/fixed_ldau/rs/C5_2_PART_D_assets/ni_rd3_ni_record/local_impurity/accepted_checkpoint_23_window4_23/scan_acceptance.py` | `Converged!`, `Not converged!` |
| `docs/validation/fixed_ldau/rs/C5_2_PART_D_assets/ni_u0_bulk_rc80_host/analysis/scan_acceptance.py` | `Converged!`, `Not converged!` |
| `docs/validation/fixed_ldau/rs/RA4_assets/ra4_analyze.py` | `Converged!`, `FIXED_U_RESIDUAL`, `Not converged!`, `Spin moment of`, `Total RMS Diff` |

`tools/run_rslmto_checked.py` parses the preflight JSON fields `status`, `preflight_stdout`, and `preflight_stderr`, and records launch `returncode`; it does not interpret solver log tokens. `docs/validation/fixed_ldau/rs/C5_4_PREFLIGHT.py` validates the deck, restart provenance, thread caps, and optional host Fermi trace; it does not parse solver log lines. The fixed-U trace readers consume named columns in `fixed_u_trace.csv`, not screen log lines. Detailed per-site energy consumers request `hubbard_log_level='detailed'`; all other listed Hubbard output markers remain available at `summary`.
