"""Focused source-contract test for the opt-in TDDFT calculation route."""
from pathlib import Path


def test_susceptibility_dispatch_is_registered() -> None:
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert "case ('susceptibility')" in source
    assert "call this%post_processing_susceptibility()" in source
    assert "post_processing /= 'susceptibility'" in source


def test_production_route_is_mpi_over_q_and_uses_exact_kq_service() -> None:
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert "tddft_plan = make_tddft_mpi_plan" in source
    assert "iq_start = tddft_plan%owner_q%first" in source
    assert "kq_workset = reciprocal_obj%k_workset%shifted" in source
    assert "calculate_eigenpairs_at_kpoints(kq_workset%points" in source
    assert "self_obj%refresh_xc_response_kernel()" in source
    assert "evaluate_goldstone(chi0_static%chi" in source
    assert "is_full_response = config%channel == 'full'" in source
    assert "build_four_component_chi_ks" in source
    assert "build_four_component_kernel" in source
    assert "config%chi0_backend == 'green'" in source
    assert "build_chi_ks_from_green_functions" in source
    assert "enhance_tddft_susceptibility" in source
    assert "MPI_ALLREDUCE(MPI_IN_PLACE, all_xi" in source


def test_response_reconstructs_signed_restart_site_moments_before_alsda_kernel() -> None:
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert "reciprocal_obj%eigenvalues = eigenvalues_k" in source
    assert "reciprocal_obj%eigenvectors = eigenvectors_k" in source
    assert "compute_kspace_spin_moments_spinor(reciprocal_obj, site_moments)" in source
    assert "set_site_spin_population(isite, abs(site_moments(3, isite)))" in source
    assert "set_site_spin_population(isite, sqrt(sum(site_moments(:, isite)**2)))" not in source


def test_transverse_production_keeps_ordered_circular_channels_separate() -> None:
    root = Path(__file__).resolve().parents[2]
    source = (root / "source" / "calculation.f90").read_text()
    config = (root / "source" / "tddft_config.f90").read_text()
    assert "circular_channel = 'both'" in config
    assert "left_channels_reverse" in source
    assert "right_channels_reverse" in source
    assert "_minus_plus_chi0.dat" in source
    assert "_minus_plus_legacy_dyson.dat" in source
    assert "_minus_plus_pair_dyson.dat" in source
    assert "use_qplus=.not. primary_minus_plus" in source
    assert "signed_site_populations(provider)" in (root / "source" / "tddft_goldstone.f90").read_text()


def test_static_ward_and_ground_state_provenance_are_not_dynamic_defaults() -> None:
    root = Path(__file__).resolve().parents[2]
    calculation = (root / "source" / "calculation.f90").read_text()
    chi0 = (root / "source" / "tddft_chi0.f90").read_text()
    xi = (root / "source" / "tddft_xi.f90").read_text()
    assert "build_static_chi_ks_from_eigenpairs" in calculation
    assert "build_static_direct_xi_from_operator_source" in calculation
    assert "response electron count does not match" in calculation
    assert "ground_state_response_electron_count" in calculation
    assert "real q=0 omega=0 Fermi divided difference; dynamic eta excluded" in calculation
    assert "dynamic_pair_raw_gamma_loss_peak_Ry" in calculation
    assert "tddft_static_divided_difference" in chi0
    assert "no dynamical eta" in chi0
    assert "no dynamical eta" in xi


def test_controlled_goldstone_correction_rescales_only_pair_potential_columns() -> None:
    root = Path(__file__).resolve().parents[2]
    source = (root / "source" / "calculation.f90").read_text()
    goldstone = (root / "source" / "tddft_goldstone.f90").read_text()
    assert "goldstone_mode=correct requires" in source
    assert "build_goldstone_column_correction(pair_xi_static%xi" in source
    assert "pair_operator_source%initialize(reciprocal_obj, signed_moments, config%q_points(:, iq), &" in source
    assert "pair_correction%scales" in source
    assert "pair_operators_corrected" not in source
    assert "_pair_corrected_dyson.dat" in source
    assert "k_perp_sum_rule" not in source[source.index("post_processing_susceptibility"):source.index("append_dynamic_gamma_peaks")]
    assert "dgesvd('S', 'S'" in goldstone
    assert "static Xi has material imaginary content" in goldstone
    assert "finite-eta inverse" in goldstone
    assert "write(chi0_filename" in source
    assert "trim(chi0_filename)" in source


def test_pair_potential_shadow_outputs_are_explicit_and_use_direct_xi() -> None:
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert "config%xi_backend == 'compare'" in source
    assert "lmto_pair_operator_tile_source" in source
    assert "build_direct_xi_from_operator_source" in source
    assert "call build_pair_potential_operators(" not in source
    assert "enhance_tddft_susceptibility_from_xi" in source
    assert '_legacy_dyson.dat' in source
    assert '_pair_dyson.dat' in source
    assert "pair_potential_raw_residual" in source
    assert "goldstone_correction_applied = " in source


def test_longitudinal_and_full_alsda_routes_are_explicitly_capability_gated() -> None:
    root = Path(__file__).resolve().parents[2]
    calculation = (root / "source" / "calculation.f90").read_text()
    four_component = (root / "source" / "tddft_four_component.f90").read_text()
    xc_kernel = (root / "source" / "xc_response_kernel.f90").read_text()
    assert "cartesian_transverse_kernel" in four_component
    assert "2.0_rp*circular_transverse_kernel" in xc_kernel
    assert "full_response_capability" in calculation
    assert "channel=''full'' is unavailable" in calculation
    assert "build_charge_longitudinal_channels" in calculation
    assert "build_charge_longitudinal_kernel" in calculation
    assert "channel=''longitudinal'' is unavailable" in calculation
    assert "longitudinal TDDFT remains a prototype" not in calculation
    assert "longitudinal_goldstone_constraint = none" in calculation
    assert "eigenpair-resolvent reference, not a native RS Green-function provider" in calculation


def test_goldstone_correction_input_error_is_rejected_before_response_setup() -> None:
    root = Path(__file__).resolve().parents[2]
    config = (root / "source" / "tddft_config.f90").read_text()
    calculation = (root / "source" / "calculation.f90").read_text()
    assert "goldstone_mode='correct' cannot be used with xi_backend='legacy_site_scalar'" in config
    assert "goldstone_mode='correct' requires output_xi=.true. or output_chi=.true." in config
    assert config.index("goldstone_mode='correct' cannot be used") < config.index("invalid frequency, band, temperature")
    assert "goldstone_mode=correct requires" in calculation  # retain driver defence for programmatic config construction


def test_response_mesh_resolves_inherited_fermi_after_eigenpairs() -> None:
    root = Path(__file__).resolve().parents[2]
    source = (root / "source" / "calculation.f90").read_text()
    namelist = (root / "source" / "include_codes" / "namelists" / "tddft.f90").read_text()
    assert source.index("calculate_eigenpairs_at_kpoints(reciprocal_obj%k_workset%points") < source.index(
        "calculate_canonical_band_energy(find_fermi=.true."
    )
    assert "response electron count does not match target: target=" in source
    assert "There is deliberately no &tddft EF input" in source
    assert "fermi_level" not in namelist


def test_eigenpair_baseline_boundary_rejects_unsupported_response_branches() -> None:
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert "eigenpair TDDFT baseline requires nsp=1" in source
    assert "requires magnetic_representation=periodic_nc" in source
    assert "rejects HOH/second-order Hamiltonians" in source
    assert "rejects CCOR-modified Hamiltonians" in source
    assert "rejects Hubbard-corrected Hamiltonians" in source
    assert "requires reciprocal_mode=ham_only" in source
    assert "rejects nonzero SOC" in source


def test_relativistic_boundary_is_explicit_and_precedes_spinor_response_work() -> None:
    root = Path(__file__).resolve().parents[2]
    source = (root / "source" / "calculation.f90").read_text()
    four_component = (root / "source" / "tddft_four_component.f90").read_text()
    noncollinear_guard = source.index("control_obj%is_noncollinear()")
    soc_mode_guard = source.index("control_obj%has_soc()")
    potential_soc_guard = source.index("if (has_soc) then")
    eigenpair_setup = source.index(
        "calculate_eigenpairs_at_kpoints(reciprocal_obj%k_workset%points"
    )
    assert noncollinear_guard < eigenpair_setup
    assert soc_mode_guard < eigenpair_setup
    assert potential_soc_guard < eigenpair_setup
    assert "full spinor response, torque terms, and a noncollinear kernel" in source
    assert "SOC response, anisotropy, and torque terms" in source
    assert "not a production capability" in four_component
    assert "build_four_component_kernel: full response unavailable" in four_component


def test_eigenpair_baseline_exposes_endpoint_aware_static_and_provenance_api() -> None:
    root = Path(__file__).resolve().parents[2]
    chi0 = (root / "source" / "tddft_chi0.f90").read_text()
    engine = (root / "source" / "tddft_transition_engine.f90").read_text()
    assert "build_static_chi_ks_from_eigenpairs_at_q" in chi0
    assert "accumulate_static_shifted" in engine
    assert "endpoint_provenance" in chi0
    assert "eta_role" in chi0
    assert "omega_grid_min_max_points" in chi0
    assert "q_direct = config%q_points(:, iq)" in (root / "source" / "calculation.f90").read_text()


def test_production_output_metadata_carries_complete_backend_provenance() -> None:
    root = Path(__file__).resolve().parents[2]
    calculation = (root / "source" / "calculation.f90").read_text()
    metadata = calculation[calculation.index("subroutine append_tddft_metadata") : calculation.index(
        "end subroutine append_tddft_metadata"
    )]
    assert "provenance_schema = rslmto.tddft.production.v1" in metadata
    assert "canonical_tddft_backend_name(config%chi0_backend)" in metadata
    for field in (
        "energy_integration",
        "eta_role",
        "green_eta_effective_Ry",
        "contour_points_per_segment",
        "realspace_rmax_request_Angstrom",
        "realspace_tail_tolerance",
        "response_convention",
        "source_vertex_provenance",
        "interaction_kernel_provenance",
        "goldstone_policy",
        "unsupported_feature_policy",
        "mpi_response_provenance",
    ):
        assert field in metadata


def test_pair_correction_compares_corrected_loss_to_raw_pair_loss() -> None:
    root = Path(__file__).resolve().parents[2]
    source = (root / "source" / "calculation.f90").read_text()
    goldstone = (root / "source" / "tddft_goldstone.f90").read_text()
    assert "spectral_weight_correction_is_acceptable" in source
    assert "pair_correction_preserves_raw_spectral_weight" in source
    assert "created or worsened" in source
    assert "negative spectral weight relative to the raw pair response" in source
    assert "corrected_weights < raw_weights - allowed_change*max(1.0_rp" in goldstone


if __name__ == "__main__":
    test_susceptibility_dispatch_is_registered()
    test_production_route_is_mpi_over_q_and_uses_exact_kq_service()
    test_response_reconstructs_signed_restart_site_moments_before_alsda_kernel()
    test_static_ward_and_ground_state_provenance_are_not_dynamic_defaults()
    test_controlled_goldstone_correction_rescales_only_pair_potential_columns()
    test_pair_potential_shadow_outputs_are_explicit_and_use_direct_xi()
    test_longitudinal_and_full_alsda_routes_are_explicitly_capability_gated()
    test_goldstone_correction_input_error_is_rejected_before_response_setup()
    test_response_mesh_resolves_inherited_fermi_after_eigenpairs()
    test_eigenpair_baseline_boundary_rejects_unsupported_response_branches()
    test_relativistic_boundary_is_explicit_and_precedes_spinor_response_work()
    test_eigenpair_baseline_exposes_endpoint_aware_static_and_provenance_api()
    test_production_output_metadata_carries_complete_backend_provenance()
    test_pair_correction_compares_corrected_loss_to_raw_pair_loss()
    print("RESULT: PASS")
