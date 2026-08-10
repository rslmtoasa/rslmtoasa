"""Focused source-contract test for the opt-in TDDFT calculation route."""
from pathlib import Path


def test_susceptibility_dispatch_is_registered() -> None:
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert "case ('susceptibility')" in source
    assert "call this%post_processing_susceptibility()" in source
    assert "post_processing /= 'susceptibility'" in source


def test_production_route_is_mpi_over_q_and_uses_exact_kq_service() -> None:
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert "region_tag='tddft-q'" in source
    assert "calculate_eigenpairs_at_kpoints(kq_points" in source
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
    assert "set_site_spin_population(isite, site_moments(3, isite))" in source
    assert "set_site_spin_population(isite, sqrt(sum(site_moments(:, isite)**2)))" not in source


def test_static_ward_and_ground_state_provenance_are_not_dynamic_defaults() -> None:
    root = Path(__file__).resolve().parents[2]
    calculation = (root / "source" / "calculation.f90").read_text()
    chi0 = (root / "source" / "tddft_chi0.f90").read_text()
    xi = (root / "source" / "tddft_xi.f90").read_text()
    assert "build_static_chi_ks_from_eigenpairs" in calculation
    assert "build_static_direct_xi_from_k_dependent_eigenpairs" in calculation
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
    assert "rescale_pair_potential_columns(pair_operators_corrected" in source
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
    assert "build_pair_potential_operators" in source
    assert "build_direct_xi_from_k_dependent_eigenpairs" in source
    assert "enhance_tddft_susceptibility_from_xi" in source
    assert '_legacy_dyson.dat' in source
    assert '_pair_dyson.dat' in source
    assert "pair_potential_raw_residual" in source
    assert "goldstone_correction_applied = " in source


def test_prototype_routes_and_full_alsda_are_explicitly_capability_gated() -> None:
    root = Path(__file__).resolve().parents[2]
    calculation = (root / "source" / "calculation.f90").read_text()
    four_component = (root / "source" / "tddft_four_component.f90").read_text()
    xc_kernel = (root / "source" / "xc_response_kernel.f90").read_text()
    assert "cartesian_transverse_kernel" in four_component
    assert "2.0_rp*circular_transverse_kernel" in xc_kernel
    assert "full_response_capability" in calculation
    assert "channel=''full'' is unavailable" in calculation
    assert "longitudinal TDDFT remains a prototype" in calculation
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
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert source.index("calculate_eigenpairs_at_kpoints(reciprocal_obj%k_points") < source.index(
        "calculate_canonical_band_energy(find_fermi=.true."
    )
    assert "response electron count does not match target: target=" in source
    assert "Remove a stale &tddft fermi_level override" in source


if __name__ == "__main__":
    test_susceptibility_dispatch_is_registered()
    test_production_route_is_mpi_over_q_and_uses_exact_kq_service()
    test_response_reconstructs_signed_restart_site_moments_before_alsda_kernel()
    test_static_ward_and_ground_state_provenance_are_not_dynamic_defaults()
    test_controlled_goldstone_correction_rescales_only_pair_potential_columns()
    test_pair_potential_shadow_outputs_are_explicit_and_use_direct_xi()
    test_prototype_routes_and_full_alsda_are_explicitly_capability_gated()
    test_goldstone_correction_input_error_is_rejected_before_response_setup()
    test_response_mesh_resolves_inherited_fermi_after_eigenpairs()
    print("RESULT: PASS")
