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


def test_sum_rule_is_not_injected_as_a_finite_eta_dyson_kernel_and_manifest_is_well_formed() -> None:
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert "kernel(isite, isite) = goldstone_result%k_perp(isite)" in source
    assert "kernel(isite, isite) = goldstone_result%k_perp_sum_rule(isite)" not in source
    assert "sum_rule_dynamic_kernel_applied" in source
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


if __name__ == "__main__":
    test_susceptibility_dispatch_is_registered()
    test_production_route_is_mpi_over_q_and_uses_exact_kq_service()
    test_response_reconstructs_signed_restart_site_moments_before_alsda_kernel()
    test_static_ward_and_ground_state_provenance_are_not_dynamic_defaults()
    test_sum_rule_is_not_injected_as_a_finite_eta_dyson_kernel_and_manifest_is_well_formed()
    test_pair_potential_shadow_outputs_are_explicit_and_use_direct_xi()
    print("RESULT: PASS")
