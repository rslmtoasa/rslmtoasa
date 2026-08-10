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


def test_response_reconstructs_restart_site_moments_before_alsda_kernel() -> None:
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert "reciprocal_obj%eigenvalues = eigenvalues_k" in source
    assert "reciprocal_obj%eigenvectors = eigenvectors_k" in source
    assert "compute_kspace_spin_moments_spinor(reciprocal_obj, site_moments)" in source
    assert "set_site_spin_population(isite, sqrt(sum(site_moments(:, isite)**2)))" in source


def test_sum_rule_is_not_injected_as_a_finite_eta_dyson_kernel_and_manifest_is_well_formed() -> None:
    source = (Path(__file__).resolve().parents[2] / "source" / "calculation.f90").read_text()
    assert "kernel(isite, isite) = goldstone_result%k_perp(isite)" in source
    assert "kernel(isite, isite) = goldstone_result%k_perp_sum_rule(isite)" not in source
    assert "sum_rule_dynamic_kernel_applied" in source
    assert "write(chi0_filename" in source
    assert "trim(chi0_filename)" in source


if __name__ == "__main__":
    test_susceptibility_dispatch_is_registered()
    test_production_route_is_mpi_over_q_and_uses_exact_kq_service()
    test_response_reconstructs_restart_site_moments_before_alsda_kernel()
    test_sum_rule_is_not_injected_as_a_finite_eta_dyson_kernel_and_manifest_is_well_formed()
    print("RESULT: PASS")
