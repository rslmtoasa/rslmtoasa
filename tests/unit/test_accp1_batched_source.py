#!/usr/bin/env python3
"""ACC-P1 source contract for the selectable cuSOLVER batch strategy."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_true_batched_strategy_and_layout_contract() -> None:
    cuda = (ROOT / "source/cuda/reciprocal_cuda.cpp").read_text(encoding="utf-8")
    header = (ROOT / "source/cuda/reciprocal_cuda.h").read_text(encoding="utf-8")
    backend = (ROOT / "source/reciprocal_backend.f90").read_text(encoding="utf-8")

    assert "RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED" in header
    assert "rslmto_reciprocal_cuda_set_solver_strategy" in header
    assert "rslmto_reciprocal_cuda_solver_strategy_supported" in header
    assert "cusolverDnCreateSyevjInfo" in cuda
    assert "cusolverDnZheevjBatched_bufferSize" in cuda
    assert cuda.count("cusolverDnZheevjBatched(") == 1
    assert "static_cast<std::size_t>(batch_size) * sizeof(int)" in cuda
    assert "solver info[" in cuda
    assert "n <= 32" in cuda
    assert "this%solver_strategy" in backend
    assert "zheevj_batched" in backend


def test_no_automatic_dispatch_threshold_or_cpu_fallback() -> None:
    cuda = (ROOT / "source/cuda/reciprocal_cuda.cpp").read_text(encoding="utf-8")
    execute = cuda[cuda.index("extern \"C\" int rslmto_reciprocal_cuda_solve_zheevd_batch") :]

    # Strategy selection is explicit; the only size check is the API's
    # unsupported ZheevjBatched limit, not a production CPU/GPU threshold.
    assert "solver_strategy" in execute
    assert "return support > 0 ? 2 : 1" in execute
    assert "cusolverDnZheevjBatched" in execute
    assert "cusolverDnZheevd" in execute
    assert "LAPACK" not in execute


def test_explicit_fp32_strategies_and_conversion_boundary() -> None:
    cuda = (ROOT / "source/cuda/reciprocal_cuda.cpp").read_text(encoding="utf-8")
    header = (ROOT / "source/cuda/reciprocal_cuda.h").read_text(encoding="utf-8")
    backend = (ROOT / "source/reciprocal_backend.f90").read_text(encoding="utf-8")
    benchmark = (ROOT / "tests/benchmarks/accp0_real_material.f90").read_text(encoding="utf-8")

    assert "RSLMTO_RECIPROCAL_CUDA_CHEEVD_SERIAL" in header
    assert "RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED" in header
    assert "cusolverDnCheevd_bufferSize" in cuda
    assert "cusolverDnCheevd(" in cuda
    assert "cusolverDnCheevjBatched_bufferSize" in cuda
    assert "cusolverDnCheevjBatched(" in cuda
    assert "host_hamiltonians_fp32" in cuda
    assert "host_conversion_seconds" in cuda
    assert "host_widen_seconds" in cuda
    assert "fp32_cheevd" in backend
    assert "fp32_cheevj_batched" in backend
    assert "H64_H32_relative_max" in benchmark
    assert "H64_to_H32_s" in benchmark
    assert "H32_to_H64_s" in benchmark
    assert "1.0e-8_rp" in benchmark
    assert "2.0e-4_rp" not in benchmark


def test_physical_probe_is_opt_in_and_reports_scf_observables() -> None:
    probe = (ROOT / "tests/benchmarks/accp1b_physical_scf.f90").read_text(encoding="utf-8")
    driver = (ROOT / "tests/benchmarks/accp1b_physical_scf.py").read_text(encoding="utf-8")
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")

    assert "self_obj%use_kspace = .true." in probe
    assert "self_obj%reciprocal_scf_cache%make_execution_backend" in probe
    assert "ACCP1B_SCF" in probe
    assert "electron_count" in probe
    assert "site_moment" in probe
    assert "near_ef_min_abs" in probe
    assert "near_ef_value_1" in probe
    assert "near_ef_value_4" in driver
    assert "comparisons_to_fp64_gpu" in driver
    assert "ReciprocalAccP1bPhysicalSCF" in cmake


if __name__ == "__main__":
    test_true_batched_strategy_and_layout_contract()
    test_no_automatic_dispatch_threshold_or_cpu_fallback()
    test_explicit_fp32_strategies_and_conversion_boundary()
    test_physical_probe_is_opt_in_and_reports_scf_observables()
    print("ACC-P1 batched source contract PASS")
