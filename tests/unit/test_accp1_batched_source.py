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


if __name__ == "__main__":
    test_true_batched_strategy_and_layout_contract()
    test_no_automatic_dispatch_threshold_or_cpu_fallback()
    print("ACC-P1 batched source contract PASS")
