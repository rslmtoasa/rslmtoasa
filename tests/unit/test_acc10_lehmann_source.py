#!/usr/bin/env python3
"""ACC-10 source contract for the backend-owned CUDA Lehmann contraction."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_lehmann_request_crosses_the_existing_backend_seam() -> None:
    reciprocal = (ROOT / "source/reciprocal.f90").read_text(encoding="utf-8")
    backend = (ROOT / "source/reciprocal_backend.f90").read_text(encoding="utf-8")
    green = (ROOT / "source/reciprocal_green.f90").read_text(encoding="utf-8")
    cmake = (ROOT / "source/CMakeLists.txt").read_text(encoding="utf-8")

    assert "procedure(backend_contract_lehmann_if), deferred :: contract_lehmann" in reciprocal
    assert "cuda_backend_contract_lehmann" in backend
    assert "contract_lehmann(lehmann_request, lehmann_result" in green
    assert "cuda/reciprocal_lehmann.cu" in cmake


def test_canonical_arrays_and_no_residency_framework_are_preserved() -> None:
    green = (ROOT / "source/reciprocal_green.f90").read_text(encoding="utf-8")
    cuda = (ROOT / "source/cuda/reciprocal_lehmann.cu").read_text(encoding="utf-8")

    assert "green_obj%gij(:, :, :, pair) = gblk" in green
    assert "green_obj%gji(:, :, :, pair) = gblk" in green
    assert "green_obj%gij_eta" in green
    assert "cudaMalloc" in cuda and "cudaMemcpyAsync" in cuda
    assert "device_resident" not in cuda.lower()
    assert "gpu_green" not in green.lower()


if __name__ == "__main__":
    test_lehmann_request_crosses_the_existing_backend_seam()
    test_canonical_arrays_and_no_residency_framework_are_preserved()
    print("ACC-10 Lehmann source contract PASS")
