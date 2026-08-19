#!/usr/bin/env python3
"""ACC-11 source contract for the narrow CUDA eigensystem handoff."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_residency_stays_backend_owned_and_token_checked() -> None:
    header = (ROOT / "source/cuda/reciprocal_cuda.h").read_text(encoding="utf-8")
    cuda = (ROOT / "source/cuda/reciprocal_cuda.cpp").read_text(encoding="utf-8")
    lehmann = (ROOT / "source/cuda/reciprocal_lehmann.cu").read_text(encoding="utf-8")
    backend = (ROOT / "source/reciprocal_backend.f90").read_text(encoding="utf-8")
    green = (ROOT / "source/reciprocal_green.f90").read_text(encoding="utf-8")

    assert "rslmto_reciprocal_cuda_get_resident_token" in header
    assert "rslmto_reciprocal_cuda_resident_eigensystem_matches" in header
    assert "rslmto_reciprocal_cuda_contract_lehmann_resident" in header
    assert "resident_valid" in cuda and "resident_token" in cuda
    assert "rslmto_reciprocal_cuda_launch_lehmann_device" in cuda
    assert "device_eigenvectors" in lehmann
    assert "cuda_backend_contract_lehmann_resident" in backend
    assert "cuda_backend_resident_eigensystem_matches" in backend
    assert "backend%contract_lehmann_resident" in green
    assert "device_eigensystem_token" in green
    assert "cuda_backend_contract_lehmann_resident" not in green


if __name__ == "__main__":
    test_residency_stays_backend_owned_and_token_checked()
    print("ACC-11 resident source contract PASS")
