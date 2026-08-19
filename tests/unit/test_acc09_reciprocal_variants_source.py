#!/usr/bin/env python3
"""ACC-09 source contract for validated reciprocal CUDA operator variants."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_cuda_variant_boundary_is_host_assembled_and_standard_hermitian() -> None:
    backend = (ROOT / "source/reciprocal_backend.f90").read_text(encoding="utf-8")
    fourier = (ROOT / "source/reciprocal_fourier.f90").read_text(encoding="utf-8")
    cuda = (ROOT / "source/cuda/reciprocal_cuda.cpp").read_text(encoding="utf-8")

    capabilities = backend[backend.index("module subroutine cuda_backend_capabilities") :
                           backend.index("end subroutine cuda_backend_capabilities")]
    assert "capabilities%standard_hermitian = this%initialized" in capabilities
    assert "capabilities%generalized_hermitian = .false." in capabilities
    assert "capabilities%first_order_assembly = .false." in capabilities
    assert "capabilities%second_order_assembly = .false." in capabilities
    assert "CCOR" in capabilities and "HOH" in capabilities

    normal_mesh = fourier[fourier.index("module subroutine execute_normal_mesh_tiles") :
                         fourier.index("end subroutine execute_normal_mesh_tiles")]
    assert "host_assembler%assemble_batch" in normal_mesh
    assert "allocate(request%input_hamiltonian" in normal_mesh
    assert "request%assemble_hamiltonian = assemble_on_backend" in normal_mesh
    assert "generalized .and. .not. backend_caps%generalized_hermitian" in normal_mesh

    # The CUDA translation unit is an eigensolver/lifecycle adapter.  Operator
    # construction remains in the established CPU reciprocal assembler.
    for forbidden in ("eeo", "eecc", "ccor_2c", "hoh"):
        assert forbidden not in cuda.lower()


def test_acc09_does_not_claim_generalized_cuda_support() -> None:
    backend = (ROOT / "source/reciprocal_backend.f90").read_text(encoding="utf-8")
    execute = backend[backend.index("module subroutine cuda_backend_execute_batch") :
                      backend.index("end subroutine cuda_backend_execute_batch")]
    assert "request%generalized" in execute
    assert "allocated(request%input_overlap)" in execute
    assert "standard-Hermitian eigensolution of host-assembled H(k) tiles" in execute


if __name__ == "__main__":
    test_cuda_variant_boundary_is_host_assembled_and_standard_hermitian()
    test_acc09_does_not_claim_generalized_cuda_support()
    print("ACC-09 reciprocal variant source contract PASS")
