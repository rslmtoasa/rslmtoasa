#!/usr/bin/env python3
"""ACC-04 source contract for the host-assembled arbitrary-k CUDA bridge."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source" / "reciprocal_fourier.f90"


def test_arbitrary_k_uses_capability_driven_host_handoff() -> None:
    source = SOURCE.read_text(encoding="utf-8")
    start = source.index("module subroutine calculate_eigenpairs_at_kpoints")
    end = source.index("end subroutine calculate_eigenpairs_at_kpoints", start)
    routine = source[start:end]

    assert "execution_backend%capabilities(backend_caps)" in routine
    assert "backend_caps%second_order_assembly" in routine
    assert "backend_caps%first_order_assembly" in routine
    assert "host_assembler%assemble_batch" in routine
    assert "request%input_hamiltonian(nmat,nmat,tile_length)" in routine
    assert "request%assemble_hamiltonian = assemble_on_backend" in routine
    assert "request%assemble_overlap = assemble_on_backend .and. use_generalized" in routine
    assert "call this%execution_backend%synchronize()" in routine
    assert "result%operator_generation /= request%operator_generation" in routine
    assert "backend name" not in routine.lower()


def test_cuda_arbitrary_k_integration_is_registered() -> None:
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    assert "ReciprocalCudaArbitraryK" in cmake
    assert "SKIP_RETURN_CODE 77" in cmake


if __name__ == "__main__":
    test_arbitrary_k_uses_capability_driven_host_handoff()
    test_cuda_arbitrary_k_integration_is_registered()
