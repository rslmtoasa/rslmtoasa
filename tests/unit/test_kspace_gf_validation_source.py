"""TDDFT-05 source contracts for the orthogonal K-space Lehmann GF."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_one_k_resolvent_is_the_shared_pure_kernel() -> None:
    kernel = (ROOT / "source/lehmann_kernel.f90").read_text(encoding="utf-8")
    provider = (ROOT / "source/tddft_chi0_green.f90").read_text(encoding="utf-8")

    assert "public :: lehmann_kspace_resolvent" in kernel
    routine = kernel.split("subroutine lehmann_kspace_resolvent", 1)[1].split(
        "end subroutine lehmann_kspace_resolvent", 1
    )[0]
    assert "allocate(" not in routine.lower()
    assert "call lehmann_kspace_resolvent" in provider


def test_orthogonal_metric_boundary_is_explicit() -> None:
    root = ROOT / "source"
    reciprocal_green = (root / "reciprocal_green.f90").read_text(encoding="utf-8")
    calculation = (root / "calculation.f90").read_text(encoding="utf-8")
    dyson = (root / "dyson_kernel.f90").read_text(encoding="utf-8")

    assert "trim(this%reciprocal_mode) /= 'ham_only'" in reciprocal_green
    assert "generalized-overlap G=(z*S-H)^(-1) with metric-consistent spectral weights" in reciprocal_green
    assert "trim(reciprocal_obj%reciprocal_mode) /= 'ham_only'" in calculation
    assert "generalized-overlap response is unsupported" in calculation
    assert "z*S - H" in dyson
    assert "This kernel hard-codes" in dyson and "S = I (the Dyson matrix" in dyson


if __name__ == "__main__":
    test_one_k_resolvent_is_the_shared_pure_kernel()
    test_orthogonal_metric_boundary_is_explicit()
    print("TDDFT-05 K-space GF source contracts PASS")
