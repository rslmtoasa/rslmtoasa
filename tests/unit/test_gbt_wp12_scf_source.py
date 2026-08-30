#!/usr/bin/env python3
"""Static source contract for the WP12 SCF frame boundary.

The numerical GBT operator gates are covered by the earlier WP tests. This
test protects the small SCF seam that is easy to regress without changing any
of those fixed-potential oracles: constrained spirals must not feed a full
Cartesian moment mix back into the frame marker, while the diagnostic trace
must retain the old/new/mixed radial and magnetic state in one record.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def section(text: str, start: str, end: str) -> str:
    begin = text.find(start)
    finish = text.find(end, begin + len(start))
    if begin < 0 or finish < 0:
        raise AssertionError(f"could not isolate {start!r} before {end!r}")
    return text[begin:finish]


def main() -> None:
    control = (ROOT / "source/control.f90").read_text(encoding="utf-8")
    control_nml = (ROOT / "source/include_codes/namelists/control.f90").read_text(encoding="utf-8")
    self_source = (ROOT / "source/self.f90").read_text(encoding="utf-8")
    calculation_source = (ROOT / "source/calculation.f90").read_text(encoding="utf-8")
    user_docs = (ROOT / "docs/source/keywords/control_parameters.rst").read_text(encoding="utf-8")
    developer_docs = (ROOT / "docs/DEVELOPER_MAP.md").read_text(encoding="utf-8")

    assert "logical :: gbt_scf_diagnostics" in control
    assert "gbt_scf_diagnostics = this%gbt_scf_diagnostics" in control
    assert "this%gbt_scf_diagnostics = .false." in control
    assert "density_policy, gbt_scf_diagnostics" in control_nml

    mixer = section(
        self_source,
        "subroutine mix_magnetic_state(this)",
        "end subroutine mix_magnetic_state",
    )
    assert "sd_constrained_spiral" in mixer
    assert "potential%mom(:) = axis(:)" in mixer
    assert "mag_mix" in mixer
    constrained = section(mixer, "if (trim(this%control%density_policy) == sd_constrained_spiral)", "else if")
    assert "mix_magnetic_moments" not in constrained
    relaxed = section(mixer, "else if (trim(this%control%density_policy) == sd_relaxed_reference)", "else\n")
    assert "mix_magnetic_moments" in relaxed

    assert "density_policy=\"relaxed_reference\" is unsupported" in self_source
    assert "gbt_single_q" in self_source
    assert "gbt_scf_diagnostics.dat" in self_source
    trace = section(
        self_source,
        "subroutine write_gbt_scf_diagnostics(this, iteration)",
        "end subroutine write_gbt_scf_diagnostics",
    )
    for token in (
        "gbt_scf_in_moment",
        "mix%qia_old",
        "mix%qia_new",
        "mix%mag_old",
        "mix%mag_new",
        "mix%mag_mix",
        "rf_density%rho",
        "cx1",
        "wx1",
        "physical_total_energy",
        "constraint_penalty_energy",
        "constraint_field_coupling_energy",
        "constraint_metric",
    ):
        assert token in trace, token

    assert "gbt_scf_diagnostics (WP12)" in user_docs
    assert "relaxed_reference" in developer_docs
    assert "GBT SCF frame contract (WP12)" in developer_docs

    probe = section(
        calculation_source,
        "function frozen_magnon_probe_energy",
        "end function frozen_magnon_probe_energy",
    )
    assert "recip_obj%fermi_level = energy_obj%fermi" in probe
    assert "e_band = recip_obj%calculate_canonical_band_energy()" in probe
    assert "calculate_canonical_band_energy(find_fermi=.true.)" not in probe

    print("WP12 SCF source contract PASS")


if __name__ == "__main__":
    main()
