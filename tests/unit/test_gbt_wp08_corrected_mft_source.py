"""Source-contract checks for the corrected constrained frozen-magnon mode.

The expensive Fe/Ni/FeCo validation campaign is documented separately.  These
checks keep the architectural boundary executable in every build: corrected MFT
must snapshot/restore the ordinary potential, iterate only the constraint field,
and report the field/energy bookkeeping without folding controller diagnostics
into the physical reference energy.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CALC = (ROOT / "source/calculation.f90").read_text()
SELF = (ROOT / "source/self.f90").read_text()
FM = (ROOT / "source/frozen_magnon.f90").read_text()
FM_NML = (ROOT / "source/include_codes/namelists/frozen_magnon.f90").read_text()


def require(text: str, needle: str) -> None:
    assert needle in text, f"missing source contract: {needle!r}"


def between(text: str, start: str, end: str) -> str:
    left = text.index(start)
    right = text.index(end, left)
    return text[left:right]


def test_mode_boundary_and_configuration() -> None:
    require(FM, "character(len=24) :: mode")
    require(FM, "if (this%mode == 'corrected_mft') this%mode = 'mft_constrained'")
    require(FM, "this%mode /= 'mft_constrained'")
    require(FM_NML, "constraint_max_iterations")
    require(FM_NML, "constraint_start_from_zero")


def test_fixed_potential_step_has_no_scf_update() -> None:
    step = between(SELF, "subroutine run_fixed_potential_constraint_step", "end subroutine run_fixed_potential_constraint_step")
    require(step, "frozen_potential(ia) = this%symbolic_atom(ia)%potential")
    require(step, "call run_recursion(this, iteration)")
    require(step, "call run_dos(this)")
    require(step, "this%symbolic_atom(ia)%potential = frozen_potential(ia)")
    assert "call run_scf(this)" not in step
    assert "mixpq" not in step


def test_corrected_branch_resets_q_and_uses_raw_field_included_band_delta() -> None:
    branch = between(
        CALC,
        "else if (corrected_mft) then",
        "else\n         ! scf: fully self-consistent spiral",
    )
    require(branch, "reset_constraint_for_fixed_potential")
    require(branch, "run_fixed_potential_constraint_step")
    require(branch, "constraint_max_iterations")
    require(branch, "eband_q(iq) = frozen_magnon_probe_energy")
    require(CALC, "omega_q(iq) = 4.0_rp*(eband_q(iq) - eband_q(1))")
    require(branch, "etot_q(iq) = etot_ref")
    assert "eband_q(iq) + constraint_penalty_q(iq)" not in CALC
    assert "eband_q(iq) + constraint_coupling_q(iq)" not in CALC


def test_reference_checksum_and_zero_field_limit() -> None:
    require(SELF, "fixed_potential_checksum")
    require(SELF, "fixed_potential_max_drift")
    require(SELF, "fixed_potential_transient_drift")
    require(SELF, "this%fixed_potential_max_drift = max_drift")
    require(CALC, "bare_mft .or. (corrected_mft .and. .not. constraints_enabled)")
    require(CALC, "mtot_q(:, iq) = mtot_q(:, 1)")


def test_q_order_and_mode_distinction_are_explicit() -> None:
    require(CALC, "reset_constraint_for_fixed_potential(fm_obj%constraint_start_from_zero)")
    require(CALC, "else\n         ! scf: fully self-consistent spiral")
    require(CALC, "mode=''mft_constrained'' is currently supported for branch_mode=''acoustic'' only.")


def test_machine_readable_audit_fields_are_present() -> None:
    for needle in (
        "ordinary_potential_frozen = T",
        "constraint_field_converged = per_q_column",
        "gauge_reference_used = F",
        "field_coupling_and_controller_penalty_are_diagnostics_only = T",
        "ordinary_potential_checksum_L1",
        "potential_max_drift",
        "DeltaE_final",
        "Bmag_1 .. Bmag_nrec",
    ):
        require(CALC, needle)
