"""WP06 source contracts for frame-aware constraint state and diagnostics."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
cfd = (ROOT / "source/include_codes/abspinlib/constrain.f90").read_text()
self_src = (ROOT / "source/self.f90").read_text()
ham = (ROOT / "source/hamiltonian_build.f90").read_text()
bands = (ROOT / "source/bands.f90").read_text()
reciprocal = (ROOT / "source/reciprocal_projection.f90").read_text()


def require(text: str, needle: str) -> None:
    if needle not in text:
        raise AssertionError(f"missing WP06 source contract: {needle}")


for field in (
    "target(:, :)",
    "actual(:, :)",
    "moment_magnitude(:)",
    "angular_error(:)",
    "transverse_residual(:, :)",
    "bfield(:, :)",
    "bfield_longitudinal(:)",
    "torque(:, :)",
):
    require(cfd, field)

require(self_src, "diagnostics_out=diagnostics")
require(self_src, "this%constraint_metric = diagnostics%convergence_metric")
require(self_src, "this%constraint_converged = this%constraint_metric <= this%constraint_tolerance")
require(self_src, "this%symbolic_atom(plusbulk)%constraint_actual = mom_in(:, ia)")
require(self_src, "this%symbolic_atom(plusbulk)%potential%mom(:) = [0.0_rp, 0.0_rp, marker_sign]")
require(self_src, "constraints_mom_ref must have shape (3,nrec)")
require(self_src, "constraints_bfield must have shape (3,nrec)")

require(ham, "The term is onsite and carries no")
require(ham, "call add_constraining_field_spin_block(this, it, pair_spin)")
require(bands, "symbolic_atoms(plusbulk)%constraint_target(:)")
require(reciprocal, "constraint_target(:)")

print("RESULT: PASS")
