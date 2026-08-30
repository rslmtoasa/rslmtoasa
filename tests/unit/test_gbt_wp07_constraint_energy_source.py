"""WP07 source contracts for explicit constraint-energy semantics."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
cfd = (ROOT / "source/include_codes/abspinlib/constrain.f90").read_text()
self_src = (ROOT / "source/self.f90").read_text()
benchmark = (ROOT / "tests/benchmarks/accp1b_physical_scf.f90").read_text()
report = (ROOT / "docs/dev/GBT_REAUDIT_WP07_CONSTRAINT_ENERGY.md").read_text()


def require(text: str, needle: str) -> None:
    if needle not in text:
        raise AssertionError(f"missing source contract: {needle}")


require(cfd, "controller_penalty_energy")
require(cfd, "field_coupling_energy")
require(cfd, "lagrange_functional_represented = .false.")
require(cfd, "field_coupling = sum(bfield*mom_in)")

require(self_src, "real(rp) :: physical_total_energy")
require(self_src, "real(rp) :: constraint_penalty_energy")
require(self_src, "real(rp) :: constraint_field_coupling_energy")
require(self_src, "this%physical_total_energy = sum(this%symbolic_atom(:)%potential%etot)")
require(self_src, "this%constraint_penalty_energy = diagnostics%controller_penalty_energy")
require(self_src, "this%constraint_field_coupling_energy = diagnostics%field_coupling_energy")
require(self_src, "Physical DFT total energy of system:")
require(self_src, "Constraint penalty/controller energy (diagnostic only):")
require(self_src, "Constraint-field coupling <V_con> (band diagnostic only):")

if "sum(this%symbolic_atom(:)%potential%etot) + this%constraint_energy" in self_src:
    raise AssertionError("constraint penalty is merged into ordinary total energy")
if "sum(self_obj%symbolic_atom(:)%potential%etot) + self_obj%constraint_energy" in benchmark:
    raise AssertionError("constraint penalty is merged into benchmark total energy")
require(benchmark, "self_obj%physical_total_energy")

for phrase in (
    "current `etcon` merit value",
    "field/eigenvalue bookkeeping diagnostic",
    "no corrected-MFT field-coupling subtraction",
    "fixed-state controller oracle",
):
    require(report, phrase)

print("RESULT: PASS")
