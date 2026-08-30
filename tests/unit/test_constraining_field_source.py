"""Source-level guard for the single physical constraining-field path."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ham = (ROOT / "source/hamiltonian_build.f90").read_text()
self_src = (ROOT / "source/self.f90").read_text()
exchange = (ROOT / "source/exchange.f90").read_text()


def require(text: str, needle: str) -> None:
    if needle not in text:
        raise AssertionError(f"missing source contract: {needle}")


# Ordinary RS/local and GBT onsite blocks all consume the same stored field.
require(ham, "call add_constraining_field_hmag(this, ntype)")
require(ham, "call add_constraining_field_hmag(this, this%charge%lattice%iz(nlim))")
require(ham, "call add_constraining_field_spin_block(this, it, pair_spin)")
require(ham, "this%hmag(i, i, 1, 1) = this%hmag(i, i, 1, 1) + bcon(1)")
require(ham, "block(i, i + spin_off) = block(i, i + spin_off) + bcon(1) - i_unit*bcon(2)")

# The SCF loop starts from the namelist seed, returns the computed energy, and
# stores the field on the atoms consumed by the next Hamiltonian build.
require(self_src, "this%symbolic_atom(plusbulk)%mag_cfield = this%control%constraints_bfield(:, ia)")
require(self_src, "call constrain(mom_in, mom_ref, bfield, this%lattice%nrec, etcon, diagnostics_out=diagnostics)")
require(self_src, "this%symbolic_atom(plusbulk)%mag_cfield = bfield(:, ia)")
require(self_src, "this%constraint_energy = etcon")
if "sum(this%symbolic_atom(:)%potential%etot) + this%constraint_energy" in self_src:
    raise AssertionError("constraint penalty is merged into ordinary total energy")
if "call constrain(mom_in, mom_ref, bfield, this%lattice%nrec)" in self_src:
    raise AssertionError("constraint field is still discarded at a call site")

# Exchange must consume the already-built operator rather than advance the
# controller and leave the operator stale.
if "call constrain(" in exchange:
    raise AssertionError("exchange still advances the constraint controller")

print("RESULT: PASS")
