# spin_response oracles

Developer-owned reference data for the `spin_response` module (B11).
Read-only to AI sessions (CLAUDE.md, "Verification integrity").

| File | Content |
|---|---|
| `make_references.py` | Generator for A1-A3 (numpy only), with built-in self-checks |
| `a1_two_level.nml` | A1: two-level model, closed form |
| `a2_single_band.nml` | A2: self-consistent single band, partially polarized |
| `a3_doubled_cell.nml` | A3: A2 on a doubled cell, random band phases |
| `tolerances.nml` | Every tolerance the spin_response tests use |
| `fe_lswt.dat` | Fe adiabatic magnon spectrum (LKAG J(r) + LSWT), meV; header gives provenance, convention and cutoff sensitivity |
| `make_fe_lswt.py` | Converts the UppASD AMS output in `lswt_inputs/` into `fe_lswt.dat` |
| `lswt_inputs/` | Raw inputs: q path, AMS output, J(r) shells (mRy) |

Regenerate with `python3 make_references.py tests/spin_response/oracles`
from the repository root. All self-checks must report `ok`. Commit the
generator, the three `.nml` files and `tolerances.nml` together.

What each model anchors (details in the generator's docstring):

- **A1**: chi0 = 5/(omega - Delta + i eta); pole at exactly delta*Delta; the
  Mills variant has residual 0.1 (mutation 11).
- **A2**: chi0(0,0) = -1/U exactly, for any mesh and temperature. Finite-q
  values are an independent numpy evaluation of the spec formula. Magnons are
  too soft to test pole finding; that is A1's job.
- **A3**: the gauge-invariant cell trace and determinant equal A2 at Q and
  Q + G_super.

Fortran tests read inputs, expected values and tolerances from these files
and contain no physics of their own. Changing any file here is a developer
decision, recorded in the commit message.
