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
| `make_references_3d.py` | Generator for A4 and R1 (numpy only), with built-in self-checks; does not touch A1-A3 |
| `a4_two_channel.nml` | A4: self-consistent two-channel Stoner ferromagnet (s-like and d-like, hybridized) with an exchange matrix |
| `r1_radial_kernel.nml` | R1: radial kernel matrix K_ll' with a closed form, on the code's exponential mesh; case 2 has m = 0 on a mesh point |

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
- **A4**: chi0(0,0) Delta = -m exactly, for any mesh and temperature; with K equal
  to the exchange matrix, 1 + chi0(0,0) K has an exact zero eigenvalue with
  eigenvector m; for K = lambda I, the eigenvalue along m is exactly 1 - lambda
  (controls lambda = 2 and 1/2); the full-line zeroth moment of each channel
  equals its moment. Finite-q chi0 and chi are an independent numpy
  evaluation. Band phases are randomized.
- **R1**: K_ll' = (1/4 pi) int dr P_l P_l' K(r) / r^2 with the SR-metric product
  P_l = sqrt(gfac_up gfac_dn) G_up G_dn + S_up S_dn and K = dV/m, in closed form;
  the Simpson value on the mesh is also given. Both anchor the spherical-vertex
  (L = 0) approximation, with the core included in m.

Fortran tests read inputs, expected values and tolerances from these files
and contain no physics of their own. Changing any file here is a developer
decision, recorded in the commit message.
