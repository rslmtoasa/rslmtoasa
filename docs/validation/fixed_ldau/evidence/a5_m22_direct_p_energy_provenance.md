# Direct P_OUT production-energy provenance check

A supplemental evaluator check tests whether the candidate logger’s `E_Hub`
values at the converged U and spectral-J endpoints materially differ from the
same production functional evaluated on the direct `P_OUT` matrices. It uses
only preserved endpoint traces and does not alter any SCF run or decision
criterion.

The final `test_density_matrix` rows are written from
`symbolic_atom%potential%ldm` at `P_OUT` in `source/self.f90:662–670`. I
reconstructed all 50 complex spin×d×d entries from these rows and passed them
to a temporary driver calling the production `hamiltonian%evaluate_hubbard_energy`
implementation. The driver is linked against the same A5 static library used
to build the campaign executable. Inputs include local U=1.95/2/2.05,
spectral U=1.95/2/2.05, and all eight spectral J stencil endpoints.

Across the 14 endpoints, the maximum absolute difference between direct-P
production `E_Hub` and the last logged candidate `E_Hub` is
`4.163336342344337e-17 Ry`. At the three local U endpoints it is no larger
than that value; its effect on the h=.05 derivative is below `1e-15 Ry/eV`.
The corresponding spectral U and J values also agree to binary roundoff.
This directly checks that the converged diagnostic density used for candidate
energy agrees with the certified P_OUT state at the precision relevant to the
A5.3 charge-only allowance.

[`a5_m22_direct_p_energy_provenance.json`](a5_m22_direct_p_energy_provenance.json)
contains the exact hashes for the production source library, campaign
executable, temporary evaluator source/binary, every reconstructed P_OUT
input, each production-evaluator output, and the corresponding trace/log
inputs. The campaign executable SHA256 is
`52eb94c6e6f6ae12fca629a4eea1b30a6800e29736c9fdfdfd836538077f2225`; the
A5 production static library SHA256 is
`cec169fc35a4e3de32f808a11053efcf88ccff6d736daaf3dccff51546807f0f`; the
summary SHA256 is `1d3627a5f6a628531d5a0a34380ca3d48c3e66064f1516daec0f5d2f3fce6c31`.
