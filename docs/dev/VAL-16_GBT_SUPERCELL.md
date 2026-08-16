# VAL-16 — GBT against commensurate magnetic supercells

## Scope

VAL-16 covers the bcc-Fe flat spirals already used by the WP9 commensurate
battery:

- `q_ss=(0,0,1/2)`: four-atom explicit supercell, 90° per (001) layer;
- `q_ss=(0,0,1/3)`: six-atom explicit supercell, 60° per layer.

The production campaign is `tests/validation/val16_gbt_supercell.py`. It
first runs the explicit supercell SCF with the current executable, requires an
SCF convergence message, and uses only the resulting `Fe*_out.nml` states for
the MFT and SCF comparisons. The old WP9 files are starting guesses only; a
GBT Hamiltonian is never paired with an old frozen potential.

The comparison reports charge, moment magnitude, Fermi energy, and band
energy per primitive atom for MFT and SCF separately. EF is a useful
diagnostic, not an energy alignment correction. No empirical energy offset is
applied.

## Diagnosis

The historical WP9 MFT rerun on the current binary, before current-state
regeneration, measured:

| case | Δcharge (e) | Δmoment (μB) | ΔEband/atom (Ry) |
|---|---:|---:|---:|
| q=1/2 | 2.53e-7 | 5.78e-4 | 7.19e-4 |
| q=1/3 | 5.71e-7 | 6.38e-3 | 5.98e-3 |

The explicit supercell directions and electron counts were correct. The GBT
and supercell EF values were also close, so this was not a charge-count or
frame-alignment failure.

As a stale-state control, the q=1/2 and q=1/3 supercell potentials were
refreshed by the current executable and then used on both sides. The residuals
fell to `6.25e-5 Ry/atom` and `2.28e-4 Ry/atom`, respectively, with moment
differences at the 1e-5–1e-4 μB level. This resolves the large historical MFT
residual as a stale-potential effect, but the force-theorem comparison still
requires the fully converged current-state campaign for its final numerical
budget.

The first remaining operator-level discrepancy is before diagonalization:
the real-space GBT run uses the primitive bcc lattice basis, while the
explicit commensurate run uses a simple-cubic supercell basis with a finite
custom cluster. Their finite-cluster pair/structure-constant sets are not the
same truncation of the infinite bcc operator. Forcing the GBT bcc case through
the VAL-01 local strux path did not change either residual. Switching the
explicit side between legacy and VAL-01 strux changed the q=1/2 residual only
from `7.19e-4` to about `5.12e-4 Ry/atom`, and left q=1/3 at about
`6.19e-3 Ry/atom`; backend agreement is therefore a check, not the GBT oracle.

This is a real-space finite-cluster representation limitation, not a scalar
energy-zero problem. The appropriate fix is to compare matched operators (or
use a periodic/k-space supercell construction), not to add an energy shift.

## Checklist

- [x] current-kernel reference-generation path added;
- [x] q=1/2 MFT comparison rerun;
- [x] q=1/3 MFT comparison rerun;
- [x] frames aligned (`mom=(0,0,1)` in the GBT rotating frame);
- [x] moments compared;
- [x] electron counts compared;
- [x] band energies compared per primitive atom;
- [x] MFT and SCF are separate campaign legs;
- [x] stale-potential hypothesis tested and large residual resolved;
- [x] remaining discrepancy traced to finite-cluster operator construction;
- [x] no empirical energy correction added.

The scoped result is: GBT/supercell equivalence is validated for the aligned
local observables after current-state regeneration, while exact band-energy
equivalence is not claimed for these unmatched finite real-space clusters.
