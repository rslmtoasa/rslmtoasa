# Mn2F2: frozen-potential screened-constant reciprocity intervention

Date: 2026-10-21. Stacked on `codex/mn2f2-gf-symmetry-audit` at
`2ef44c63fc64cc23cc3a914866627f8588f299f4` (PR #10).

## Status and scope

The completed first-pair gate supports the causal interpretation: restoring
real-basis Sbar reciprocity restores the tested Hamiltonian's Hermiticity and
reduces time-odd Dx from `3.8867415e-7` to `-2.455e-15` mRy on that bond.
The full 24-pair MPI endpoint runs are in progress. Full-run tables and the final
conclusion will be added after their manifests and complete dumps are verified.
No incomplete, interrupted or file-replacement-only run is endpoint evidence.

This is a fixed-potential intervention, not an SCF calculation or a study of
screening convergence. Screening cutoffs, real-space sums, HoH, SOC, recursion,
energy grid, atom files and exchange formulas are unchanged. Neither Dx nor any
GF component is set to zero. Diagonal domains and lambda interpolation have not
been run in this task.

## Intervention and binary validation

In the real Cartesian orbital basis stored in the binary, each directed block
is matched by source type, destination type and displacement to its unique reverse:

$$X_{ij}(R)=\frac{\bar S_{ij}(R)+\bar S_{ji}(-R)^T}{2},\qquad
\bar S^{\rm rec}_{ij}(R)=X_{ij}(R),\quad
\bar S^{\rm rec}_{ji}(-R)=X_{ij}(R)^T.$$

Only a genuinely onsite self-reversed key is averaged with its own transpose.
There is no spherical-orbital C operation and no offsite self-symmetrization.
The original binary is read-only; `sbar.reciprocal` is a separate exclusive-create
file. Four-byte little-endian record markers, nine float64 entries per record,
block order and counts are retained. Rounded preview values are never used as
matrix data. Vector labels identify physical bonds with absolute tolerance
`2e-6` in their Cartesian units; matching must be unique.

| Validation over all 86 directed blocks | Original | Reciprocal |
|---|---:|---:|
| max element of Sij-Sji^T | 2.1474105637e-5 | 0 exactly |
| concatenated Frobenius norm | 1.9715686202e-4 | 0 exactly |
| missing/ambiguous reverse partners | 0 | 0 |

There are 41 distinct offsite reverse pairs and four self-reversed onsite keys
(45 unique pairs including onsite). Maximum `|delta S|=1.0737052818e-5` and
concatenated `||delta S||F=9.8578431010e-5`. Parse/encode roundtrip preserves all
bytes; applying the projection twice gives identical bytes. The two axial jobs
give the same projected binary hash.

Original SHA256:
`951f587460c1ade64d6676e1a6ab591f3bd7de831fe5ff54c716cfa25f818c4a`.
Reciprocal SHA256:
`ab18c763931f6f794b06abd23a3c4b96c959e46c81a3bc8a47a21fe30cc509eb`.

## Why an explicit diagnostic input hook is necessary

The actual exchange entry point calls `structb(.true.)`, which recomputes Sbar
in memory. It does **not** read saved `sbar.`. Replacing that file alone left the
captured Hamiltonian unchanged in an initial negative control. Consequently,
such a run cannot test the proposed intervention.

The optional `RSLMTO_AUDIT_SBAR` hook loads the projected copy after the unchanged
screening call and before Hamiltonian assembly. It requires GF audit mode and
`RSLMTO_AUDIT_SBAR_VIEW`, checks freshly generated versus saved vector/type/order
labels, checks record counts/finite entries, and announces the load. An unset
variable leaves the ordinary path intact. The runner now rejects a successful
process that did not announce the override. All copied file hashes are checked
before execution: only `sbar.` differs in the intervention; full-run pair lists
are not edited. Saved output-potential hashes are also checked after execution.

## Hamiltonian gate

The input intervention has a direct operator interpretation:

$$h_{ij}-h_{ji}^{\dagger}
=W_i[1_2\otimes(\bar S_{ij}^{T}-\bar S_{ji})]W_j.$$

The actual full finite neighbour graph contains 34,400 directed blocks. Its
read-only check compares each bare-hopping edge with its unique reverse and
checks the weighted hopping used in recursion against `hO=h*O`. It excludes
allocation padding beyond the actual 1600 atoms.

| Full-graph check, [100] | Original | Reciprocal |
|---|---:|---:|
| max bare-h Hermiticity defect (Ry) | 4.6260833085e-7 | 2.7755575616e-17 |
| full bare-h Frobenius defect (Ry) | 1.6038935243e-4 | 4.4961696711e-15 |
| max hO identity residual (dimensionless) | 1.1102230246e-16 | 1.1102230246e-16 |
| O / energy-centre / SOC Hermiticity defects | 0 / 0 / 0 | 0 / 0 / 0 |

These exhaustive graph checks establish the Hermitian ingredients of

$$H=h-hOh+E_\nu+H_{\rm SOC}.$$

Actual first-action captures independently check the effective Hamiltonian,
including the two-hop term, without replacing H by its Hermitian part:

| First bond r/alat=(-1/2,+1/2,0), [100] | Original max (Ry) | Reciprocal max (Ry) |
|---|---:|---:|
| bare h pair | 2.8921200027e-8 | 9.1940344227e-17 |
| HoH pair | 5.5173235879e-8 | 8.3266726847e-17 |
| effective H pair | 5.5026866503e-8 | 9.3675074565e-17 |
| effective H onsite, maximum of i/j | 4.2646425377e-8 | 4.1733285233e-17 |
| direct mz*T covariance over captured rows | 1.2035797359e-10 | 1.2035838992e-10 |

For the reversed domain, the reciprocal effective-pair and onsite maxima are
`8.3266726847e-17` and `5.5511151231e-17` Ry. Direct mz*T covariance remains
`1.2045175968e-10` Ry. Bare-h spatial covariance over all 86 type/displacement
blocks stays controlled: the largest tested direct residual remains approximately
`2.34e-10` Ry. Restoring Hermiticity does not erase the saved inputs' small
magnetic/crystal symmetry differences.

The graph check certifies all bare-hopping ingredients; the action dumps measure
all output rows of the selected effective-H columns. No dense 28800-by-28800
effective Hamiltonian was constructed. For an intervention, a full-graph ingredient
failure above `1e-12` aborts before recursion.

## Completed first-pair Green-function and exchange gate

The previous audit derives mirror signs and the reciprocal orbital transpose
`R(A)=C A^T C`. Its sector definitions and all exchange contractions are retained.
These are maxima of `||M G M - s G||F` across the 64 sampled energies, ij direction:

| Sector, [100] first pair | Original | Reciprocal |
|---|---:|---:|
| G01 | 3.7688176659e-6 | 4.7565269240e-10 |
| Gx1 | 3.6717589946e-6 | 5.0783074214e-10 |
| Gx0 control | 6.3844985739e-8 | 1.9802792208e-8 |
| Gz1 control | 3.0725215504e-8 | 1.5907801344e-9 |

For a sector of mirror sign s, the forbidden-block norm is half this residual.
The small remaining leakage should not be called exact symmetry restoration.

Using the **same ordered bond** in reversed domains, in mRy:

| Quantity | Original | Reciprocal |
|---|---:|---:|
| Dx([100]) | +3.886604e-7 | -1.378187e-11 |
| Dx([-100]) | -3.886879e-7 | -1.377696e-11 |
| Dx_odd = (Dx+ - Dx-)/2 | +3.8867415e-7 | -2.455e-15 |
| Dx_even = (Dx+ + Dx-)/2 | -1.375e-11 | -1.3779415e-11 |
| Dz([100]) | -1.700948e-3 | -1.701022e-3 |
| Dz([-100]) | -1.700950e-3 | -1.701023e-3 |

The odd signal falls by about eight orders of magnitude on this bond, while the
small even remainder persists. Dz is finite and changes under the intervention;
it must not be claimed to be unaffected. Printed-output precision limits the last
digits of the parity differences. Full-run CCx/SCx orbital-block decomposition,
all-bond/shell comparisons and summed even/odd tables are pending the complete runs.

## Validation and reproducibility

The 21-test suite passes, including shuffled reverse pairs, onsite, Mn-F and F-F,
missing/ambiguous partners, format errors, exclusive writes, marker preservation,
byte idempotence, actual Fortran input loading/rejection and graph checks.
The GNU Fortran 13.3 build passes. Production exchange, Green-function,
recursion, Hamiltonian and screening formulas have no source diff.

See [diagnostic run instructions](diagnostics/gf_symmetry/README.md#frozen-input-reciprocity-intervention).
The full endpoints use a single MPI executable, two ranks, four OpenMP threads
per rank and one OpenBLAS thread, with original 2000-channel grid, `lld=250`,
`hoh=T`, `local_axis=F`, geometry, pair list and saved `*_out.nml` potentials.
Raw binaries, potentials, matrices and detailed manifests remain private.

## Causal interpretation and limits

The first-pair intervention changes the independently identified upstream
reciprocity defect, restores Hermitian Hamiltonian ingredients, collapses the
forbidden G01/Gx1 blocks and removes the dominant time-odd Dx. This supports the
nonreciprocity mechanism rather than a physical higher-order interpretation of
that signal. Full-run confirmation is pending.

It does not identify why the original screening construction produces the
nonreciprocal coefficients. No cutoff, radius, real-space sum or convergence
parameter was varied, and no production correction is proposed. This fixed-input
counterfactual is not a prediction from a newly self-consistent electronic state.

## Source references

| Role | Source |
|---|---|
| Physical reverse-pair projection and non-overwriting writer | [reciprocal_sbar.py](diagnostics/gf_symmetry/reciprocal_sbar.py#L18) |
| Binary reader and type/displacement identification | [sbar_audit.py](diagnostics/gf_symmetry/sbar_audit.py#L22) |
| Saved-input loader with vector-order validation | [gf_audit_io.f90](source/gf_audit_io.f90#L72) |
| Exhaustive actual bare-hopping graph and weighted-h check | [gf_audit_io.f90](source/gf_audit_io.f90#L11) |
| Exchange-only load after unchanged structb call | [calculation.f90](source/calculation.f90#L816) |
| Original screening construction | [lattice.f90](source/lattice.f90#L2425) |
| Bare noncollinear hopping assembly | [hamiltonian.f90](source/hamiltonian.f90#L2689) |
| Actual HoH action, unchanged | [recursion.f90](source/recursion.f90#L2317) |
| GF reciprocal split, unchanged | [green.f90](source/green.f90#L386) |
| Physical DMI contractions/integration, unchanged | [exchange.f90](source/exchange.f90#L1176) |
| Fresh-copy runner and hashes | [run_saved.py](diagnostics/gf_symmetry/run_saved.py) |
| Complete-run parity and trace analysis | [causal_summary.py](diagnostics/gf_symmetry/causal_summary.py) |
