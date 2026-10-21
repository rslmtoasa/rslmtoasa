# Mn2F2: frozen-potential screened-constant reciprocity intervention

Date: 2026-10-21. Stacked on `codex/mn2f2-gf-symmetry-audit` at
`2ef44c63fc64cc23cc3a914866627f8588f299f4` (PR #10).

## Status and scope

The completed first-pair gate supports the causal interpretation: restoring
real-basis Sbar reciprocity restores the tested Hamiltonian's Hermiticity and
reduces time-odd Dx from `3.8867415e-7` to `-2.455e-15` mRy on that bond.
The reciprocal [100] 24-pair MPI endpoint completed and was analyzed successfully.
The execution service subsequently disconnected and could not be recovered.
The full reversed-domain endpoint and same-executable original baselines remain
unverified; this PR is intentionally draft and the full-domain test is incomplete.
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
These are maxima of `||M G M - s G||F` in inverse Ry across the 64 sampled energies, ij direction:

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

## Confirmed full reciprocal [100] endpoint

All 24 original pairs completed with return code 0 in 698.28 s. Both MPI ranks
announced loading all 86 projected blocks, the original pair-list input hash
remained `c0ad589c764f88148c0a832d89539b5be60e4d2b84f42389c1745543deb1b583`,
and all saved output-potential hashes were unchanged. The complete GF and
H-action dumps were parsed and the spatial partner maps validated before the
execution connection was lost.

The calculation executable SHA256 is
`3ba06baf0852f56449e4826e263b06e5954a343aa19b5439982c39f29ed890ed`
(GNU Fortran 13.3, OpenBLAS 0.3.26, Open MPI 4.1.6; two ranks,
four OpenMP threads per rank, BLAS one thread).

For the 12 Mn1-to-Mn2 bonds, in mRy:

| Sum | Reciprocal [100] |
|---|---:|
| Dx | -4.3910e-16 |
| Dy | -3.8661e-15 |
| Dz | 0 |
| CCx | -7.0730e-16 |
| SCx | +2.6500e-16 |

The supplied original [100] benchmark has sum Dx `+1.82060157e-6` mRy.
This supplied benchmark is **not** presented as a completed same-executable
baseline. Small discrepancies between the independently printed D and CC+SC
sums are at their rounding level.

Across all 24 H-action captures, maximum element residuals are:

| Quantity | Reciprocal [100] maximum (Ry) |
|---|---:|
| bare-h pair Hermiticity | 1.8962018370e-16 |
| HoH pair Hermiticity | 1.1796130168e-16 |
| effective-H pair Hermiticity | 1.8471416906e-16 |
| effective-H onsite Hermiticity | 5.6182374711e-17 |
| direct mz*T covariance | 1.2035838992e-10 |

Across all 24 pairs, 64 sampled energies and both directions, mirror-residual
Frobenius maxima, in inverse Ry, are:

| Sector | Reciprocal [100] mirror residual |
|---|---:|
| G00 | 3.4928303045e-8 |
| G01 | 1.2467169606e-9 |
| Gx0 | 3.8552353477e-8 |
| Gx1 | 1.3207890962e-9 |
| Gz0 | 5.9193153088e-10 |
| Gz1 | 1.5907809497e-9 |

The forbidden-block norm is half the corresponding mirror residual. Residual
forbidden blocks remain; claiming they vanish exactly would be incorrect.

Compared with the **supplied** original [100] output, maximum bondwise changes
are `3.88674184e-7`, `5.62589424e-10`, `7.4e-8` mRy in Dx, Dy, Dz,
and `9e-6` mRy in Jijso. Thus allowed exchange values shift under the
intervention; they are not numerically unchanged.

The reversed-domain reciprocal job was launched, passed the full graph gate
(max bare-h defect `2.7755575616e-17` Ry; full Frobenius defect
`4.7410690851e-15` Ry), and was last observed progressing through its bond
recursions. Its completion and final outputs could not be verified. No full-run
summed odd/even Dx, final all-bond CC/SC decomposition or matching original
endpoint regression is claimed from that incomplete observation.

## Remaining work after execution access is restored

1. Recover the existing `endpoints-reciprocal` and `endpoints-original`
   manifests if available; the shell was scheduled to run both original
   baselines after both reciprocal domains. Do not infer completion from output
   files that may have been copied from the fixtures.
2. Validate both cases' return codes, override announcements, immutable-input
   hashes and identical executable SHA. Recreate fresh runs if the interrupted
   workspace is unavailable; never reuse a partial destination.
3. Run `causal_summary.py` on complete 24-pair endpoint roots. It checks
   manifests, input integrity, pair coverage and matching energy grids.
4. Add the full summed odd/even Dx and CCx/SCx tables, bond/shell Dz/Jij
   comparisons, and complete before/after orbital trace decomposition.
5. Preserve the private numerical archive, update the scientific conclusion,
   and only then mark this stacked PR ready for review. Do not merge.

## Validation and reproducibility

The 23-test suite passes, including shuffled reverse pairs, onsite, Mn-F and F-F,
missing/ambiguous partners, format errors, exclusive writes, marker preservation,
byte idempotence, actual Fortran input loading/rejection, graph checks,
failed/modified-input rejection and known even-scalar/odd-vector leakage parity.
CI commit `d3f5dd328bede51b2d7d2d86f6eca1101df98006` passes all three
platform builds; the diagnostic test step passes on Linux and macOS and is
explicitly skipped on Windows.
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
that signal. The full reciprocal [100] endpoint also shows collapse of summed
Dx to the numerical residual level. Full two-domain confirmation is blocked by
the execution-service disconnection and is not claimed.

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
