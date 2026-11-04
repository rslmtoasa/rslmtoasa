# Conductivity output units: unit-cell trace

This revision supersedes the population-weighted version of PR #17.
For `cond_calctype='per_type'`, the physical result is

    response = conversion * sum(atomic_subspace_responses) / cell_measure

Each selected atomic subspace enters exactly once. Neither cluster populations
`iz`, cluster size `kk`, nor replication counts `n1,n2,n3` enter these weights.
The selected subspaces must describe the unit cell whose measure is supplied.
The implementation does not infer crystallographic multiplicities or invent
missing atomic contributions.

## Input

```fortran
&conductivity
 response_geometry = 'bulk'
 response_measure = 0.0
/
```

- `bulk`: cell volume in Angstrom^3; default.
- `sheet`: cell in-plane area in Angstrom^2; selected explicitly even for n3=1.
- `response_measure > 0`: explicit **unit-cell** volume/area override.
- Zero: volume `abs(det(a))*alat^3`, area `norm(cross(a1,a2))*alat^2`.
  These are cell measures, also for an open surrounding cluster.
- `per_atom` was removed: averaging over atoms is not the atomic-subspace trace.

## Output

Legacy `cond_*.out` and individual atomic raw files retain their conventions:
energy relative to EF in Ry; the legacy total is an arithmetic average over
selected subspaces. Physical `response_total.out` instead reports their SUM,
divided by the cell measure. Do not convert the legacy averaged total without
restoring its number-of-subspaces factor. Individual `response_type_N.out`
files are additive atomic contributions with the SAME cell denominator.
`response_orbitals_real.out` and `response_orbitals_imag.out` resolve the 18
spin-orbitals. New energy columns are eV. `response_units.txt` records the cell
measure, conversion, trace weights, operators and sign convention.

For charge input, with angular multipoles Qphys=hbar^r Qcode:

| Output | Bulk | Sheet |
|---|---|---|
| Charge xx/Hall | S/cm | e^2/h |
| Spin/orbital current | (hbar/e) S/cm | e/(2*pi) |
| Quadrupole current | (hbar^2/e) S/cm | e*hbar/(2*pi) |
| Octupole current | (hbar^3/e) S/cm | e*hbar^2/(2*pi) |
| Rank-r accumulation | hbar^r / [(V/m) m^3] | hbar^r / [(V/m) m^2] |
| Dipole torque | (hbar/s) / [(V/m) m^3] | (hbar/s) / [(V/m) m^2] |

The spin matrices already include sigma/2. These are angular-momentum powers,
not electric position multipoles. The existing positive-e sign is retained.
For a non-charge input operator the output is labeled generalized raw response;
no electric-field interpretation is imposed.

For random_vec, the original seed average is retained and divided by the cell
measure. No factor kk is introduced. Equivalence of that normalized stochastic
estimator to the full cell trace is NOT established by this units patch;
its absolute conductivity normalization needs a separate estimator audit.

## Numerical scope and limits

The historical `16/(pi*DeltaE^2)` prefactor and `(2-0.3)` rescaling are restored.
The proposed `4/(pi*a*a)` correction and full-precision rescaling are deferred.
No Hamiltonian, velocity or recursion implementation changes remain in this PR.

The former integration replacement is removed. `conductivity_integral` retains
simpson_f's uniform-grid Simpson panels and zero-temperature Fermi factor
(including one-half at EF). It bounds panels by the actual array extent; the
old `nv1+10` extent could read beyond the allocated energy grid. It discards an
unpaired final interval exactly as a bounded complete-panel Simpson loop does.
Thus raw values agree with the legacy algorithm on valid complete panels, but
undefined out-of-bounds contributions cannot be preserved. This small safety
repair is explicitly distinct from units conversion. No global SCF integrator
is modified. A separate numerical audit should address the prefactor and
integration accuracy; physical labels do not certify those legacy choices.

## Validation

The compiled fixture checks all 16 operator names/aliases for bulk and sheet,
cell-area independence of the third vector, explicit measures, malformed
inputs, and unchanged weights/output when cluster populations, size, boundary
flag and replication counts change. Independent Python Simpson/Fermi weights
check the retained numerical convention, with debug bounds checking.

Run:

```sh
python -m pytest -q tests/regression/test_response_units.py tests/regression/test_hoh_transport.py tests/regression/test_octupole_operators.py tests/regression/test_quadrupole_operators.py
```

The previous fccPt smoke CSV describes the superseded implementation and has
been removed. The smoke script remains reproducible for the revised executable.

Revised validation: **78 tests passed and 15 subtests passed** (GNU Fortran
13.3, Debug bounds checks, OpenMP, reference BLAS/LAPACK). This includes direct
agreement with the unmodified production simpson_f on valid complete panels.
The four fccPt smoke runs reported for the previous version were not rerun for
this revision and are not claimed as evidence for it.
