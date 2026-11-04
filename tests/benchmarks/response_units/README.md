# Conductivity normalization and physical response outputs

This change builds on the periodic velocity correction (PR #16). No Fortran
module is added. The H/velocity loops, operator definitions and positive-e
Bastin response-sign convention are unchanged.

## Numerical changes

The common factor is now `4/(pi*a*a)`, where `a=DeltaE/1.7` is the actual
Chebyshev scaling. The previous `16/(pi*DeltaE**2)` assumed `a=DeltaE/2`.
The corresponding prefactor-only change is 0.7225. The three conductivity
scaling expressions now use `2.0_rp-0.3_rp`, avoiding a single-precision
intermediate. Other SCF/DOS scaling paths are untouched.

A local, assumed-shape quadrature replaces conductivity's calls to
`simpson_f`. It integrates quadratic panels to each output energy, including
partial panels, without applying a discretized Fermi step at that endpoint.
It is exact for quadratic polynomials on increasing nonuniform grids. The
first output integral is zero. Both endpoint parities are supported, with a
trapezoid fallback for a two-point grid. The old call could read element 2511
of a 2510-element array at `channels_ldos=2500`. The global SCF integration
routine is unchanged. This conductivity path remains zero-temperature.

Lorentz damping, lambda=6, is retained; misleading Jackson comments are
corrected. No global factor of two for spin is added: all 18 spin-orbitals
are already traced and the spin matrix is sigma/2.

## Inputs and geometry

Optional input group (bulk is the default):

```fortran
&conductivity
  response_geometry = 'bulk'  ! 'bulk', 'sheet', or 'per_atom'
  response_measure = 0.0     ! automatic; positive overrides TOTAL sample measure
/
```

A positive `response_measure` is the total sampled volume in Angstrom^3 for
bulk, or total sampled area in Angstrom^2 for sheet. It is not a per-atom
volume/area. `per_atom` requires zero (no override).

Automatic measures require `pbc=.true.` and a complete generated supercell
`kk=ntot*n1*n2*n3`:

- Bulk: `abs(det(a))*alat**3*n1*n2*n3`.
- Sheet: `norm(cross(a(:,1),a(:,2)))*alat**2*n1*n2`.
- Per atom: the actual number of sites, `kk`, irrespective of periodicity.

For an open/cut/custom sample, explicitly supply the total physical measure
or select `per_atom`. A bulk zero-volume cell is rejected. A planar cell with
zero third lattice vector is valid with sheet normalization.

For a genuinely isolated layer one can use:

```fortran
&lattice
  pbc = .true.
  n3 = 1
  b3 = .false.
  ! Other geometry and in-plane wrapping parameters as appropriate.
/
&conductivity
  response_geometry = 'sheet'
/
```

`n3=1` alone never switches units: one periodically repeated cell may represent
bulk. Conversely, a multilayer slab can use sheet normalization with `n3>1`;
the result includes the entire slab and grows with the number of conducting
layers. Changing vacuum height changes bulk density normalization but not
sheet normalization. Selecting sheet units does not disable interlayer
hopping or make a Hamiltonian two-dimensional. Use in-plane directions for
ordinary sheet conductivity; the code does not silently project directions.

## Trace normalization

`per_type` supplies one trace contribution from each representative site.
The physical output sums these with weights `count(iz==type)/sample_measure`.
This is exact for equivalent sites in each type and remains a representative-
site approximation for inhomogeneous finite samples; weights cannot repair
an inadequate trace estimator.

`random_vec` seeds are normalized by `1/sqrt(kk)`, so they estimate Tr/kk.
Weights are `kk/(sample_measure*random_vec_num)`. Dividing by both the full
sample volume and an extra number of atoms would double-normalize the result.
`per_atom` produces a population-weighted trace average, not an additional
spatial projection onto an arbitrarily chosen output atom.

## Files and units

- Existing `cond_*.out` and atom-symbol `_cond*` files remain bare numerical
  kernels, with energy E-Ef in Ry. Their historical unweighted averaging over
  representative types is retained for comparison. They include the corrected
  prefactor/quadrature, so they are not byte-compatible with old curves.
- `response_total.out`: E-Ef in eV, real and imaginary normalized response.
- `response_orbitals_real.out`, `response_orbitals_imag.out`: E-Ef in eV and
  all 18 orbital/spin trace contributions, whose sum gives the total.
- `response_type_N.out` (per_type): additive population-weighted type
  contributions in the same units as the total, with unambiguous numeric IDs.
- `response_units.txt`: operator names, directions, geometry, trace weights,
  normalization measure, raw dimensions, physical-unit label, conversion
  factor and conventions. New numerical files have comment headers.

Physical tables below assume `linear_in='charge'`. Let q0=hbar^r for an angular
multipole of rank r. Spin and orbital dipoles have r=1, quadrupoles r=2, and
pure/legacy octupole components r=3. This is a polynomial of angular momentum,
not an electric position multipole. The raw matrices remain dimensionless;
q0 is absorbed into the named output unit, never multiplied into tiny floating
point hbar^2 or hbar^3 numbers.

| Output | Bulk unit | Sheet unit | Per-atom unit |
|---|---|---|---|
| Charge, xx or Hall | S/cm | e^2/h | S m^2 per atom |
| Spin/orbital current | (hbar/e) S/cm | e/(2 pi) | (hbar/e) S m^2 per atom |
| Quadrupole current | (hbar^2/e) S/cm | e hbar/(2 pi) | (hbar^2/e) S m^2 per atom |
| Octupole current | (hbar^3/e) S/cm | e hbar^2/(2 pi) | (hbar^3/e) S m^2 per atom |
| Rank-r accumulation | hbar^r / [(V/m) m^3] | hbar^r / [(V/m) m^2] | hbar^r / (V/m) per atom |
| Implemented dipole torque | (hbar/s) / [(V/m) m^3] | (hbar/s) / [(V/m) m^2] | (hbar/s) / (V/m) per atom |

Headers use the equivalent sheet label `(hbar^r) e/h`. The spin unit is
hbar/e, not hbar/(2e); multiply its printed coefficient by two when reporting
in hbar/(2e) S/cm. The physical operator itself must not be divided by two
again.

For two current vertices, raw dimensions are Angstrom^2. An accumulation
output with charge input has Angstrom/Ry; torque with charge has Angstrom.
The code maps an operator to energy/length powers: current (1,1), accumulation
(0,0), torque (1,0). Any pair therefore has raw powers
Ry^(p_out+p_in-2) Angstrom^(l_out+l_in).

An input other than charge remains allowed. Without specifying the conjugate
physical driving field, it is not automatically an electrical response.
`response_*.out` then contains a spatially normalized generalized RAW kernel,
explicitly labeled with these powers and geometry. No fictitious S/cm label,
charge factor, or angular-momentum conversion is applied.

## Dimensional derivation

With E0=1 Ry in joules, l0=1 Angstrom in meters, numerical W=-i[r_c,H_c]
represents physical velocity `(E0*l0/hbar)*W`. A physical observable is
`Q=q0*Q_c`. The common response scale for physical vertices A=a0*A_c and
B=b0*B_c is `hbar*a0*b0/(E0^2*measure)`.

This gives, retaining the existing response sign convention:

- Charge current / charge drive: `(e^2/hbar)*l0^2/measure`.
- Multipole current / charge drive: `(e*q0/hbar)*l0^2/measure`.
- Accumulation / charge drive: `(e*q0/E0)*l0/measure`.
- Torque / charge drive: `(e*q0/hbar)*l0/measure`.

The Ry conversion cancels for current-current response but remains for
accumulation. Bulk current coefficients in the units listed above multiply
the spatially normalized Angstrom^-1 kernel by `(e^2/hbar)*1e8`; sheet current
coefficients multiply it by `2*pi`. Existing math_mod constants are used
consistently for SI/eV conversion. No claim of a newly changed Hall sign is
made: this patch retains the positive-e convention used by the existing kernel.

## Validation and reproduction

Set FC, CMAKE and optionally RS_LMTO_TEST_CMAKE_ARGS as in the periodic-velocity
report. Then run:

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q \
  tests/regression/test_response_units.py \
  tests/regression/test_hoh_transport.py \
  tests/regression/test_octupole_operators.py \
  tests/regression/test_quadrupole_operators.py
```

The compiled production fixture tests all unit families and geometry choices,
unequal type populations, normalized random vectors, raw prefactor, orbital
and type decomposition, gamma/Lorentz/half-weights, nonuniform polynomial
quadrature, explicit measure overrides, generic inputs, planar cells, and
invalid namelists. Debug builds enable bounds checking. The existing transport
suite exercises HOH, impurities, periodic bonds, currents, local operators,
SOC torques, masks, and multipole definitions.

Validation completed: **93 tests passed and 15 subtests passed** with the
bounds-checked production fixtures. The release build also completed four
fccPt smoke runs at `cond_ll=50`: bulk charge xx, bulk spin Hall, bulk orbital
accumulation, and sheet charge xx. All produced 2510 finite energy rows and
passed conversion, orbital-sum, and type-sum checks. The maximum conversion
error relative to each signal scale was below 2.9e-7, consistent with the
limited precision of the legacy raw text files. See [smoke_results.csv](smoke_results.csv).

Reproduce using `python tests/benchmarks/response_units/smoke_fccpt.py`
with the absolute release executable path and output directory as arguments.
Set `OPENBLAS_NUM_THREADS=1`.

These runs reuse frozen fccPt potentials; they do not perform a new SCF cycle.
The sheet case checks geometry and units and is not a converged physical
calculation of a free-standing Pt layer. Timings are smoke-run observations,
not a controlled performance comparison.

