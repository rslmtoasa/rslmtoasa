# Stage 1 validation record

Baseline c6e1bb4668a48b17c4f37468cefd949919bedf01, clean checkout.
The legacy FP64 moment loop stores T_m Z for every order. The fast_sp kernel
constructs type-level bond blocks separately. The GPU moments route is guarded.
Physical output integrates each site's contribution before applying the existing
conversion. The historical raw output averages trace seeds.

## Build and reference

Used GNU Fortran 14.2.0, Release (-O3 -fbacktrace -g), GNU OpenMP, Debian
BLAS/LAPACK 3.12.1, MPI off. Both material executables used the same solver
flags, libraries, geometry, energy grid and potential files in each comparison.
The HOH regression builds its own Debug Fortran driver. Stage 1 uses the
historical energy scaling (the default-kind literal in 2 − 0.3 matters at
raw-moment precision). The current conductivity source uses Lorentz damping,
including for the feature. No kernel choice was changed.

## Independent validation

* The dense finite model checks two independent block recurrences and the
  single-recurrence adjoint against spectral polynomials and direct complex
  traces. Both meet the absolute 1e−12 gate. Positive and negative exchange
  fields give Tx=−Bσy, Ty=Bσx and Tz=0. A spin-dependent hopping gives the
  extra full-torque commutator of the hopping itself.
* The compiled test uses native ee, hall, SOC, overlap, H and velocity actions
  on a two-layer slab periodic in x/y only. It checks both total and diagonal
  orbital storage and separate site IDs for identical labels. The largest
  raw block error over 49 order pairs, two representatives and both HOH
  settings was 5.07e−17. The largest complete-cell moment-sum error was
  6.30e−17. The maximum seed/local site difference reached 3.78e−4 in the
  same model units. Moving the representative in-plane preserved the answer.
  A site-specific transverse field in hall also passes with and without HOH.
* Across all 41 reconstructed energies, the largest discrepancy of the
  total divided by the largest seed layer signal was 7.50e−16 without HOH
  and 4.55e−16 with it. Without HOH, the spin-scalar layer had zero local
  torque and a finite seed spectrum. With HOH, the return paths generated
  a finite onsite magnetic term on that layer. The CSV contains the actual
  complex layer spectra in sheet units. Full-grid tests also pass for spin,
  orbital, quadrupole, OAP and octupole accumulations.
* The open slab tests measure both layers and the total. At order 7, the
  total moment residual was 1.28e−2 at a lateral corner, 4.93e−4 at the
  center of width 3, 1.01e−5 at width 5, and 1.17e−7 at width 7. The exact
  reference includes the top and bottom surfaces. Increasing order at fixed
  width brings the artificial edges back into the response.
* The compiled input checks reject invalid projection, random trace, a
  non-charge input, full spin torque and orbital current under local projection.
  B2FeCo local smoke runs executed nsp=1,2,3,4 for spd and nsp=2 for sp.
  spdf output is guarded because the existing formatted writer has an
  18-orbital record. These smoke checks exercise dispatch, not a separate
  reference for every radial spin configuration.

## Cluster size and reciprocity

The 512-site (8³) frozen Pt run uses the repository's legacy structure
constants, 20 moments, 160 LDOS channels and all three periodic flags.
The charge, spin Hall and orbital Hall raw and physical grids reproduce main
exactly across all output files checked by the compact benchmark launcher.
Ratios of peak imaginary to peak real raw conductivity are respectively
1.67e−16, 1.83e−14 and 2.06e−13. For the new onsite Tx the raw seed/local
real discrepancy decreases from 4.15e−25 on 8³ to 1.61e−27 on 12³
(1,728 sites); peak imaginary discrepancies decrease from 2.09e−26 to
3.23e−28. Both spin-z projections were identically zero. Pt is a null
test for the onsite response. Use the magnetic slab for a finite signal.

The reduced-cutoff 4³ Pt and FeCo legacy structure-constant decks failed
the directed Hermiticity check at 3.90e−5 and 1.25e−4 Ry. On the reduced
6³ Pt deck the first rejected pair differed by 1.39e−9 Ry. A nominal 5³
Pt geometry produced only 124 sites and was excluded. The paired small-cell
resource series uses strux_lib for both binaries. Both 8³ Pt structure
constant backends passed the adjoint guard. Periodic flags alone do not
establish reciprocity, a complete representative trace, or size convergence.

The reduced Ni(001) HOH + U deck (U=3 eV, J=0.7 eV, rc=16, ct=2.8, r2=8)
exposed an effective non-Hermitian onsite block at representative 365. The
directed check detects a defect earlier at sites 1 and 99 (9.19e−11 Ry,
legacy; 9.23e−11 with strux_lib). No symmetrization or tolerance relaxation
was used. Quantitative Ni surface torque and its size convergence remain
unverified. The model with a local transverse field tests the assembled
operator algebra; it does not establish independent material-level
provenance for every Hubbard or constraint configuration.

## Cost and architecture

The compiled two-layer HOH fixture at order 250 (18 sites, 18 spin-orbitals,
two traced sites) took 7.02–7.06 s for historical full torque, 4.17–4.32 s
for onsite seed, and 4.10–4.20 s for onsite local. Peak RSS was 60.6–60.9
MiB across these runs. These two samples are indicative; other validation
processes shared the host. The benchmark measures the full test driver,
including reconstruction and writing. A separately timed build measured the
explicit recurrence Hamiltonian calls and input-velocity calls with
system_clock. Those timings exclude output-vertex applications, contractions,
and setup. At order 250 compact storage consumes 36,000,000 bytes for the
two-type moment tensor. The local method reuses the historical left-vector
cache. The optional full orbital storage preserves the adjoint of the entire
orbital block.

Control owns parsing, defaults and input validation. Calculation rejects
reciprocal backends before their moment dispatcher. Recursion owns a single
site-indexed cache for onsite seed, the Hermiticity check, the onsite probe
and the cyclic loop. It reuses select_spin_operator, kubo_ham_vec_matmul,
apply_kubo_operator and transport_moment_diagonal. Conductivity owns the unit
classification and type CSV. The Hamiltonian assembler, SCF, Gamma generation
and energy integration have no feature changes. The old operator-count test
now checks the historical seed loop independently of the new local branch.

## Acceptance status

Draft. The CPU fast_sp spin comparison on 12³ Pt failed its existing 1e−6
threshold; the untouched baseline reproduces the 2.51e−6 normalized
difference. The feature uses FP64 moments. Resolving a pre-existing fast
backend discrepancy or changing the baseline numerical algorithm is outside
Stage 1. The guarded surface case and incomplete independent material-level
Hubbard/constraint checks also limit acceptance. The PR is unmerged.

Reproduction and machine-readable evidence are in
tests/benchmarks/local_projection, tests/regression/test_local_projection*,
tests/regression/test_local_spin_projection.py and docs/validation/local_spin_torque.
