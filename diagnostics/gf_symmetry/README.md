# Magnetic Green-function diagnostics

Read `../../MN2F2_GF_SYMMETRY_AUDIT.md` for the derivation and measured findings.
These tools read data and report violations; they never project production matrices
or change exchange formulas. Python requires NumPy (already a repository dependency).

## Tests and build

```bash
python -m unittest discover -s diagnostics/gf_symmetry -p 'test_*.py' -v
cmake -S . -B build_gf_audit -DCMAKE_Fortran_COMPILER=gfortran -DENABLE_MPI=ON
cmake --build build_gf_audit -j 4
```

Use a separate baseline worktree at the SHA recorded in the audit for an unmodified
build. Some installations require explicit BLAS/LAPACK paths; when linking with a
combined OpenBLAS library, pass its path as both `BLAS_LIBRARIES` and
`LAPACK_LIBRARIES`. No repository build-system workaround beyond registering the
new diagnostics module is included here.

## Keep unpublished data outside Git

Unpack the supplied archive outside the repository, for example in
`/absolute/private/mn2f2`. `scratch/` and `gf-audit-*.dat` are ignored, but an external
working directory is preferable. Do not add the original fixture, potentials, raw
output, dump files, detailed JSON reports, or execution manifests to a public PR.

The runner creates fresh copies and refuses to reuse a destination. It verifies
`post_processing='exchange'` and refuses an active `pre_processing` assignment.
It preserves the 2000-point grid, recursion length, saved potentials and Hamiltonian
flags. `--one-pair` and `--pair-count` are explicitly reduced *pair-list* smoke tests;
full runs omit both options.

```bash
# Serial first-pair gate. Repeat without --audit/--h-action for regression.
python diagnostics/gf_symmetry/run_saved.py \
  --fixtures /absolute/private/mn2f2 \
  --binary build_gf_audit/bin/rslmto.x \
  --destination /absolute/private/smoke \
  --cases 100 1-00 --one-pair 1 --audit 1 --h-action --threads 1

# All four saved jobs, all 24 original pairs, four OpenMP threads per MPI rank.
# Choose rank count to fit memory; this code has substantial per-rank allocations.
python diagnostics/gf_symmetry/run_saved.py \
  --fixtures /absolute/private/mn2f2 \
  --binary build_gf_audit/bin/rslmto.x \
  --destination /absolute/private/audited \
  --audit all --h-action --mpi 2 --threads 4

# Same full jobs, instrumentation disabled, in another fresh destination.
python diagnostics/gf_symmetry/run_saved.py \
  --fixtures /absolute/private/mn2f2 \
  --binary build_gf_audit/bin/rslmto.x \
  --destination /absolute/private/control --mpi 2 --threads 4
```

For direct execution from a copied `exchange_jobs/Mn1` directory:

```bash
RSLMTO_GF_AUDIT=all RSLMTO_GF_AUDIT_H=1 \
  OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  mpiexec -n 2 /absolute/build_gf_audit/bin/rslmto.x input.nml
```

On an existing SLURM allocation, the runner can use `--mpiexec srun`, with
`--mpi "$SLURM_NTASKS" --threads "$SLURM_CPUS_PER_TASK"`. Keep OpenBLAS at one
thread. Recursion buffers are still partly replicated, so increasing MPI ranks
can increase total memory even though GF arrays are distributed.

Controls:

| Variable | Meaning |
| --- | --- |
| `RSLMTO_GF_AUDIT` | Unset disables. Positive global pair index selects one pair; `all` selects the first 24. |
| `RSLMTO_GF_AUDIT_H=1` | Also capture actual first Hamiltonian action on the real +/- seeds. Used in Chebyshev exchange recursion. |
| `RSLMTO_GF_AUDIT_STRIDE` | Increase energy sampling stride. Sampling is capped at 64 full snapshots per pair. All energy integrands are always logged. |

The sole pair owner writes `gf-audit-gf-r<R>-p<P>.dat` and, optionally,
`gf-audit-h-r<R>-p<P>.dat`. All routines run outside OpenMP regions. Different
processes do not share a diagnostic file. Production output still follows the
repository's existing I/O behavior. `P` metadata records actual cluster indices,
types and displacement; file order is never a symmetry map.

`I` records contain energy followed by **CCx, CCy, CCz, SCx, SCy, SCz**, before
integration and its `2*1000/(4*pi)` conversion. `M` records give index, energy,
label and dimension followed by rows of `(real,imaginary)` entries. Raw 18x18 GFs,
raw 9x9 Pauli components, all 16 split matrices, and two Delta matrices are sampled.
`Hp/Hm` are the actual scaled first actions; `Op/Om`, `Ep/Em`, `Sp/Sm` separately
capture HoH, energy-centre and SOC contributions when HoH is enabled. `A` records
supply the full cluster geometry for mapping Hamiltonian output rows.

## Analysis

```bash
python diagnostics/gf_symmetry/fixtures.py /absolute/private/mn2f2 \
  --output /absolute/private/inventory.json
python diagnostics/gf_symmetry/analyze.py /absolute/private/audited/100 \
  --minus /absolute/private/audited/1-00 \
  --output /absolute/private/domain100.json
python diagnostics/gf_symmetry/analyze.py /absolute/private/audited/110 \
  --minus /absolute/private/audited/1-1-0 --neel 110 \
  --output /absolute/private/domain110.json
python diagnostics/gf_symmetry/sbar_audit.py \
  /absolute/private/mn2f2/100/exchange_jobs/Mn1 \
  --h-dump /absolute/private/audited/100/gf-audit-h-r0-p1.dat \
  --output /absolute/private/sbar100.json
python diagnostics/gf_symmetry/fixtures.py /absolute/private/control/100 \
  --regression /absolute/private/audited/100 \
  --output /absolute/private/regression100.json
```

The analyzer reports absolute max/mean/Frobenius and relative residuals, leaking
orbital indices and mirror parities, raw/Pauli/split consistency, Delta commutators,
CC/SC traces by allowed/forbidden factor and s/p/d parity transition, inter-domain
matrix parity, and mapped C2x/C2zT/inversion covariance. Use `--neel 110` for the diagonal domain; C2x is then explicitly excluded as a same-domain symmetry. Missing mapped bonds are
reported explicitly. Domain comparisons reject unmatched grids and Fermi levels.
Hamiltonian checks cover full output rows of the two captured columns, including
component and pair Hermiticity; passing them does not prove symmetry of every
Hamiltonian column. Default failure thresholds are absolute `1e-12` plus relative
`1e-9` of the compared matrix scale. Inspect raw residuals rather than treating a
threshold label as a physical conclusion.

`fixtures.py --regression` compares matched directed-bond values at the existing
printed precision. It does not mistake a compiler/parallel-reduction difference
between supplied outputs and a new build for an instrumentation regression.

`sbar_audit.py` independently reads the **binary** `sbar.` matrices. Rounded
`view.sbar.` vector labels only identify the geometric map. It checks reciprocity,
reconstructs bare spinor hopping from saved widths/moments, and optionally compares
it with the first-action capture. The accepted binary format is checked, not guessed.

Run only on complete output files. Detailed reports remain private. No script
changes potentials, reruns SCF, symmetrizes structure constants, or discards Dx.
