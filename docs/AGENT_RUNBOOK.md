# Running RS-LMTO-ASA — runbook for AI agents

How to set up, run and check a calculation. This file is about *using* the
code. For changing the code, read [`CLAUDE.md`](../CLAUDE.md) and
[`DEVELOPER_MAP.md`](DEVELOPER_MAP.md).

The keyword reference is [`docs/source/keywords/`](source/keywords/). Every
namelist key in `source/` has an entry there (a quick-tier test checks this).
Use it for meaning, type, default and units. This runbook doesn't repeat it.

## 1. Build

```bash
cmake -S . -B build -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The executable is `build/bin/rslmto.x`. Options and MPI/Intel builds are in
the [README](../README.md#build-from-source). The executable has no `--help`.
Its only argument is the input file name (default `input.nml`).

## 2. Pick the workflow and a starting deck

Always start from the closest deck in `example/`. Don't write a deck from
scratch. The workflow is selected by `&calculation` and `&control`:

| Goal | `&calculation` | `calctype` | Starting deck | Doc page |
|---|---|---|---|---|
| Bulk SCF, recursion | `pre_processing='bravais'` | `'B'` | `example/bulk/bccFe`, `example/bulk/Si` | `user_guide/examples/bulk_bccfe.rst`, `bulk_si.rst` |
| Bulk SCF, noncollinear | `pre_processing='bravais'` | `'B'` | `example/bulk/Mn3Sn` | `bulk_mn3sn.rst` |
| Bulk SCF with k-space parts | `pre_processing='bravais'` | `'B'` | `example/k_space_scf/bccFe` | `kspace_scf_bccfe.rst` |
| Surface SCF | `pre_processing='buildsurf'` | `'S'` | `example/surface/fccCu001` | `surface_fcccu001.rst` |
| Interface / layered SCF | `pre_processing='buildinterface'` | `'L'` | `example/interface/fccCu111_AA` | `interface_fcccu111.rst` |
| Impurity in bulk | `pre_processing='newclubulk'` | `'I'` | `example/impurity/B2FeCo` | `impurity_b2feco.rst` |
| Exchange J_ij, D_ij | `post_processing='exchange'` | `'B'` | `example/exchange/bccFe` | `exchange_bccfe.rst` |
| Conductivity | SCF, then `post_processing='conductivity'` | `'B'` | `example/conductivity/fccCu` (+ `sigma_*/`) | `conductivity_fcccu.rst` |
| Band structure, DOS (k-space) | `post_processing='band_structure'` / `'density_of_states'` | `'B'` | `example/band_structure/bccFe`, `example/density_of_states/bccFe` | `keywords/reciprocal.rst` |
| Spin spirals / frozen magnons | `post_processing='frozen_magnon'` | `'B'` | `example/k_space_spiral/` (see its README) | `keywords/frozen_magnon.rst` |

Doc pages are under `docs/source/`. The full list of `pre_processing` and
`post_processing` values, including the removed ones, is
[`DEVELOPER_MAP.md` §1](DEVELOPER_MAP.md#1-execution-skeleton).

Not every directory in `example/` is a curated starting point. Use the decks
in the table. Of these, the regression tests read only `example/bulk/bccFe`,
`example/surface/fccCu001` and `example/interface/fccCu111_Avac`.

### Starting values

When a deck needs changing, start from these values (developer guidance):

| Key | Namelist | Baseline | Production |
|---|---|---|---|
| `beta` | `&mix` | 0.05–0.2 for stable systems | same |
| `nstep` | `&self` | 25–50 for bulk metals, magnetic or not (metals converge slower); about 100 for surfaces | same |
| `lld`, block recursion | `&control` | 21 | about 31 |
| `lld`, Chebyshev | `&control` | 150 | about 200 |

For production bulk runs, use a periodic cluster: `pbc`, `b1`, `b2` and `b3`
all `.true.`, with `n1`, `n2`, `n3` set (`&lattice`, see
[`keywords/lattice_geometry.rst`](source/keywords/lattice_geometry.rst)).

## 3. Files in a calculation directory

| File | Role |
|---|---|
| `input.nml` | `&calculation`, `&lattice`, `&atoms`, `&control`, `&self`, `&energy`, `&mix`, ... |
| `lattice.nml` | Geometry, read when `crystal_sym = 'file'` |
| `<label>.nml` | `&element` and `&par` for each label in `&atoms` (looked up in `database`) |
| `<label>_out.nml` | Atom state from a previous run; **read instead of `<label>.nml` when present**, unless `fresh_start = .true.` (`&control`) |

`alat` is in Ångström. `nsp = 1` means collinear scalar-relativistic, not
nonmagnetic (2 = collinear + SOC, 3 = noncollinear, 4 = noncollinear + SOC).

## 4. Run protocol

1. Copy the whole example directory out of the repository. Never run inside
   `example/`, because the run writes files there.
   ```bash
   mkdir -p runs && cp -r example/bulk/bccFe runs/bccFe && cd runs/bccFe
   ```
2. Change one thing at a time from the starting deck. Record what you changed.
3. Check for stale `*_out.nml` files: they override `<label>.nml` (§3).
4. Run and capture all output:
   ```bash
   OMP_NUM_THREADS=4 ../../build/bin/rslmto.x input.nml > run.log 2>&1
   echo "exit $?"
   ```
5. Check the result (§5) before using any number.

MPI: `mpirun -np N` with the MPI runtime of the build (README, "Parallel
execution").

## 5. Was the run successful?

**The exit code is 0 even when SCF did not converge.** A non-zero exit code
means a fatal error (`<fatal>` line in `run.log`). Check convergence in the log:

```bash
grep -a 'Converged!' run.log          # present: SCF converged (final delta printed)
grep -a 'Not converged! Diff=' run.log | tail -3   # one line per unconverged iteration
grep -a '<fatal>' run.log             # reason for a non-zero exit
```

`Converged!` is case-sensitive and doesn't match `Not converged!`. Log lines
may contain ANSI colour codes (build option `COLOR=ON`), hence `grep -a`.

If `Converged!` is missing, the run used up `nstep` iterations. Report it as
not converged and give the last `Diff=` value. Don't present its moments or
energies as results.

After SCF, `report.out` holds the summary:

```bash
grep -E 'Total energy of system|Fermi energy|Spin moment of atom|Charge transfer at atom|Total RMS Diff' report.out
```

Output files per workflow: [`user_guide/output_files.rst`](source/user_guide/output_files.rst).
Their columns and units depend on the writing routine. Read that page before
converting values.

`Converged!` refers to the SCF criterion only (`conv_thr`, `&self`).
Convergence with respect to cluster size, recursion depth, energy mesh and
k-mesh is a separate study; do it only on request.

## 6. Post-processing runs

Post-processing (`exchange`, `conductivity`, k-space routes) doesn't repeat
the SCF. It reads the atom state from `database`:

1. Converge the SCF first (§5).
2. Copy the converged `<label>_out.nml` into the post-processing directory.
3. Set `fermi` (`&energy`) to the converged Fermi energy from `report.out`.
4. Keep geometry and labels identical to the SCF run.

Exchange pair indices (`njij`, `ijpair` in `&lattice`) refer to the generated
cluster. Recheck them after any geometry change.

## 7. Common fatal messages

| `run.log` message | Cause |
|---|---|
| `file <name> not found` | Input, `lattice.nml` or `<label>.nml` missing from the run or `database` directory |
| `lattice%calctype must be one of: 'B', 'S', 'I', 'L'` | Invalid `calctype` |
| `post_processing='...' was removed` | Value from the archived campaign; see `DEVELOPER_MAP.md` §1 |
| `Invalid &<group> namelist; iostatus = ...` | Misspelled key, wrong type or bad syntax in that group; compare against `docs/source/keywords/` |

A namelist group that is missing from the file isn't an error: all of its
keys keep their defaults silently.

## 8. Reporting a run

Give: commit, build options, deck directory, the changes from the starting
deck, exit code, converged yes/no with the final `Diff=`, and the requested
numbers with their `report.out` or output-file source. Don't call results
"validated" or "converged" beyond what §5 shows.
