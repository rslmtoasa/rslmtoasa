# RS-LMTO-ASA: bcc Fe from SCF to exchange interactions

**First draft · 2 hours · limited-compute edition · based on `main`**

This hands-on tutorial follows one material through self-consistent electronic
structure, a density-of-states plot, and two magnetic exchange interactions.
The aim is to understand inputs, checkpoints, convergence, and outputs before
attempting larger calculations.

Code basis: [`main` at `2fccbffc`](https://github.com/rslmtoasa/rslmtoasa/tree/2fccbffcde6f0f8eb1bd997c47a839af0f9fc887).
Use the [`tutorial/codespaces` environment](https://codespaces.new/rslmtoasa/rslmtoasa/tree/tutorial/codespaces),
which adds installation and plotting tools to that code. This draft does not
use `main_rc` features. Record `git rev-parse HEAD` with saved results.

## Learning goals

By the end, participants should be able to:

- Edit Fortran namelist input files and run the executable from a terminal.
- Distinguish an SCF stopping criterion from an iteration limit and from numerical-parameter convergence.
- Restart using a saved atom state and identify the role of the Fermi energy.
- Plot the DOS and interpret the geometry and units of two selected `Jij` values.
- Export inputs, logs, and outputs so the calculation can be inspected later.

Prerequisites: basic magnetism/DFT concepts and elementary terminal use. Python
experience helps with plotting, but there is no Python interface to the solver:
Python reads its output files.

## Schedule

| Minutes | Activity | Result |
| --- | --- | --- |
| 0–10 | Physical picture: ASA, real-space recursion, SCF, magnetic force theorem | Know what is being approximated |
| 10–20 | Open Codespace; inspect namelists; verify the executable | A working environment |
| 20–50 | SCF from the supplied checkpoint and from the original atom input | Residual history and a saved state |
| 50–65 | Plot convergence and the DOS; discuss restart files | Two figures and a convergence decision |
| 65–75 | Introduce exchange pairs, geometry, and the reference electronic state | Understand the exchange input |
| 75–100 | Transfer the SCF state; compute and plot two `Jij` values | A short exchange table and figure |
| 100–110 | Compare results; discuss numerical sensitivity and limitations | A plan for a convergence study |
| 110–120 | Download results; discuss extensions; questions | A reproducible result bundle |

Calculation waiting time is used for discussion. Do not add a large parameter
scan to fill a faster group's remaining time.

## Instructor preparation and compute budget

1. Have students create their own GitHub accounts and launch the Codespace
   **before the session**. Initial image download and compilation are outside
   the two-hour teaching plan. Confirm they use the personal billing owner and
   have available free quota.
2. Use the 2-core environment, `OMP_NUM_THREADS=2`, and
   `OPENBLAS_NUM_THREADS=1`. Run one calculation at a time per student. No GPU,
   MPI, or desktop display is needed.
3. Pilot the exact handout on a fresh Codespace. Set a **three-minute wall-time
   cap per run**; move to prepared results if a run takes longer. Do not reduce
   cluster size or change atom-pair indices ad hoc during class.
4. Prepare a downloadable fallback containing the actual inputs, `run.log`,
   `report.out`, `Fe_scf.nml`, DOS, and exchange outputs. Include the code commit
   and settings. The fallback archive is an instructor preparation item, not
   an assumed file already present in this repository.
5. Keep the supplied input and atom files intact. Work in new `runs/tutorial-*`
   directories. Commands below assume these directories do not already exist;
   choose new names for a second attempt.

Budget: a seeded check, at most 20 fresh SCF iterations, one short restart check,
then two exchange pairs. The optional mixing comparison adds at most 20 SCF
iterations. Allow **5–10 minutes total calculation time per student** as a
planning allowance, to be confirmed by the Codespace pilot. A Codespace left
running for the whole two hours consumes about **4 core-hours**, plus setup
usage, even when students are only editing or discussing.

## 1. Find your way around (10–20 min)

Open **Terminal → New Terminal** in Codespaces. From the repository root:

```bash
cd "$(git rev-parse --show-toplevel)"
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
command -v rslmto.x
python3 -c "import numpy, matplotlib; print('Plotting tools ready')"
git rev-parse HEAD
```

If the executable is missing, wait for setup to finish or inspect the creation
log. The environment guide is `.devcontainer/README.md`. Do not switch the
Codespace to `main`: that branch does not contain the tutorial setup files.

Open `example/bulk/bccFe/input.nml` in the editor, without changing it.

| Namelist | What to inspect |
| --- | --- |
| `&calculation` | `pre_processing = 'bravais'` selects the SCF workflow |
| `&lattice` | bcc structure, `alat = 2.86120` Å, and real-space environment |
| `&atoms` | `Fe` label and the location of the atom inputs |
| `&self` | `nstep` is the iteration cap; `conv_thr` is the SCF threshold |
| `&control` | `recur = 'block'`, `lld = 21`, and `nsp = 2` |
| `&energy` | Energy interval, Fermi energy, and sampling |
| `&mix` | Broyden mixing and its `beta` parameter |

Retain `nsp = 2` for this first lesson: collinear magnetism with spin–orbit
coupling. **`nsp = 1` means collinear scalar-relativistic, not nonmagnetic.**
Keep the supplied geometry, recursion depth, and energy sampling throughout
the main exercise.

## 2. See what the supplied checkpoint does (20–25 min)

The bulk example includes both an atom input (`Fe.nml`) and a saved state
(`Fe_out.nml`). The saved state is loaded automatically when present.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs
cp -r example/bulk/bccFe runs/tutorial-seeded
cd runs/tutorial-seeded
timeout 180 rslmto.x input.nml > run.log 2>&1
```

Open `run.log` and `report.out`. The supplied state can satisfy convergence
in a single iteration; this is a warm-start check, not an SCF-from-scratch
exercise. Ask: what would this tell us about a new structure or a changed
numerical parameter? What would it leave untested?

If a command returns code 124, the three-minute limit was reached. Use the
instructor's prepared results instead of repeatedly restarting a timed-out run.

## 3. Follow an SCF trajectory (25–50 min)

Create a separate copy and move its supplied checkpoint aside. This falls back
to the existing `Fe.nml` starting parameters; it is **not** a new atomic-density
or `cold = .true.` initialization.

```bash
cd "$(git rev-parse --show-toplevel)"
cp -r example/bulk/bccFe runs/tutorial-scf
cd runs/tutorial-scf
mv Fe_out.nml Fe_out.supplied.nml
```

In this directory's `input.nml`, replace the existing `&self` and `&mix`
sections with the following. Do not append duplicate namelist sections.

```fortran
&self
nstep = 20
conv_thr = 5.0d-8
/
&mix
beta = 0.1
mixtype = 'broyden'
/
```

All other settings stay as supplied. Save the file, then run:

```bash
cp input.nml input.initial.nml
timeout 180 rslmto.x input.nml > run.log 2>&1
```

This teaching threshold is ten times looser than the code's default
`5.0d-9`. It limits classroom work; it is not a production-accuracy recommendation.

Read the log. `Not converged! Diff= ...` reports a residual of mixed band-moment
variables, **not a total-energy difference in Ry**. This build computes it as
`sqrt(sum((qia_old - qia_new)**2)) / (6*nrec)` over the first 12 components.
`Converged!` means this residual is below `conv_thr`. A normal process exit or
20 completed iterations alone does not establish convergence.

Record:

| Quantity | Your result |
| --- | --- |
| Number of SCF iterations | |
| Last residual and whether `Converged!` appears | |
| Total energy from `report.out` (Ry) | |
| Spin moment of atom 1 (mu_B) | |
| Orbital moment of atom 1 (mu_B) | |
| Final Fermi energy in `input.nml` (Ry) | |

The code updates the Fermi energy in `input.nml` during SCF, which is why the
initial input was saved separately. `report.out` summarizes the final state;
it is not a table of energies for every iteration.

**Discussion:** Does the residual decrease monotonically? Why might Broyden
mixing be faster than using a very small mixing parameter? What distinguishes
SCF convergence at fixed settings from convergence with respect to `lld`, the
energy mesh, or the real-space environment?

**Optional comparison:** half the class can run a second fresh copy with
`beta = 0.01`, the same checkpoint removed, and the same 20-step cap. Compare
iteration histories; do not increase the iteration limit to force every run
to converge. Use the converged `beta = 0.1` state for the following exercise.

### Restart check

If the main run converged, preserve its last SCF checkpoint and original log,
then rerun from that state. If it did not converge, perform at most one further
20-step block; otherwise move to the instructor's converged fallback.

```bash
cp Fe_scf.nml Fe_out.nml
timeout 180 rslmto.x input.nml > restart.log 2>&1
```

The restart loader reads `Fe_out.nml`; it does not automatically select an
arbitrary `_scf.nml` file. The same geometry, atom labels, and compatible
settings must accompany the checkpoint. Compare `restart.log` with `run.log`.

## 4. Plot the convergence and DOS (50–65 min)

Run this in `runs/tutorial-scf`. It reads the original run's residual history
and the current DOS, strips terminal colour codes, and saves two PNGs.

```bash
python3 - <<'PYTHON'
from pathlib import Path
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

text = re.sub(r'\x1b\[[0-9;]*m', '', Path('run.log').read_text())
pattern = r'(?:Not converged!\s*Diff=\s*|Converged!\s*)([0-9.EeDd+\-]+)'
residuals = [float(x.replace('D', 'E').replace('d', 'e'))
             for x in re.findall(pattern, text)]
if not residuals:
    raise RuntimeError('No SCF residuals found; inspect run.log')
fig, ax = plt.subplots()
ax.semilogy(range(1, len(residuals) + 1), np.maximum(residuals, 1e-12), 'o-')
ax.axhline(5e-8, color='black', linestyle='--', label='Teaching threshold')
ax.set(xlabel='SCF iteration', ylabel='Printed mixing residual')
ax.legend()
fig.tight_layout()
fig.savefig('scf-convergence.png', dpi=150)

energy, dos = np.loadtxt('totaldos.out', unpack=True)
fig, ax = plt.subplots()
ax.plot(energy, dos)
ax.axvline(0, color='black', linestyle='--')
ax.set(xlabel='E - E_F (Ry)', ylabel='Total DOS (states/Ry)')
fig.tight_layout()
fig.savefig('dos.png', dpi=150)
PYTHON
```

Click the images in the file explorer. In this main-branch writer,
`totaldos.out` has two columns: energy already shifted by the Fermi energy at
DOS writing, and the spin-summed total DOS. Do not subtract `E_F` again.
The Fermi estimate may change during the iteration; use a converged state for
interpretation. The plot is not spin-resolved.

Ask: what does nonzero DOS near zero energy suggest? Does a spin-summed DOS
alone reveal the spin moment? Why should a converged moment not guarantee a
converged exchange interaction?

## 5. Prepare the exchange calculation (65–75 min)

Explain `Jij` as the response to small rotations around a reference electronic
state, rather than a new SCF calculation for two rotated spins. The two
selected pairs provide a small, useful result without computing an entire
exchange-shell table.

Use the complete input from `tests/postproc/cases/exchange/bccFe`. The shorter
public `example/exchange/bccFe/input.nml` omits the explicit `ntype`, `ct(1)`,
and `r2` entries present in this tested deck; use the complete deck for class.

```bash
cd "$(git rev-parse --show-toplevel)"
cp -r tests/postproc/cases/exchange/bccFe runs/tutorial-jij
cp runs/tutorial-scf/Fe_scf.nml runs/tutorial-jij/Fe_out.nml
cd runs/tutorial-jij
```

Only transfer a checkpoint from a run that actually satisfied the selected
SCF threshold. If using an instructor fallback, copy its matched checkpoint
and input settings instead and label the result as a fallback calculation.

Open the exchange `input.nml`. Its key settings are:

```fortran
&calculation
post_processing = 'exchange'
verbose = T
/
```

Inside its existing `&lattice` section, retain these entries and all other
supplied geometry settings:

```fortran
ntype = 1
ct(1) = 3.0d0
r2 = 9.00d0
njij = 2
ijpair(1, :) = 1, 2634
ijpair(2, :) = 1, 2635
```

The pair numbers index this particular generated cluster. They are not
universal neighbour labels. Keep the supplied `ndim = 50000`, `rc = 080`, and
bcc geometry. Confirm the pair displacements in the log after running.

Make the electronic settings consistent with your SCF state. The following
helper copies the Fermi energy and energy-mesh settings from its final input.
Both decks already use `nsp = 2`, block recursion, and `lld = 21`; retain those.

```bash
python3 - <<'PYTHON'
from pathlib import Path
import re
scf = Path('../tutorial-scf/input.nml').read_text()
p = Path('input.nml')
text = p.read_text()
for key in ('fermi', 'energy_min', 'energy_max', 'channels_ldos'):
    value = re.search(r'^\s*' + key + r'\s*=\s*([^!\n]+)', scf, re.M).group(1).strip()
    text = re.sub(r'(?m)^(\s*' + key + r'\s*=\s*)[^!\n]+',
                  lambda match: match.group(1) + value, text)
text = text.replace('nstep = 100', 'nstep = 1')
p.write_text(text)
PYTHON
```

Here the `exchange` post-processing path uses the saved potentials directly;
`nstep` does not request another SCF cycle. Merely setting `nstep = 1` would
not make an SCF input into an exchange input.

## 6. Compute and interpret two exchange interactions (75–100 min)

From `runs/tutorial-jij`:

```bash
cp input.nml input.used.nml
timeout 180 rslmto.x input.nml > run.log 2>&1
cat jij.out
```

Check the run completed and both output rows exist. The routine can also write
DMI, anisotropic, and perturbation-order outputs. Use **`jij.out`** for this
exercise; `jijso.out` and `jijfo.out` label different perturbation orders, not
SOC-on and SOC-off runs.

For this `main` calculation path, the columns of `jij.out` are:

| Column | Meaning |
| --- | --- |
| 1–2 | Symbolic atom-type indices (`iz`), not the cluster pair numbers |
| 3–5 | Displacement from i to j in units of the lattice parameter |
| 6 | Exchange value in mRy |
| 7 | Pair distance in units of the lattice parameter |

The writing routine is `calculate_exchange` in `source/exchange.f90`.
For the selected geometry, expect distances close to `sqrt(3)/2` and `1` in
units of `a`: nearest and next-nearest neighbour distances in bcc Fe. Multiply
by `2.86120` to obtain Å. These are two representatives, not all eight first
and six second neighbours of a bcc site.

```bash
python3 - <<'PYTHON'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rows = np.atleast_2d(np.loadtxt('jij.out'))
r_angstrom = rows[:, 6] * 2.86120
j_mry = rows[:, 5]
j_mev = j_mry * 13.605693
print('Pair     r (Angstrom)       J (mRy)       J (meV)')
for pair, (r, j, mev) in enumerate(zip(r_angstrom, j_mry, j_mev), start=1):
    print(f'{pair:4d} {r:16.6f} {j:13.6f} {mev:13.6f}')
fig, ax = plt.subplots()
ax.plot(r_angstrom, j_mev, 'o', linestyle='none')
ax.axhline(0, color='black', linewidth=0.8)
ax.set(xlabel='Pair distance (Angstrom)', ylabel='J from jij.out (meV)')
fig.tight_layout()
fig.savefig('jij.png', dpi=150)
PYTHON
```

Treat the sign using an explicitly stated model convention. For discussion,
`H = -sum_(i<j) J_ij e_i dot e_j` has positive `J` favouring parallel unit
vectors. Before exporting the code's values into another spin model, verify
its pair-counting, sign, and spin-length normalization; this handout does not
establish that mapping or an extra factor of two. Do not estimate a Curie
temperature from these two pairs alone.

**Questions:** Which representative pair has larger exchange? What does the
bcc coordination imply for a shell sum? Which numerical parameters might
change `J` even if the SCF residual is small? Why would changing the Fermi
energy independently of the checkpoint be questionable?

## 7. Compare, save, and wrap up (100–120 min)

Compare the class's iteration counts, moments, DOS, and two exchange values.
If two results differ, first compare the code commit, starting state,
convergence status, Fermi energy, `lld`, energy interval, and sampling. Do not
assume an installation error or silently copy another student's result.

From the repository root, bundle the working directories:

```bash
cd "$(git rev-parse --show-toplevel)"
tar -czf tutorial-results.tar.gz -C runs \
    tutorial-seeded tutorial-scf tutorial-jij
```

Right-click the archive in the file explorer and choose **Download**. Include
the recorded commit identifier with the archive. Stop the Codespace at
https://github.com/codespaces when finished. Download important files before
deleting it; stopped Codespaces still use storage allowance.

**Exit questions:** Did your SCF calculation converge? Which file would you
use to restart it? Which state did your exchange calculation use? What further
checks would you need before using the exchange values in a publication?

## Suggested further topics, ordered by compute cost

| Topic | Suggested format | Compute impact |
| --- | --- | --- |
| DOS, magnetization, and orbital character | Discuss the computed DOS; optionally inspect `Fe_orbital_dos.out` | No new solver run |
| Restart reproducibility | Compare the short warm restart with the final fresh run | One short run; already included |
| Inversion symmetry and DMI in bulk bcc Fe | Inspect `dij.out`; discuss why symmetry suggests small DMI and why numerical noise matters | No new solver run |
| Mixing strategy | Compare fresh `beta = 0.01` and `0.1` runs with the same iteration cap | At most 20 additional SCF steps |
| Recursion-depth sensitivity of `Jij` | Give different groups one `lld` value; compare against instructor-prepared data | One short exchange run per group; an incomplete convergence study |
| Energy-mesh sensitivity of exchange | One assigned mesh per group, retaining the same converged reference state and energy interval | One short exchange run per group |
| SOC and orbital moments | Compare independently converged `nsp = 1` and `2` reference data | Precomputed in the first session; live later |
| Lattice parameter and magnetovolume effects | Discuss a small prepared volume series | Precomputed; independent SCF runs needed |
| Surfaces and impurities | Preview the existing Cu surface and B2 FeCo impurity examples | Follow-up session; additional geometry and convergence work |
| Transport, DMI tensors, and noncollinearity | Explain the questions addressed, using prepared figures | Follow-up session; more theory and convergence checks |

For a first two-hour session, keep the core to **SCF + restart + DOS + two
exchange pairs**. Use the class comparison for numerical sensitivity, and
save parameter scans, surfaces, impurities, and transport for later.

## Instructor pilot notes (not acceptance tolerances)

These are measurements from a local GNU Fortran 14/OpenBLAS Release build of
the pinned `main` code, with two OpenMP threads, one OpenBLAS thread, MPI off,
and native CPU optimization off. They are not timings measured in GitHub
Codespaces, performance promises, or production reference values.

The complete command sequence, including the documented manual input edits,
was also exercised locally: SCF, restart, exchange, all three PNG plots, and
the result archive completed. The fresh-Codespace rehearsal remains an
instructor preparation step.

| Pilot | Observed result | Elapsed time |
| --- | --- | --- |
| Supplied bulk checkpoint, `conv_thr = 5e-8` | Converged in 1 step; residual about `4.0e-9` | 6.7 s |
| Original atom input, `beta = 0.1`, cap 20, `conv_thr = 5e-8` | Converged in 19 steps; residual `4.82e-8` | 30.7 s |
| Same original atom input, `beta = 0.01`, cap 20 | Not converged; final residual `6.7864e-6` | 32.0 s |
| Two-pair exchange using the 19-step checkpoint, before the optional restart | `J = 0.729460, 0.483801` mRy; `r/a = 0.866025, 1.000000` | 12.6 s |

The SCF pilot reported a spin moment about `2.103997 mu_B` and an orbital
moment about `0.041702 mu_B`. These are orientation values for this setup,
not accuracy targets. The pilot's peak resident memory was approximately
1.0 GB for SCF and 0.8 GB for exchange; the full Codespace also needs memory
for its editor, services, and build. The extra restart can change last digits.

The classroom workflow has not established convergence of exchange with
respect to recursion depth, energy sampling, real-space environment, or the
SCF threshold. It does not determine a full exchange-shell sum, spin-model
normalization, magnon dispersion, or Curie temperature.

### Sources within the pinned main branch

- `README.md`: build, examples, spin modes, and restart conventions.
- `example/bulk/bccFe/`: SCF inputs and supplied atom states.
- `tests/postproc/cases/exchange/bccFe/`: complete exchange input deck.
- `source/self.f90`, `source/mix.f90`: SCF criterion and printed residual.
- `source/symbolic_atom.f90`: checkpoint loading and saving.
- `source/bands.f90`: total DOS writing and Fermi-energy handling.
- `source/calculation.f90`: exchange post-processing workflow.
- `source/exchange.f90`: `jij.out` columns and energy scaling.
