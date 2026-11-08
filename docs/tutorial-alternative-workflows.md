# Alternative two-hour tutorial: electronic structure to a spin model

**Draft · `main_rc` · two-core Codespaces · Fe, Si, MnTe, and UppASD export**

This alternative follows five tasks: (1) bcc Fe SCF, (2) periodic Si with
150-moment Chebyshev recursion and the built-in reciprocal-space solver,
(3) Fe exchange from task 1, (4) Mn-only exchange in multisublattice MnTe,
and (5) export to UppASD using RS-utils.

[Launch the alternative Codespace](https://codespaces.new/rslmtoasa/rslmtoasa/tree/tutorial/alternative-workflows).
The code basis is `main_rc` at `d11823829ac678e9c0a862c2e7099e54a8b4e11e`;
this branch adds the environment and tutorial files. The original `main`-based
tutorial remains separate. **KS below means reciprocal/k-space evaluation of
the same RS-LMTO model**, not a second DFT program or different functional.

## What fits in two hours

| Time | Exercise | Student deliverable |
| --- | --- | --- |
| 0–10 min | Open the prepared Codespace; discuss approximations | Working executable and provenance |
| 10–30 min | 1. bcc Fe SCF and restart | Residual, moment, and checkpoint |
| 30–65 min | 2. Si Chebyshev/PBC and k-space Gaussian DOS | Two convergence decisions and DOS overlay |
| 65–80 min | 3. Fe exchange from the SCF state | Two representative Fe–Fe couplings |
| 80–100 min | 4. Multisublattice MnTe, selecting Mn–Mn pairs | Audited magnetic pair list and a small exchange job |
| 100–115 min | 5. RS-utils → UppASD | Reviewed position, moment, and exchange files |
| 115–120 min | Save results and discuss next steps | Input/output archive |

Use one calculation at a time, two OpenMP threads, one OpenBLAS thread, and
no GPU or MPI. Each solver command below has a three-minute cap. Have fallback
results available; do not let repeated unconverged runs consume the session.

**Preparation is essential.** Students can converge the lean Si models live
when timings permit. Converging MnTe from a new atomic guess, checking ASA
sphere choices, and establishing numerical convergence of all observables
belong to instructor preparation. A two-hour session cannot promise fresh,
production-converged calculations for all five tasks.

## Instructor preparation

- Launch and build Codespaces before class. Test the exact student account
  workflow, free quota, internet access, and per-run timings.
- Stage matched Si and MnTe reference states and fallback results, with code
  revision, full inputs, SCF residuals, energy meshes, and output provenance.
  Converged SCF and converged numerical settings are separate requirements.
- Use the instructor-supplied MnTe lattice and four checkpoints in
  `tutorials/mnte/`. It retains both Mn and both Te sites, with opposite Mn
  starting directions along x and `nsp = 4` (SOC). The local full-size replay met the SCF criterion on the pinned code.
  The reduced-size replay did not; prepare a converged reduced-size state
  before class if the exercise requires self-consistent smaller-cell J values.
- Install the instructor-provided RS-utils package before class, rather than
  building it on every student's first launch. The uploaded source identifies
  itself as `rsutil 0.11.0rc1`. Its source repository is private; it is not copied
  into this public tutorial branch. Students need an instructor-authorized
  archive/wheel or repository access. No access tokens belong in this repository.
- Rehearse the complete export with the chosen UppASD version. The supplied
  exporter is a starter converter, not a substitute for checking the target
  spin Hamiltonian, geometry, units, and simulation settings.

The environment includes NumPy, Matplotlib, gnuplot, and `f90nml`. Additional
RS-utils dependencies are installed in a virtual environment. From the
repository root, after placing the approved source at `work/RS-utils-main`:

```bash
python3 -m venv --system-site-packages .venv-rsutils
source .venv-rsutils/bin/activate
python -m pip install ./work/RS-utils-main
rsutil --help
```

The source installation builds native Fortran utilities and installs its
runtime dependencies. The existing compiler and CMake are used. An appropriate
prebuilt wheel avoids the compile step. Do not assume a public PyPI package
with the same name is the instructor's version.

## 0. Start and identify the environment

Run from the repository root in a fresh terminal:

```bash
cd "$(git rev-parse --show-toplevel)"
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
command -v rslmto.x
python3 -c 'import numpy, matplotlib, f90nml; print("Python tools ready")'
git rev-parse HEAD
mkdir -p runs
```

Keep the starting inputs as well as the final ones: SCF updates the Fermi
energy in `input.nml`. Commands assume the new `runs/alternative-*` directories
do not exist. Choose fresh names when repeating the lesson.

### Editing inputs in the browser

Click a file in the Codespaces explorer to open its text editor. Change the
entries **inside their existing namelist sections**, then save with **Ctrl+S**
(**Cmd+S** on macOS). Add missing entries before that section's terminating
`/`. Do not append a second `&self`, `&control`, or `&lattice` section.

The lesson uses the editor for all input changes. Terminal commands copy
working files, run the solver, and invoke plotting or pair-selection tools.
After saving an edited input, preserve a copy before starting SCF:
`cp input.nml input.initial.nml`. Use `input.used.nml` for an exchange run. Commands that create working directories are intended to be run once.

## 1. bcc Fe SCF (10–30 min)

Use the supplied explicit-lattice Fe deck, Broyden mixing, and a bounded
SCF run. This deck holds atom inputs rather than the old supplied bulk
checkpoint, so do not assume the old main-branch iteration count.

Both Fe starting inputs use **`crystal_sym = 'file'`**. Open
`lattice.nml` beside `input.nml`: it contains the one-site primitive bcc
vectors, an Fe site at the origin, and the site/type mapping. The lattice
parameter remains `alat = 2.86120` Angstrom in `input.nml`; vectors and
positions in `lattice.nml` are scaled by `alat`. This primitive cell has
one atom and volume `a^3/2`, not two atoms in a conventional cubic box.

Keep `pre_processing = 'bravais'`: that selects the bulk SCF workflow,
while `crystal_sym = 'file'` selects the explicit geometry source. The same
`lattice.nml` will be carried into exchange and UppASD export.

```bash
cd "$(git rev-parse --show-toplevel)"
cp -r tutorials/bccFe runs/alternative-fe-scf
```

Open **`runs/alternative-fe-scf/input.nml`** in the editor. Replace its
existing `&self` and `&mix` sections with:

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

Save, then run in the terminal:

```bash
cd runs/alternative-fe-scf
cp input.nml input.initial.nml
timeout 180 rslmto.x input.nml > run.log 2>&1
```

Inspect `run.log`, `report.out`, and `Fe_scf.nml`. Record the printed mixing
residual, total energy, spin/orbital moments, and Fermi energy. Only an explicit
`Converged!` at a residual below `5e-8` satisfies this teaching criterion;
reaching the iteration cap or returning exit code zero does not.

If needed, allow **one** further block using the saved checkpoint:

```bash
cp Fe_scf.nml Fe_out.nml
timeout 180 rslmto.x input.nml > restart.log 2>&1
```

Use a matched instructor checkpoint if this still does not converge. Do not
pass an unconverged state into task 3. The threshold is a classroom criterion,
not a statement of exchange-parameter accuracy. `nsp = 2` retains the supplied
collinear SOC treatment; `nsp = 1` would not mean nonmagnetic.

## 2. Silicon: RS versus k-space (30–65 min)

### Same physical model, two numerical routes

We use the existing two-Si diamond primitive cell and a deliberately small
**sp basis (`lmax = 1`)**. This is a teaching model, not a quantitative claim
about the experimental Si band gap. The fixture's original one-step checks
are not converged-material benchmarks; the lesson increases the SCF cap.

| Setting | RS | K-space |
| --- | --- | --- |
| Geometry, basis, atoms, XC, screened structure constants | Identical | Identical |
| `nsp` | 1; both Si types forced nonmagnetic | Same |
| Hamiltonian order | `hoh = .false.` | `kspace_ham_order = 'first'`, same `hoh` |
| SCF selector | `use_kspace = .false.` | `use_kspace = .true.` |
| Spectral solver | `recur = 'chebyshev'`, **`lld = 150`** | Built-in eigenvalue workflow; `lld` does not set k-space accuracy |
| Sampling | PBC, `n1 = n2 = n3 = 12` | Full Gamma-centred `8 × 8 × 8` mesh |
| DOS smoothing | Chebyshev reconstruction/kernel | Gaussian, `gaussian_sigma = 0.02` Ry |
| SCF cap / threshold | 25 / `5e-8` | Same |

The input key is **`lld`**, rather than `ll`. PBC replication counts and a
k-point mesh are different controls. A 150-moment reconstruction has finite,
energy-dependent resolution; a Gaussian is not the same kernel. The starting
sigma `0.02 Ry ≈ 0.272 eV` gives a readable classroom comparison without
claiming an exact kernel match. It must not be tuned solely to hide a difference.

### Prepare and edit independent SCF inputs

Copy the same starting fixture into two separate working directories:

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs/alternative-si
cp -r tests/scf/cases/bulk/diamondSi runs/alternative-si/rs
cp -r tests/scf/cases/bulk/diamondSi runs/alternative-si/ks
```

Open **`runs/alternative-si/rs/input.nml`**. Set the following entries in the
existing sections; retain the atom labels, geometry, starting Fermi energy,
and any entries not mentioned here.

| Existing section | Entries to set or add |
| --- | --- |
| `&lattice` | `pbc = .true.`, `morton_sfc = .false.`, `n1 = 12`, `n2 = 12`, `n3 = 12`, `strux_backend = 'strux_lib'` |
| `&self` | `nstep = 25`, `conv_thr = 5.0d-8`, `cold = .true.`, `use_kspace = .false.`, `force_nonmagnetic(:) = 2*.true.` |
| `&control` | `nsp = 1`, `lmax = 1`, `lld = 150`, `recur = 'chebyshev'`, `cheb_backend = 'legacy'` |
| `&energy` | `energy_min = -1.5`, `energy_max = 1.0`, `channels_ldos = 1000` |
| `&mix` | `beta = 0.1`, `mixtype = 'broyden'` |
| `&hamiltonian` | `hoh = .false.` |

Each comma-separated item in the table is a separate namelist assignment;
placing one assignment per line is easiest to read. Keep the existing
periodic-direction flags `b1`, `b2`, and `b3` set to `.true.`. Save the file.

In the editor, select and copy the entire edited RS input, then open
**`runs/alternative-si/ks/input.nml`**, select its entire contents, and paste.
This copies the same starting settings before either run changes its input.
In the KS file, change **only** `use_kspace` to `.true.` inside `&self`, and
add this new section after the other namelists:

```fortran
&reciprocal
nk1 = 8
nk2 = 8
nk3 = 8
use_symmetry_reduction = .false.
use_time_reversal = .true.
use_shift = .false.
n_energy_points = 1000
dos_method = 'gaussian'
gaussian_sigma = 0.02
kspace_ham_order = 'first'
/
```

Save both files. If a `&reciprocal` section already exists in a reused input,
edit it instead of adding a duplicate. Both routes use the same starting
atom inputs with `cold = .true.` and are solved independently.

Run from the repository root:

```bash
cd "$(git rev-parse --show-toplevel)"
cd runs/alternative-si/rs
cp input.nml input.initial.nml
timeout 180 rslmto.x input.nml > run.log 2>&1
cd ../ks
cp input.nml input.initial.nml
timeout 180 rslmto.x input.nml > run.log 2>&1
```

Read the last residual, final energy, and moments for **each** route. Confirm
the forced-nonmagnetic Si moments remain negligible. The preparation helper
`prepare_si.py` remains available for instructors, but students do not need
it for this editor-based workflow.

If a cap is reached, compare the saved state with the instructor reference.
Do not label the result converged because a DOS file exists. To continue a
run in a controlled second block, copy each site's `_scf.nml` to its
`_out.nml`, set `cold = .false.`, retain that route's other settings, and
save its next log under a new name. Avoid repeatedly cold-starting.

### Overlay `totaldos.out`

```bash
cd "$(git rev-parse --show-toplevel)"
python3 tools/tutorial/compare_totaldos.py \
    runs/alternative-si/rs/totaldos.out \
    runs/alternative-si/ks/totaldos.out \
    --rs-energy-min -1.5 --ks-energy-min -1.5 \
    --output runs/alternative-si/dos-comparison.png
```

Both SCF writers print `E - E_F` in Ry, but their Fermi estimates need not
coincide; an insulator's Fermi position inside the gap is not a unique band
alignment reference. The command restores the **absolute grid** using the
known first energy point, `-1.5 Ry`, for each of these generated inputs. This
avoids relying on a possibly different iteration's logged Fermi value. Do
not use these two options for another grid without checking its first point.

The script also accepts explicit `--rs-fermi` and `--ks-fermi` (in Ry), or
header values when available. For a separate plot relative to each route's
own Fermi zero, use `--alignment fermi` and label that choice. The script does
not fit a shift, smooth the data, rescale areas, or set a pass/fail tolerance.

It plots states/(eV atom), dividing the two-site cell DOS by two and converting
states/Ry to states/eV. It prints full sampled-window integrals, minimum DOS,
and an RMS difference over the common energy range. For this spin-summed,
two-site sp model, a window covering the entire represented spectrum contains
**16 states/cell**, while the valence occupation is **8 electrons/cell**.
These are different counts. Confirm the window really covers the spectrum;
do not normalize a truncated DOS to 16 by hand.

**Questions:** Are peak positions and integrated weights consistent? What
comes from smoothing, sampling, or independently converged potentials? Does
a common small SCF residual imply identical DOS? Would separate Fermi zeros
misrepresent a comparison of two insulating spectra?

### Convergence: what is and is not done in class

The lesson converges each SCF state at its specified numerical settings.
**Numerical convergence of the DOS is a further task.** Before class, the
instructor should provide small comparisons such as:

- RS: compare 150 versus 200 moments, and a larger PBC replication, keeping
  the structure, basis, Hamiltonian, and energy window fixed.
- KS: compare `8^3`, `12^3`, and if necessary denser meshes at fixed sigma.
  Reduce sigma only after the mesh supports the finer spectral resolution.
- Separate a shared-potential solver comparison from the final comparison of
  independently converged states. Use the code's existing DOS-equivalence
  regression as a guide, without changing its references or tolerances.

Give each group at most one extra numerical setting, or use the prepared
curves. Agree on an observable, energy window, and acceptable change before
claiming numerical convergence. The two-hour exercise does not establish a
production-converged gap or prove equality of the two broadening kernels.

## 3. Fe exchange from task 1 (65–80 min)

Use the explicit-lattice exchange starting input and the converged Fe
checkpoint. Retain the same lattice, spin treatment, and recursion settings. In this code path,
exchange is post-processing of the saved electronic state, not a new SCF.

```bash
cd "$(git rev-parse --show-toplevel)"
cp -r tutorials/bccFe runs/alternative-fe-jij
cp runs/alternative-fe-jij/input.exchange.nml runs/alternative-fe-jij/input.nml
cp runs/alternative-fe-scf/lattice.nml runs/alternative-fe-jij/lattice.nml
cp runs/alternative-fe-scf/Fe_scf.nml runs/alternative-fe-jij/Fe_out.nml
```

Open **`runs/alternative-fe-jij/input.nml`** alongside the **final**
`runs/alternative-fe-scf/input.nml` from task 1. In the exchange input:

- In `&energy`, copy the final SCF values of `fermi`, `energy_min`,
  `energy_max`, and `channels_ldos`. Copy the actual final Fermi value rather
  than a number printed in this handout.
- In `&control`, match the SCF values of `nsp`, `recur`, and `lld`.
- In `&self`, set `nstep = 1`.
- Retain `post_processing = 'exchange'`, **`crystal_sym = 'file'`**, the
  copied `lattice.nml`, and the two supplied `ijpair` entries.

Save the exchange input, then run:

```bash
cd runs/alternative-fe-jij
cp input.nml input.used.nml
timeout 180 rslmto.x input.nml > run.log 2>&1
cat jij.out
```

Check the actual pair displacements in the log; do not reuse cluster indices
for a changed geometry. In this `jij.out`, columns 1–2 are symbolic-type IDs,
3–5 are displacement/a, 6 is J in mRy, and 7 is distance/a. The two example
pairs are representatives of the first two bcc shells. They are not a
complete exchange model. `1 mRy = 13.605693 meV`.

## 4. MnTe: electronic sites versus magnetic sites (80–100 min)

### Supplied reference and electronic/magnetic site mapping

The instructor's working deck is included in `tutorials/mnte/`, retaining
its original crystallographic attribution and checkpoint parameters.

| Electronic type / central cluster site | Label | Atomic number | Starting direction |
| --- | --- | --- | --- |
| 1 | Mn1 | 25 | +x |
| 2 | Mn2 | 25 | -x |
| 3 | Te1 | 52 | +z orientation marker |
| 4 | Te2 | 52 | +z orientation marker |

Here `izp = 1,2,3,4`, and the four primitive sites are the first four entries
of the generated cluster. **The `mom` checkpoint vector is a direction, not
its moment magnitude.** In particular, a unit-length Te direction does not
mean Te carries a 1 mu_B moment. Inspect actual magnitudes in `report.out`.

Retain the supplied `alat = 4.1501998901` Angstrom, `wav`, sphere radii,
lattice, energy interval/sampling, Chebyshev `lld = 200`, and `hoh = .false.`. The teaching copy uses
**`11x11x7` PBC (3,388 sites)** instead of the supplied `15x15x15`
(13,500 sites), reducing the site count by about 75%. Its c vector is
longer than its basal vectors, so the reduced physical box lengths are
roughly balanced. This changes finite-size sampling; it is a classroom
profile, not an established numerical-convergence setting.
The MnTe setting is **`nsp = 4`, with SOC**, rather than
the no-SOC example sketched in the earlier draft. Do not change that spin
treatment between the reference replay and exchange.

**Keep all four sites and all four checkpoints in the electronic Hamiltonian.**
Te contributes to bonding and Mn–Mn exchange even when it is omitted from
the later spin model. `force_nonmagnetic` changes the electronic calculation;
it is not an exchange-pair filter and is not added by this lesson.

### Bounded checkpoint replay

Copy the tracked starting files to a working directory. The original input
allows 100 SCF steps and contains `njij = 162` without an explicit pair list.
For this SCF check, use one step and no exchange pairs; do not inherit that
unused count into an exchange calculation.

```bash
cd "$(git rev-parse --show-toplevel)"
cp -r tutorials/mnte runs/alternative-mnte-scf
```

Open **`runs/alternative-mnte-scf/input.nml`**. In its existing `&lattice`
section, change or add:

```fortran
njij = 0
morton_sfc = .false.
n1 = 11
n2 = 11
n3 = 7
```

In its existing `&self` section, change `nstep = 100` to **`nstep = 1`**.
Retain all other settings, including `nsp = 4`, `lld = 200`, and the energy
and atom sections. Save, then run:

```bash
cd runs/alternative-mnte-scf
cp input.nml input.initial.nml
timeout 180 rslmto.x input.nml > run.log 2>&1
```

The replay creates `clust` and, if completed, a new report and SCF states.
Inspect the residual, opposite Mn moments, and Te magnitudes. Do not call a
one-step replay converged unless the log says it met the threshold. Use the
instructor's checked reference or prepared output when the wall-time cap is
reached. This does not require students to converge MnTe from atomic guesses.

### Select only Mn endpoints

In the reduced teaching geometry, `clust` contains 3,388 sites. The first two central
sites are Mn; types 3 and 4 are Te. Keep `morton_sfc = .false.`: the helper
uses the printed Bravais cluster order, which precedes optional Morton
reordering. Retain this same geometry, replication, and ordering in exchange.

```bash
cd "$(git rev-parse --show-toplevel)"
python3 tools/tutorial/select_mn_pairs.py \
    runs/alternative-mnte-scf/clust \
    --mn-types 1 2 --centres 1 2 --neighbours-per-centre 2 \
    --output runs/alternative-mnte-scf/mn-pairs.nml.inc
```

For the reduced teaching copy of this deck the four selected directed pairs are:

| Centre | Neighbour | Endpoint types | Distance/a |
| --- | --- | --- | --- |
| 1 | 2 | Mn1–Mn2 | 0.80951520 |
| 1 | 1273 | Mn1–Mn2 | 0.80951520 |
| 2 | 1 | Mn2–Mn1 | 0.80951520 |
| 2 | 428 | Mn2–Mn1 | 0.80951521 |

Indices are valid only for this generated cluster and settings. The last
printed distance differs by rounding of stored coordinates. This selects
two nearest neighbours per centre; it is not a complete shell list, and
includes directed counterparts. Never count them as independent undirected
bonds without checking the spin-model convention.

### Stage the exchange job without dropping Te

The next commands copy the reduced working state and its matching input into
an exchange directory. Only use it as a converged smaller-cell result after
checking its SCF criterion. A one-step nonconverged replay can support an
explicitly labelled **frozen-state demonstration**, not a converged exchange
reference. Instructor preparation should provide a smaller-cell checkpoint
if that distinction would otherwise consume the lesson.

The pair snippet replaces the unused `njij = 162` with four explicit Mn-only
pairs. The staged job keeps all four electronic types and matching energy
settings; it does not drop Te or revert to the original 15x15x15 geometry.

```bash
cd "$(git rev-parse --show-toplevel)"
cp -r runs/alternative-mnte-scf runs/alternative-mnte-jij
```

Open **`runs/alternative-mnte-jij/input.nml`**. Replace the existing
`&calculation` section with:

```fortran
&calculation
pre_processing = 'none'
post_processing = 'exchange'
verbose = T
/
```

Open **`runs/alternative-mnte-scf/mn-pairs.nml.inc`** in another editor tab.
Copy its `njij` and `ijpair` assignments into the exchange input's **existing
`&lattice` section**, replacing its old `njij` assignment and any old
`ijpair` assignments. Do not replace the whole lattice section: keep all its
geometry settings. For the reduced teaching profile, the assignments are:

```fortran
njij = 4
ijpair(1, :) = 1, 2
ijpair(2, :) = 1, 1273
ijpair(3, :) = 2, 1
ijpair(4, :) = 2, 428
```

Use the helper's newly generated indices if the geometry has changed.
Retain `morton_sfc = .false.`, the `11x11x7` replication, `nstep = 1`, all
four atom labels, and the energy/Fermi settings copied with the working
state. Save, then run:

```bash
cd runs/alternative-mnte-jij
cp input.nml input.used.nml
timeout 180 rslmto.x input.nml > run.log 2>&1
```

This is a frozen-reference exchange exercise. Record whether its smaller-cell
reference met the selected SCF threshold. To use a newly converged MnTe state,
transfer **all four** final checkpoints and the same run's final energy/Fermi
settings, retaining their labels and geometry.

Check the exchange log's pair indices/displacements. Every `jij.out` row
must have types `(1,2)` or `(2,1)`, and both Mn centres must be represented.
The electronic input still contains types 1–4 and loads both Te checkpoints.
There should be no Te exchange endpoints. Do not decide element identity
from the magnitude of its moment.

**Discussion:** Why does Te-mediated exchange remain after Te spins are
omitted? Why does AFM/SOC symmetry require more care than chemical symmetry?
How would you extend the pair list to in-plane Mn–Mn shells without changing
the electronic reference? The next shells require new pair selection and a
convergence study, not changing Te into a magnetic atom.

## 5. Export to UppASD using RS-utils (100–115 min)

Use **Fe** for the first export. Its one-site magnetic primitive cell avoids
mixing multisublattice projection with a first introduction to the file
formats. Keep the MnTe projection as the extension described below.

The supplied RS-utils source supports:

```bash
rsutil --export --data /path/to/staged-fe \
    --jij /absolute/path/to/jij.out \
    --report /absolute/path/to/fe-scf/report.out \
    --out uppasd_out
```

Here `--out` is relative to `--data`; `--jij` accepts the first six numeric
columns in the supplied exchange file, so select `jij.out` deliberately.
Do not silently substitute another perturbation-order output. The CLI reads
`lattice.nml` and obtains moment magnitudes from `report.out`.

### Reuse the Fe calculation's lattice file

Copy the **same `lattice.nml` used in task 1** into a separate export
directory. No geometry reconstruction or new lattice file is needed. It
represents the one-site primitive bcc cell, not the generated real-space
cluster or a conventional two-site cubic cell.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir runs/alternative-fe-export
cp runs/alternative-fe-scf/lattice.nml runs/alternative-fe-export/lattice.nml
```

Then run the exporter from the repository root:

```bash
rsutil --export --data runs/alternative-fe-export \
    --jij "$(pwd)/runs/alternative-fe-jij/jij.out" \
    --report "$(pwd)/runs/alternative-fe-scf/report.out" \
    --out uppasd_out
```

The exported files are `inpsd.dat`, `posfile`, `momfile`, and `jfile`. In the
uploaded version, moments come from the SCF report, but **directions default
to +x**. The lattice-only CLI path uses a single default type and does not
apply `input.nml`'s `alat`. Coordinates and cell are in a common `a` unit;
configure the intended physical lattice parameter in the target simulation
where required. Do not interpret the raw cell as Å without that conversion.

### Inspect the output before simulating

- `posfile`: one bcc primitive-basis spin, no spurious cluster sites.
- `momfile`: the Fe spin-moment **magnitude** from task 1, not a fallback 1 mu_B;
  the default direction is +x, not inherited from the SCF orientation.
- `jfile`: displacement and J values correspond to the chosen input rows,
  with J retained in mRy. Check sign, pair counting, coordinate units, and the
  unit-vector versus spin-length convention for your UppASD Hamiltonian.
- `inpsd.dat`: review its cell, boundary conditions, `Sym`, `ncell`, and file
  references. The starter output is not a complete production dynamics run.

**The exporter copies the exchange rows; it does not expand the two Fe
representatives into a full neighbour table.** Use an instructor-provided,
adequately converged exchange-shell dataset or a separately verified symmetry
expansion before meaningful UppASD dynamics. Do not derive a Curie temperature
from this two-pair teaching export. An UppASD executable or dynamics run is
not required for the core export exercise.

### MnTe export extension

Project the full electronic reference onto a **Mn-only magnetic basis** while
retaining the original lattice vectors and an explicit old-site → new-spin
map. Preserve distinct Mn sublattice/site identities even though both are Mn.
Remap both exchange endpoints and the displacement/translation consistently;
copy the AFM moment directions explicitly. Te and empty spheres are absent
from the spin model, not from the electronic state used to obtain its J values.

The supplied starter CLI has no Mn-only selection option, collapses its
lattice-only path to one default type, and writes +x moment directions.
Therefore it must not be presented as a one-command, AFM-preserving MnTe
exporter. The first session exports Fe; a reviewed Mn-only projection and
AFM-preserving export are a follow-up step. Do not silently export Te with
a fabricated default moment or merge the two Mn sublattices.

## Instructor pilot evidence and remaining work

Local pilots use the pinned RC executable, GNU Fortran 14/OpenBLAS, two
OpenMP threads, one OpenBLAS thread, and no GPU/MPI. They are not fresh
Codespaces timings or production-accuracy claims.

- The generated Si inputs (150 moments/PBC for RS, `8^3` Gaussian sigma 0.02 Ry
  for KS) reached the `5e-8` SCF threshold in 13 iterations each. Printed final
  residuals were `1.56e-8` (RS) and `4.7e-9` (KS).
- The two total DOS files integrated to approximately `16.000006` and
  `15.999993` states/cell. The common-grid relative RMS difference was about
  `0.0969`, without fitted shifts, renormalization, or extra smoothing.
- The plot script was also checked on independently specified triangular
  spectra with different Fermi offsets: each integrated to 0.5 states/cell
  and the absolute-energy RMS difference was zero after origin restoration.
- The magnetic-pair helper excluded deliberately closer Te-type sites in a
  small synthetic cluster and selected only the two declared Mn endpoints.
  This checks indexing/filtering, not MnTe material physics.
- The RC Fe pilot converged after a second bounded SCF block (final printed
  residual `1.28e-8`). Its two-pair exchange output contained `0.731162` and
  `0.483901` mRy. These are setup-specific observations, not reference targets.
- The explicit `crystal_sym = 'file'` Fe deck reproduced the built-in bcc
  preset's 5,984 cluster positions and types. A checkpoint replay met the
  teaching threshold (printed residual `1.16e-8`), and the exchange rows
  matched the built-in pilot at printed precision. The same `lattice.nml`
  was used for the SCF, exchange, and one-site UppASD export.
- The uploaded RS-utils Python CLI path wrote all four Fe export files. The
  moment magnitude was `2.104052 mu_B`; the two exchange rows retained the
  source values and displacement components. Those equalities hold by
  construction of the serialization and do not test the physical model.
- The native RS-utils wheel build failed in this local extracted toolchain
  when a static, non-PIC `libgfortran.a` was linked into its shared library.
  Export was exercised through `python -m rsutil.cli` with the supplied
  source on `PYTHONPATH` and its Python runtime dependencies installed; it
  did not use the native setup library. A normal Codespace wheel/source
  installation still needs the instructor's rehearsal.
- The uploaded full-size MnTe state (`15x15x15`, 13,500 sites) met the SCF
  criterion in one replayed step, residual `1.7e-9`, taking about 170 s locally.
  Printed Mn spin magnitudes were `4.126320 mu_B`; Te magnitudes printed as
  `0.000000 mu_B`. The checkpoint direction vectors are not those magnitudes.
- Reducing PBC to `11x11x7` (3,388 sites) lowered the one-step replay time to
  about 51 s. Its residual was `1.34003e-5`, so the smaller state was **not**
  SCF converged. Its printed Mn magnitudes were `4.126274 mu_B` and its
  Fermi energy was `-0.224518 Ry`.
- The four-pair reduced-size exchange run completed in about 147 s, using
  the updated one-step state and its matched Fermi energy. All rows were
  Mn1–Mn2 or Mn2–Mn1, with J about `-2.157862 mRy`. This is a workflow and
  frozen-state teaching observation, not a converged-material reference.
- Numerical Si convergence curves, a converged **reduced-size** MnTe reference,
  full exchange-shell coverage, and an actual UppASD run remain preparation
  or follow-up work. No MnTe production result is claimed here.

## Wrap-up

Save the initial/final inputs, checkpoints, logs, DOS, exchange data, plots,
exported files, and both software versions. Download the result archive before
deleting the Codespace; stop it when finished. Ask each group to identify
which approximations were held fixed, what actually converged, and what must
be checked before their electronic calculation becomes a useful spin model.

References: the pinned RS-LMTO source/fixtures; RS-utils `0.11.0rc1` supplied
source, `docs/export_submodule.rst`, `docs/notes/UppASD_format.md`, and
`rsutil/export/uppasd.py`; [RS-utils manual](https://aibergman.github.io/RS-utils/).
