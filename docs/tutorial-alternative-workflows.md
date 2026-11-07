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
- For MnTe, choose the phase and magnetic reference explicitly. Suggested:
  hexagonal NiAs-type alpha-MnTe with two opposite Mn sublattices, retaining
  Te and any justified ASA empty spheres. Supply a checked lattice and Mn/Te
  atom inputs, suitable sphere radii, and a converged reference. The checked
  RS-LMTO checkout contains no MnTe example or Te potential: this draft does
  not invent one or relabel another element's potential.
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

## 1. bcc Fe SCF (10–30 min)

Use the familiar bulk example, Broyden mixing, and a bounded SCF run. In this
`main_rc` checkout the example holds atom inputs rather than the old supplied
bulk checkpoint, so do not assume the old main-branch iteration count.

```bash
cd "$(git rev-parse --show-toplevel)"
cp -r example/bulk/bccFe runs/alternative-fe-scf
cd runs/alternative-fe-scf
python3 - <<'PYTHON'
from pathlib import Path
import f90nml
p = Path('input.nml')
f90nml.patch(p, {'self': {'nstep': 20, 'conv_thr': 5e-8},
                'mix': {'beta': 0.1, 'mixtype': 'broyden'}}, 'input.lesson.nml')
Path('input.lesson.nml').replace(p)
Path('input.initial.nml').write_text(p.read_text())
PYTHON
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

### Prepare and run independent SCF calculations

```bash
cd "$(git rev-parse --show-toplevel)"
python3 tools/tutorial/prepare_si.py
cd runs/alternative-si/rs
timeout 180 rslmto.x input.nml > run.log 2>&1
cd ../ks
timeout 180 rslmto.x input.nml > run.log 2>&1
```

`prepare_si.py` copies the checked two-site fixture into **independent**
folders and refuses to overwrite existing results. Both start from the same
atom inputs with `cold = .true.`; both electronic states are then solved
independently. Read the last residual, final energy, and moments for **each**
route. Confirm the forced-nonmagnetic Si moments remain negligible.

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

Use the complete exchange fixture and the converged Fe checkpoint. Retain the
same lattice, spin treatment, and recursion settings. In this code path,
exchange is post-processing of the saved electronic state, not a new SCF.

```bash
cd "$(git rev-parse --show-toplevel)"
cp -r tests/postproc/cases/exchange/bccFe runs/alternative-fe-jij
cp runs/alternative-fe-scf/Fe_scf.nml runs/alternative-fe-jij/Fe_out.nml
cd runs/alternative-fe-jij
python3 - <<'PYTHON'
from pathlib import Path
import f90nml
scf = f90nml.read('../alternative-fe-scf/input.nml')
p = Path('input.nml')
f90nml.patch(p, {
    'energy': {k: scf['energy'][k] for k in
               ('fermi', 'energy_min', 'energy_max', 'channels_ldos')},
    'control': {k: scf['control'][k] for k in ('nsp', 'recur', 'lld')},
    'self': {'nstep': 1},
}, 'input.lesson.nml')
Path('input.lesson.nml').replace(p)
PYTHON
timeout 180 rslmto.x input.nml > run.log 2>&1
cat jij.out
```

Check the actual pair displacements in the log; do not reuse cluster indices
for a changed geometry. In this `jij.out`, columns 1–2 are symbolic-type IDs,
3–5 are displacement/a, 6 is J in mRy, and 7 is distance/a. The two example
pairs are representatives of the first two bcc shells. They are not a
complete exchange model. `1 mRy = 13.605693 meV`.

## 4. MnTe: electronic sites versus magnetic sites (80–100 min)

### Reference and preparation contract

Use the instructor's alpha-MnTe calculation with a declared AFM state and
validated atom inputs. A minimal NiAs crystallographic sketch has Mn at
fractional `(0,0,0)` and `(0,0,1/2)`, and Te at `(1/3,2/3,1/4)` and
`(2/3,1/3,3/4)`. This is a structural sketch, **not a complete runnable ASA
deck**: lattice constants, radial potentials, sphere packing/empty spheres,
XC/+U choice, magnetic direction, and convergence still need to be supplied.
Do not convert fractional positions into `crd` without the lattice-vector
transformation; RS uses Cartesian coordinates scaled by `alat` there.

Use separate Mn labels for opposite magnetic sublattices and preserve their
mapping to the primitive sites. A noncollinear/no-SOC reference (`nsp = 3`)
can encode opposite in-plane moments; retain the instructor's actual reference
settings rather than changing the spin treatment in the exchange step.

**Keep all Mn, Te, and required empty-sphere sites in SCF and the Hamiltonian.**
Te contributes to bonding and Mn–Mn exchange even if it is not represented as
a spin in the later model. An explicit `force_nonmagnetic` constraint changes
the electronic calculation; it is not an exchange-pair filter. Do not impose
that constraint merely to remove Te from a spin model.

### Select only Mn endpoints

Copy the prepared full calculation to `runs/alternative-mnte-scf`, retaining
`lattice.nml`, `input.nml`, atom files, final checkpoints, report, and generated
`clust`. Inspect its residuals and Mn/Te moments. Verify the mapping between
input labels, `izp`, symbolic types, and generated cluster sites.

For a small demonstration, the tutorial helper selects nearby **Mn-only**
endpoints from that generated Bravais cluster. It does not truncate the
Hamiltonian. It assumes `morton_sfc = .false.` and the same geometry, replication,
and ordering in the SCF and subsequent exchange setup. Rebuild and inspect
indices if any of those change. Do not use the printed `clust` order with
Morton reordering enabled.

The following invocation is valid **only after verifying** the assumed mapping
`Mn types = 1,2`, `central Mn cluster sites = 1,2` in the supplied data:

```bash
cd "$(git rev-parse --show-toplevel)"
python3 tools/tutorial/select_mn_pairs.py \
    runs/alternative-mnte-scf/clust \
    --mn-types 1 2 --centres 1 2 --neighbours-per-centre 2 \
    --output runs/alternative-mnte-scf/mn-pairs.nml.inc
```

It prints endpoint types and distances and writes a snippet containing
`njij` and `ijpair`. Merge that snippet into the **existing** `&lattice` in a
separate exchange job; retain all electronic atom types and checkpoints.
Remove any old `njij`/`ijpair` entries first. The example selects two neighbours
per Mn centre, not a complete shell list; degeneracies and directed-pair
coverage must be considered for a quantitative model.

Set `pre_processing = 'none'`, `post_processing = 'exchange'`, and retain
compatible reference potentials, Fermi energy, geometry, nsp, and Hamiltonian
settings. Run one small job with the same three-minute cap. Inspect every
`jij.out` row: **both type indices must map to Mn**. Check that each Mn
sublattice appears as a centre. Never decide an element's identity merely
from the numerical size of its moment.

**Discussion:** Why can Te-mediated exchange survive after Te spins are
omitted? How does an AFM magnetic reference affect symmetry reduction? Why
must one avoid applying ferromagnetic or purely chemical symmetry reduction
without checking the magnetic state?

If a checked MnTe reference is not ready, use this slot as a worked input and
pair-selection demonstration with prepared results. An existing Mn3Sn example
is a possible later substitute, but its noncollinear magnetism adds theory;
it is not a drop-in MnTe replacement.

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

### Stage the one-site bcc primitive cell

From the repository root, prepare a separate export directory. This is the
bcc primitive-cell geometry in the code's lattice-parameter units; it is not
a 5984-site real-space cluster or a conventional two-site cubic cell.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir runs/alternative-fe-export
cat > runs/alternative-fe-export/lattice.nml <<'NML'
&lattice
nbulk_bulk = 1
ntot = 1
nbas = 1
nrec = 1
a(:, 1) = -0.5, 0.5, 0.5
a(:, 2) = 0.5, -0.5, 0.5
a(:, 3) = 0.5, 0.5, -0.5
crd(:, 1) = 0.0, 0.0, 0.0
/
NML
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
- A checked MnTe input/ASA reference, numerical Si convergence curves, full
  exchange-shell coverage, and an actual UppASD run remain preparation or
  follow-up work. No MnTe production result is claimed here.

## Wrap-up

Save the initial/final inputs, checkpoints, logs, DOS, exchange data, plots,
exported files, and both software versions. Download the result archive before
deleting the Codespace; stop it when finished. Ask each group to identify
which approximations were held fixed, what actually converged, and what must
be checked before their electronic calculation becomes a useful spin model.

References: the pinned RS-LMTO source/fixtures; RS-utils `0.11.0rc1` supplied
source, `docs/export_submodule.rst`, `docs/notes/UppASD_format.md`, and
`rsutil/export/uppasd.py`; [RS-utils manual](https://aibergman.github.io/RS-utils/).
