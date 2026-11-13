# RS-LMTO-ASA tutorial for ASP 2026

**Draft · `main_rc` · two-core Codespaces · Fe, Si, NiMnSb, Co2FeSi, and UppASD export**

This tutorial covers the following tasks:

1. Converge bcc Fe and investigate how recursion levels and cluster size affect its density of states (DOS).
2. Converge periodic Si with the Chebyshev method, investigate the effect of the number of polynomials, and compare its DOS with the reciprocal-space solver.
3. Calculate selected and automatically generated exchange interactions for bcc Fe, then export the results to UppASD with RS-utils.
4. Select Mn–Mn exchange pairs in half-Heusler NiMnSb and compare DOS results across calculation methods.
5. Generate RS-LMTO-ASA inputs from a Co2FeSi CIF file, automate exchange calculations, and export the results to UppASD.
6. Set up and converge a material of your choice from a CIF file.

[Launch the Codespace](https://codespaces.new/rslmtoasa/rslmtoasa/tree/tutorial/belem-26).

## Practical details

### Navigating files and viewing plots

Use the **Explorer** icon in the left sidebar to browse the repository. Expand a
folder to see its files; for example, `tutorials/bccFe` contains the supplied
Fe inputs, while the tasks create working copies under `runs/`. Click a text
file to open it in the editor. To locate the directory used by the terminal,
run `pwd`; use `ls` to list its contents and `cd ..` to move up one directory.
The command `cd "$(git rev-parse --show-toplevel)"` returns to the repository
root from any subdirectory.

After running `plot_dos.py`, expand the `plots` folder in the current run
directory and click a `.png` file to view it in the editor. `plot_exchange.py`
saves `exchange.png` in the current directory by default. If that image already
exists, choose a new filename with, for example,
`plot_exchange.py jijso.out --output exchange-second.png`.

### Editing inputs in the browser

Click a file in the Codespaces explorer to open its text editor. Change the
entries **inside their existing namelist sections**, then save with **Ctrl+S**
(**Cmd+S** on macOS). Add missing entries before that section's terminating
`/`. If a namelist section is missing, add it once. Do not duplicate an
existing `&self`, `&control`, `&lattice`, or `&reciprocal` section.

The lesson uses the editor for all input changes. Terminal commands copy
working files, run the solver, and invoke plotting or other tools.
After saving an edited input, preserve a copy before starting SCF:
`cp input.nml input.initial.nml`. Use `input.used.nml` for an exchange run. Commands that create working directories are intended to be run once.


## 1. Bcc Fe 

### A) Convergence to self-consistent solution

Use the supplied Fe deck in the `tutorials/bccFe` folder.

Open `lattice.nml` and `input.nml`. Inspect their contents and consider how the bcc Fe crystal is defined. The lattice
parameter remains `alat = 2.86120` Angstrom in `input.nml`; vectors and
positions in `lattice.nml` are scaled by `alat`. 

Keep `pre_processing = 'bravais'`: that selects the bulk SCF workflow,
while `crystal_sym = 'file'` selects the explicit geometry source. The same
`lattice.nml` will be used for the exchange calculation and UppASD export. To avoid overwriting the original setup, copy the files to a run folder as shown below:

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs
if [ ! -e runs/bccfe-scf ]; then
    cp -r tutorials/bccFe runs/bccfe-scf
fi
```

Open **`runs/bccfe-scf/input.nml`** in the editor. If needed, replace its
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
cd runs/bccfe-scf
rslmto.x input.nml > run.log 2>&1
```

Inspect `run.log`, `report.out`, and `Fe_out.nml`. Record the printed mixing
residual, total energy, spin/orbital moments, and Fermi energy. Only an explicit
`Converged!` at a residual below `5e-8` satisfies this teaching criterion;
reaching the iteration cap or returning exit code zero does not.

If needed, allow **one** further block using the saved checkpoint:

```bash
rslmto.x input.nml > restart.log 2>&1
```	
Check again whether the run has converged.

### B) Investigation of recursion levels and cluster size

Once the bcc Fe calculation has converged to an accepted residual, generate the density of states (DOS) plots by running

```
plot_dos.py
```
and inspect the PNG images in the `plots/` subdirectory. What does the total
DOS look like? What about the spin-projected DOS?

Next, copy the converged potential and input files to a new directory and move to that directory

```bash
mkdir ../fe-dos-comparison
cp input.nml lattice.nml ../fe-dos-comparison
cp Fe_out.nml ../fe-dos-comparison/Fe.nml
cd ../fe-dos-comparison

```

Then, open the `input.nml` in the new directory


Note: Make sure you have the correct `input.nml` open, i.e. not the one from the previous task. 

Now, change `nstep` to `1`, and change the number of recursion levels `lld` to another value (keep it between 15 and 50). Then run the `rslmto.x` binary again, and generate the new DOS files.

```bash
rslmto.x > out.log 2>&1
plot_dos.py
```

How does the new DOS differ from the previous one? Why?

If you want, you can test several different choices of `lld` or even the cluster size by changing the `rc` parameter and continue to investigate the resulting DOS. But remember to remove the output potential file `Fe_out.nml` between each run in this task. This ensures that all calculations start from the same, converged, potential in `Fe.nml`.


## 2. Si

### A) Convergence to self-consistent solution

Copy the supplied Si tutorial deck into a working directory and move into it.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs/silicon
cp -r tutorials/diamondSi runs/silicon/scf
cd runs/silicon/scf
```

Open **`runs/silicon/scf/input.nml`** and investigate the inputs. What looks
different from bcc Fe? The tutorial deck uses Chebyshev recursion with
`lld = 200` and `nstep = 20`.

This system uses the Chebyshev method instead of block recursion because Si is
non-metallic. Chebyshev calculations need more recursion levels (`lld`).

In the existing `&self` section, set `conv_thr = 5.0d-8` and keep `nstep = 20`
as the iteration cap. Save a copy of the input before running:

```bash
cp input.nml input.initial.nml
rslmto.x input.nml > run.log 2>&1
```

Inspect `run.log`, `report.out`, and the generated `Si1_out.nml` and
`Si2_out.nml` checkpoint files. As for Fe, accept the run only if the log explicitly says
`Converged!` with a residual below `5e-8`; reaching the iteration cap or
returning exit code zero is not sufficient. If it has not converged, allow one
further block from the saved checkpoint by rerunning:

```bash
rslmto.x input.nml > restart.log 2>&1
```

Check the convergence criterion again before continuing. How many iterations
were needed?


### B) Effect of number of Chebyshev polynomials

Follow the procedure you used for changing the number of recursion levels `lld` in the bcc Fe example: First, generate the original DOS with `plot_dos.py`, then make a new directory and copy your converged output there.

```bash
mkdir ../dos-comparison
cp input.nml lattice.nml ../dos-comparison
cp Si1_out.nml ../dos-comparison/Si1.nml
cp Si2_out.nml ../dos-comparison/Si2.nml
cd ../dos-comparison
```

In the copied `input.nml`, set `nstep = 1` and change `lld` to compare the
number of Chebyshev polynomials while holding the converged potentials fixed.
Investigate the effect `lld` has on the DOS for Si.

### C) Comparison with k-space DOS

Copy the supplied Si tutorial deck to a new directory:

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs/kspace-si
cp -r tutorials/diamondSi runs/kspace-si/scf
cd runs/kspace-si/scf
```

Open `input.nml` in the new directory. In `&self`, add
`use_kspace = .true.`, set `nstep = 20`, and set `conv_thr = 5.0d-8`. Add the
following `&reciprocal` section, or edit the existing section if one is already
present:

```fortran
&reciprocal
nk1 = 12
nk2 = 12
nk3 = 12
use_symmetry_reduction = .false.
use_time_reversal = .true.
use_shift = .false.
n_energy_points = 1000
dos_method = 'tetrahedron'
kspace_ham_order = 'first'
/
```

This fixes the mesh at 12 × 12 × 12 and uses the linear tetrahedron method. Save `input.nml` and perform the calculation.

```bash
rslmto.x input.nml > run.log 2>&1
```

Check for `Converged!` below the `5e-8` residual threshold, as in the Si
Chebyshev run. If needed, allow one restart block from the saved checkpoint:

```bash
rslmto.x input.nml > restart.log 2>&1
```

Once converged, generate the DOS with `plot_dos.py`.

How does the DOS differ from the Chebyshev DOS?

## 3. Fe exchange interactions

Now we can go back to bcc Fe and calculate exchange interactions between neighbouring spins. For this we need a new directory to run things in. Copy the **converged** bcc Fe files to the new directory.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir runs/bccfe-jij
cp runs/bccfe-scf/Fe_out.nml runs/bccfe-jij/
cp runs/bccfe-scf/input.nml runs/bccfe-jij/
cp runs/bccfe-scf/lattice.nml runs/bccfe-jij/
cp runs/bccfe-scf/report.out runs/bccfe-jij/
cp runs/bccfe-scf/fort.805 runs/bccfe-jij/clust.xyz
cd runs/bccfe-jij/
```

In addition to the input files, we also copied `fort.805` to `clust.xyz` for reasons that soon will be explained. We will also need the `report.out` file for later tasks.

In the `runs/bccfe-jij` directory, open `input.nml` and change the `&calculation` namelist to

```
&calculation
post_processing = 'exchange'
verbose = T
/
```

This sets the code to calculate exchange interactions. Next, decide which pairs of atoms to calculate `Jij` for by adding the following lines to the `&lattice` namelist.

```
nij = 2
ijpair(1, :) = 1, 2634
ijpair(2, :) = 1, 2635
```

These lines specify two interactions between the atoms listed on each `ijpair` line. The atoms are identified by their indices in the cluster. Therefore, use the same cluster-size settings in the SCF and exchange runs. You can find the atom indices in `clust.xyz`, which lists the coordinates of every atom in the cluster.

Save `input.nml` and run the code. 

Check the results in the `jijso.out` file. Here, columns 1–2 are symbolic-type IDs,
3–5 are displacement/a, 6 is J in mRy, and 7 is distance/a. 

### Plot the selected exchange rows

From the run directory, plot the actual exchange output:

```bash
plot_exchange.py jijso.out --csv jijso_selected.csv
```

The plot shows J in **mRy** versus distance/a, grouped by **directed
symbolic-type pair**. This command also writes `jijso_selected.csv`, which
retains the seven source columns, including J in **mRy** and Cartesian
displacement/a.

### Automated Jij-setup with `rsutils`

The helper package `rsutils` can generate exchange lists and inputs automatically from converged calculations. To try this, run the following command:

```bash
rsutil --exchange --cutoff 2.0 --unique
```

As the output says, this command results in the generation of a new folder `Atom_1`, ready for exchange interaction calculations. Try to run the new setup and plot the exchange.

```
cd Atom_1
rslmto.x | tee out.jij | grep atoms
plot_exchange.py jijso.out
```

### Automated export to UppASD

The `rsutil` package can also export `rslmto` calculations to a minimal UppASD input. This requires `input.nml`, `lattice.nml`, the exchange file `jijso.out`, and the log file `report.out`. The log is needed because UppASD also requires moment magnitudes, in addition to the structure and exchange data. Run:

```
cd ..
pwd
rsutil --export --jij Atom_1/jijso.out
```

This creates a minimal `UppASD` setup in the folder `uppasd_out`, unless another name is specified with the `--out OUTDIR` flag.

Inspect the resulting UppASD input. Does it look correct?

## Half-Heusler NiMnSb


### Selective exchange-pair selection

While the `rsutil` exchange setup tool is convenient, sometimes you need more control over pair selection. For example, the half-Heusler compound NiMnSb has only one strongly magnetic ion, Mn. Ni also has a small moment, but we can treat it as negligible here.

Here, we want to set up an exchange calculation for Mn–Mn pairs only, so we will select the pairs manually. Start by copying the supplied NiMnSb files to a run folder.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs/exchange-nimnsb
cp ./tutorials/MnNiSb/* runs/exchange-nimnsb
cd runs/exchange-nimnsb
```
Open `clust.xyz` in the new folder and identify relevant Mn–Mn couplings. Enter them in the existing `input.nml` as `ijpair` lines. Set `nij` to the number of exchange pairs you have defined. Also comment out the `pre_processing` line and add the appropriate `post_processing` line.

Run the calculation and plot the resulting exchange couplings with `plot_exchange.py`. Does it look as though you have captured the Mn–Mn couplings?

### Advanced: DOS features across methodologies

NiMnSb is half-metallic, meaning it has a gap in one spin channel. This can be difficult to resolve with block recursion, and the Chebyshev method requires many recursion levels.

To see how these different methods reconstruct the DOS, perform single-shot DOS calculations using the converged potentials. Edit `input.nml` by changing the values of `use_kspace`, `lld`, and `recur`. Set `nstep` to `1` to run a single iteration.

Make the changes in `input.nml` and run the code. After each run, generate the DOS with `plot_dos.py` and save the relevant images so they are not overwritten. What differences do you see across the methods?

## Heusler Co2FeSi

Another well-studied material in the Heusler family is `Co2FeSi`. Here, we will use the `rsutil` tools with minimal manual intervention, moving toward a high-throughput workflow.

To prepare the tasks, copy the reference data to a run directory and list the copied files.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs/ht-Co2Fesi
cp  ./tutorials/Co2FeSi/* runs/ht-Co2Fesi
cd runs/ht-Co2Fesi
ls .
```

### From .cif to RS-LMTO-ASA input

The run directory contains only a `.cif` file, a crystal-structure description. To run `rslmto.x`, first process the CIF file with `rsutil` using this command:

```
rsutil setup bulk --cif EntryWithCollCode221028.cif --empty_spheres
```

Open the generated `input.nml` and `lattice.nml` and inspect them. If everything looks good, run the code to converge the system. *Hint: If slow, you can run a `use_kspace` enabled run to speed things up.*

### Automatic exchange pair setup

Once the system has converged, use `rsutil` again to generate the relevant exchange couplings. Before setting them up, ensure that `use_kspace = .false.` in `input.nml`.

```
rsutil --exchange --cutoff 1.0 --unique
```

This time, the script generates four folders, one for each atom in the unit cell. Here, we use a limited cutoff.

### Running the exchange calculations

Now run `rslmto.x` in each directory to calculate the exchange interactions, then combine the results in a single output file.

```bash
for dir in Atom_?
do
cd $dir 
echo Processing $dir
rslmto.x | tee out.jij | grep atoms
cd ..
done
cat Atom_?/jijso.out > jij_all.txt
```

### Exporting to UppASD

With all exchange interactions prepared, we can run `rsutil` a final time to export our system to a minimal `UppASD` setup.

```
rsutil --export --jij ./jij_all.txt
```

Unless you specify an output directory, the generated files are placed in an `uppasd_out` subdirectory. Examine the outputs. Do they look as expected?

## Material X

Now you are ready to run the code on your own. Choose a CIF file from an online source, set up the calculation, and see whether you can converge the system.


### End of tutorial

References: the RS-LMTO tutorial decks in `tutorials/`; the bundled RS-utils
release-prep wheels ([package notes](../vendor/rsutils/README.md)); and the
[RS-utils manual](https://aibergman.github.io/RS-utils/).
