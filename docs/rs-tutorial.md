# RS-LMTO-ASA tutorial for ASP 2026

**Draft · `main_rc` · two-core Codespaces · Fe, Si, MnTe, and UppASD export**

This tutorial contains the following tasks: 
(1) bcc Fe SCF and 
(2) periodic Si with
Chebyshev recursion and the reciprocal-space solver,
(3) Fe exchange from task 1, 
(4) Exchange in multisublattice MnTe,
(5) Export calculation to UppASD using RS-utils.
(6) Import your own .cif file to RS-LMTO-ASA

[Launch the Codespace](https://codespaces.new/rslmtoasa/rslmtoasa/tree/tutorial/alternative-workflows).

## Practical details

### Editing inputs in the browser

Click a file in the Codespaces explorer to open its text editor. Change the
entries **inside their existing namelist sections**, then save with **Ctrl+S**
(**Cmd+S** on macOS). Add missing entries before that section's terminating
`/`. Do not append a second `&self`, `&control`, or `&lattice` section.

The lesson uses the editor for all input changes. Terminal commands copy
working files, run the solver, and invoke plotting or other tools.
After saving an edited input, preserve a copy before starting SCF:
`cp input.nml input.initial.nml`. Use `input.used.nml` for an exchange run. Commands that create working directories are intended to be run once.


## 1. Bcc Fe 

### A) Convergence to self-consistent solution

Use the supplied Fe deck in the `tutorials/bccFe`folder. 

Open
`lattice.nml` and `input.nml`: Inspect their contents and try to understand how the bcc Fe crystal is defined. The lattice
parameter remains `alat = 2.86120` Angstrom in `input.nml`; vectors and
positions in `lattice.nml` are scaled by `alat`. 

Keep `pre_processing = 'bravais'`: that selects the bulk SCF workflow,
while `crystal_sym = 'file'` selects the explicit geometry source. The same
`lattice.nml` will be carried into exchange and UppASD export. In order to not overwrite the original setup, copy the files to a run folder as per below:

```bash
cd "$(git rev-parse --show-toplevel)"
cp -r tutorials/bccFe runs/alternative-fe-scf
```

Open **`runs/alternative-fe-scf/input.nml`** in the editor. If needed, replace its
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
rslmto.x input.nml > run.log 2>&1
```

Inspect `run.log`, `report.out`, and `Fe_scf.nml`. Record the printed mixing
residual, total energy, spin/orbital moments, and Fermi energy. Only an explicit
`Converged!` at a residual below `5e-8` satisfies this teaching criterion;
reaching the iteration cap or returning exit code zero does not.

If needed, allow **one** further block using the saved checkpoint:

```bash
rslmto.x input.nml > restart.log 2>&1
```	
Check again if the run has converged. 

### B) Investigation of recursion levels and cluster size

Once the bcc Fe calculation has converged to an accepted norm, generate the density of states (DOS) plots by running

```
plot_dos.py
```
and investigate the images that has been created in the current folder. How doest the total DOS look, and how about the spin-projected DOS?

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
rslmto.x | out.log
plot_dos.py
```

How does the new DOS differ from the previous one? Why?

If you want, you can test several different choices of `lld` or even the cluster size by changing the `rc` parameter and continue to investigate the resulting DOS. But remember to remove the output potential file `Fe_out.nml` between each run in this task. This ensures that all calculations start from the same, converged, potential in `Fe.nml`.


## 2. Si

### A) Convergence to self-consistent solution

Copy the starting fixture into a working directories and move to that directory.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs/alternative-si
cp -r tests/scf/cases/bulk/diamondSi runs/alternative-si/scf
cd runs/alternative-si/scf
```

Open **`runs/alternative-si/scf/input.nml`** and investigate the inputs. 
What looks different from bcc Fe?

This system uses the Chebyshev method instead of block recursion since we expect Si to be non-metallic. The Chebyshev method needs a larger number of recursion levels `lld`. 

Converge the system as you did for bcc Fe. How many iterations were needed?


### B) Effect of number of Chebyshev polynomials

Follow the procedure you used for changing the number of recursion levels `lld` in the bcc Fe example: First, generate the original DOS with `plot_dos.py`, then make a new directory and copy your converged output there.

```bash
mkdir ../dos-comparison
cp input.nml lattice.nml ../dos-comparison
cp Si1_out.nml ../dos-comparison/Si1.nml
cp Si2_out.nml ../dos-comparison/Si2.nml
cd ../dos-comparison
```

Investigate the effect `lld` has on the DOS for Si.

### C) Comparison with k-space DOS

Copy the starting Si files to a new directory:

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs/kspace-si
cp -r tests/scf/cases/bulk/diamondSi runs/kspace-si/scf
cd runs/kspace-si/scf
```

and open `input.nml` in the new directory. In the `%self` namelist, add the line `use_kspace = .true.` and change `nstep` to `100`. Then run the calculation and generate the DOS plots. 

How does the DOS differ from the Chebyshev DOS? 


## 3. Fe exchange interactions

Now we can go back to bcc Fe and calculate exchange interactions between neighbouring spins. For this we need a new directory to run things in. Copy the **converged** bcc Fe files to the new directory.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir runs/alternative-fe-jij
cp runs/alternative-fe-scf/Fe_out.nml runs/alternative-fe-jij/
cp runs/alternative-fe-scf/input.nml runs/alternative-fe-jij/
cp runs/alternative-fe-scf/lattice.nml runs/alternative-fe-jij/
cp runs/alternative-fe-scf/report.out runs/alternative-fe-jij/
cp runs/alternative-fe-scf/fort.805 runs/alternative-fe-jij/clust.xyz
cd runs/alternative-fe-jij/
```

In addition to the input files, we also copied `fort.805` to `clust.xyz` for reasons that soon will be explained. We will also need the `report.out` file for later tasks.

In the `runs/alternative-fe-jij` directory, open `input.nml` and change the `&calculation` namelist to

```
&calculation
post_processing = 'exchange'
verbose = T
/
```

This sets the task for the code to calculating exchange interactions. Next we need to decide which pairs of atoms to calculate `Jij` for. We do this by adding the following lines in the `&lattice` namelist.

```
nij = 2
ijpair(1, :) = 1, 2634
ijpair(2, :) = 1, 2635
```

I.e. two interactions between the atoms that are given on each `ijpair` line. Note that the atoms are given by their index in the cluster. This means that 
a) it is advised to have the same cluster size settings in the SCF run as in the exchange run, and b) you can find the indices for the atoms in the `clust.xyz` file which shows the coordinates for each atom in the cluster.

Save `input.nml` and run the code. 

Check the results in the `jijso.out` file. Here, columns 1–2 are symbolic-type IDs,
3–5 are displacement/a, 6 is J in mRy, and 7 is distance/a. 

### Plot the selected exchange rows

From the run directory, plot the actual exchange output:

```bash
plot_exchange.py jijso.out
```

The plot shows J in meV versus distance/a, grouped by **directed symbolic-type
pair**. CSV retains the seven source columns, including J in **mRy** and
Cartesian displacement/a. 

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

The `rsutil` package can also export `rslmto` calculations to minimal `UppASD`input. For this we need the input files `input.nml` and `lattice.nml` as well as the exchange file `jijso.out` and the log file `report,out`. The latter is needed since `UppASD` needs moment magnitude information in addition to structure and exchange inputs. So this can be done as follows:

```
cd ..
pwd
rsutil --export --jij Atom_1/jijso.out
```

This creates a minimal `UppASD` setup in the folder `uppasd_out`, unless another name is specified with the `--out OUTDIR` flag.

Inspect the resulting `UppASD` input. Does it look correct`

## half-Heusler NiMnSi


### Selectiv exchange pair selection

While the `rsutil` exchange setup tool is convenient, one sometimes want more control over the selection. As an example, the half-Heusler compound NiMnSi only contains one magnetic ion, Mn. Ni does have a small moment here, but can be approximated to not be relevant here. 

So here we want to setup an exchange calculation between Mn-Mn pairs only and this we do manually. Start by copying the converged NiMnSi system to a run folder.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs/exchange-nimnsi
cp  ./tutorials/MnNiSb/* runs/exchange-nimnsi
cd runs/exchange-nimnsi
```
Now, open `clust.xyz` in the new folder and identify relevant Mn-Mn couplings. Enter these into the existing `input.nml` as `jijpair` lines. And remember to set `njij` to the number of exchange pairs you have defined. Also comment out the `pre_processing` line and add the relevant `post_processing` line.

Run the calculation and plot the resulting exchange couplings using the `plot_exchange.py` script. Does it seem as you have captured Mn-Mn couplings?

### Advanced: DOS features across methodologies

Among its interesting properties, NiMnSi is half-metallic, i.e. it has a gap in one of the spin channels. This is difficult to resolve with block recursion and the Chebyshev method does need a lot of recursion levels. 

To get a feel of how the DOS is reconstructed from these different approachs, one can do single-shot DOS calculations from the converged potentials. To do this, edit `input.nml` accordingly, i.e. change the values of `use_kspace`, `lld`, and `recur`. To ensure a single iteration, also change `nstep` to `1`.

Perform the `input.mml` changes and run the code. After each run, genreate the DOS with `plot_dos.py` and store the relevant images so you don't overwrite them. What differnces do you see across the different methods?

## Heusler Mn2FeSi

Another well studied material of the Heusler family is `Mn2FeSi`. Here we can try to use our `rsutil` functionalities to the max with minimal human interaction, i.e. moving towards a high-troughput methodology.

To prepare the tasks, copy the reference data to a run directory and list the copied files.

```bash
cd "$(git rev-parse --show-toplevel)"
mkdir -p runs/ht-mn2fesi
cp  ./tutorials/Mn2FeSi/* runs/ht-mn2fesi
cd runs/ht-mn2fesi
ls .
```

### From .cif to RS-LMTO-ASA input

As you can see, the run directory only contains a `.cif` file, i.e a crystal structure descriptor. In order to run `rslmto.x` we thus need to process the `.cif` file. This can be done with `rsutil` using the following syntax:

```
rsutil setup bulk --cif EntryWithCollCode221028.cif --empty_spheres
```

Open the generated `input.nml` and `lattice.nml` and inspect them. If everything looks good, run the code to converge the system. *Hint: If slow, you can run a `use_kspace` enabled run to speed things up.*

### Automatic exchange pair setup

Once converged, you can use `rsutil` tool again to generate the relevant exchange couplings for the system. Note: Before the setup, ensure that `use_kspace = .false` is set in the `input.nml` file.

```
rsutil --exchange --cutoff 1.0 --unique
```

This time, the script will generate four folders, one for each atom in the unit cell. Here we use a limited cutoff 

### Running the exchange calculations

Now you can run the `rsutil.x` code in each directory to calculate the exchange interactions, and then paste the results to a single output file.

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

Unless you specify an output directory name, the generated files are placed in an `uppasd_out` subdirectory. Examine the outputs, do the look as expected?

### End of tutorial

References: 
The pinned RS-LMTO source/fixtures; RS-utils `0.11.0rc1` supplied
source, `docs/export_submodule.rst`, `docs/notes/UppASD_format.md`, and
`rsutil/export/uppasd.py`; [RS-utils manual](https://aibergman.github.io/RS-utils/).
