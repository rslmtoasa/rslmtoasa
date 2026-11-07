# Run RS-LMTO-ASA in Codespaces

[Open the tutorial in Codespaces](https://codespaces.new/rslmtoasa/rslmtoasa/tree/tutorial/codespaces)

Sign in with your own GitHub account. Check that the branch is
`tutorial/codespaces`, choose the 2-core machine, and create your codespace.
Each student gets an independent environment and uses their own available
personal Codespaces allowance. Check the billing owner before creating it.
Initial setup downloads the environment and compiles the code; allow several
minutes. No local compiler, Docker installation, or Python interface is needed.

## Edit inputs and run

Setup copies the supplied bcc Fe example into `runs/bccFe`, leaving the original
example intact. Click `runs/bccFe/input.nml` in the file explorer, edit it, and
save with Ctrl+S (Cmd+S on macOS).

Open **Terminal > New Terminal** and run from the repository root:

```bash
cd runs/bccFe
rslmto.x input.nml > run.log 2>&1
```

This uses two OpenMP threads and one OpenBLAS thread. MPI is disabled.
Calculations run only when you start them; opening the codespace does not launch
a calculation. The supplied input allows 100 SCF steps. For a quick installation
check, change `nstep = 100` to `nstep = 1`; this is not a converged calculation.
Runtime and memory depend on the example and input settings.

Open `run.log`, `report.out`, and other output files in the editor. Check the
convergence information before interpreting scientific results. To monitor a
running calculation, use `tail -f run.log` in a second terminal in the same
directory (Ctrl+C stops the monitor).

## Plot inside Codespaces

Python 3, NumPy, Matplotlib, and gnuplot are installed. Save plots as image files,
then click them in the file explorer to preview or download them. A graphical
desktop or a Python interface to the Fortran code is not needed.

After running the bcc Fe example, use its first two DOS columns as a plotting
example (from `runs/bccFe`):

```bash
python3 - <<'PYTHON'
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("totaldos.out")
plt.plot(data[:, 0], data[:, 1])
plt.xlabel("Energy (Ry)")
plt.ylabel("DOS (first data column)")
plt.tight_layout()
plt.savefig("dos-python.png", dpi=150)
PYTHON
```

Or use gnuplot in the same directory:

```bash
gnuplot <<'GNUPLOT'
set terminal pngcairo size 1000,600
set output 'dos-gnuplot.png'
set xlabel 'Energy (Ry)'
set ylabel 'DOS (first data column)'
plot 'totaldos.out' using 1:2 with lines title 'DOS'
GNUPLOT
```

For additional pip packages, create a virtual environment from the repository
root. This keeps Ubuntu's managed Python installation intact and makes the
preinstalled NumPy and Matplotlib available:

```bash
python3 -m venv --system-site-packages .venv
source .venv/bin/activate
python -m pip install scipy  # example additional package
```

To include more packages for everyone, add apt packages to
`.devcontainer/Dockerfile`, or add a virtual-environment setup and
`python -m pip install -r ...` step to `.devcontainer/setup.sh`.
After changing the Dockerfile, use **Codespaces: Rebuild Container** from the
command palette. Existing codespaces need to pull the updated branch first;
new codespaces use the updated configuration automatically.

## Download results for local plotting

Right-click a file in the browser file explorer and select **Download**.
For the whole calculation, run from the repository root:

```bash
tar -czf bccFe-results.tar.gz -C runs bccFe
```

Download `bccFe-results.tar.gz`, extract it locally, and plot with your preferred
tools. Files under `runs/` and result archives are ignored by Git.

## Stop, resume, and rebuild

Stop the codespace at https://github.com/codespaces when finished to stop compute
usage. Stopped codespaces still consume storage allowance. Download important
files before deleting a codespace or its automatic retention period expires.

To rebuild after editing source files, run from the repository root:

```bash
cmake --build build-codespaces --parallel 2
```

If setup fails, inspect the creation log, then retry from the repository root:

```bash
bash .devcontainer/setup.sh
```

Setup preserves an existing `runs/bccFe` directory. To start another calculation,
copy an example into a new directory under `runs/` and run from that directory.

The same configuration can be opened locally with VS Code's Dev Containers
extension and a working Docker installation. Local use consumes local resources.
