# Run RS-LMTO-ASA in Codespaces

For the guided session, see the [two-hour bcc Fe SCF and exchange tutorial](../docs/tutorial-bcc-fe-2h.md).

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

## Bundled RS-utils (release-prep)

RS-utils installs automatically during Codespaces setup. Students do not
upload a package or need access to the private source repository.
`vendor/rsutils/` contains the author-provided x86_64 and ARM64 wheel archives;
the installer selects the wheel matching the current Python/CPU/platform.
The distribution is `rslmto-utils` 0.11.0rc1; command/import remain `rsutil`.

After setup, open a new terminal and run:

```bash
rsutil --help
```

Installation uses an isolated `.venv-rsutils` with the `bands` extras
(including `seekpath`). It does not inherit Ubuntu Python packages, avoiding
mixed Matplotlib/`mpl_toolkits` versions. The installer checks archive checksums,
package version, CLI, and the loaded native Fortran library. It replaces any
older installation even if the version number is unchanged. Python 3.11+
is required; the Ubuntu 24.04 image uses Python 3.12. Python dependencies are
resolved from PyPI, so initial setup still needs internet access.

To reinstall or repair an existing mixed Matplotlib environment, run from
the repository root (a container rebuild is not needed):

```bash
bash .devcontainer/install-rsutils.sh
source .venv-rsutils/bin/activate
```

For an instructor override, pass an individual compatible Linux `.whl`, one
of the pinned wheel bundles, or the supplied `RS-utils-release-prep.zip`:

```bash
bash .devcontainer/install-rsutils.sh "work/rsutils/path-to-instructor-package.whl"
```

Use the actual provided filename. The source ZIP requires compilation;
bundled wheels do not. A different source/bundle snapshot requires updating
the pinned checksums (and expected version if changed) together. Local uploads
and venvs are ignored by Git; bundled wheels are deliberately versioned.

## Recover an existing Codespace

In a working terminal, run `git pull --ff-only` on the tutorial branch.
Open the command palette (F1) and choose **Codespaces: Rebuild Container**.
If normal rebuild continues using an old image, use
**Codespaces: Full Rebuild Container**. Results under `runs/` are preserved;
download important outputs first.

The alternative branch previously requested `python3-f90nml` from apt, but
Ubuntu 24.04 has no package by that name. That image build fails before setup
runs. The corrected setup installs `f90nml==1.5.0` with pip in
`.venv-tutorial`, keeping Ubuntu's system Python intact.

If recovery mode persists, open **Codespaces: View Creation Log** and locate
the **first failing Docker/build/post-create command**, including the error
lines above it. A later log that only starts the recovery container can end
successfully without showing the original configuration failure.

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
