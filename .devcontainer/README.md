# Run RS-LMTO-ASA main_rc in Codespaces

[Open the tutorial in Codespaces](https://codespaces.new/rslmtoasa/rslmtoasa/tree/tutorial/codespaces-main-rc)

Sign in with your own GitHub account. Check that the branch is
`tutorial/codespaces-main-rc`, choose the 2-core machine, and create your codespace.
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

This branch is based on `main_rc`. The build includes spglib for symmetry
operations and uses two OpenMP threads and one OpenBLAS thread. MPI is disabled.
Calculations run only when you start them; opening the codespace does not launch
a calculation. The supplied input allows 100 SCF steps. For a quick installation
check, change `nstep = 100` to `nstep = 1`; this is not a converged calculation.
Runtime and memory depend on the example and input settings.

Open `run.log`, `report.out`, and other output files in the editor. Check the
convergence information before interpreting scientific results. To monitor a
running calculation, use `tail -f run.log` in a second terminal in the same
directory (Ctrl+C stops the monitor).

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
