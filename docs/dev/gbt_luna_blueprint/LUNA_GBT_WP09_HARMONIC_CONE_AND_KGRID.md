# Luna Prompt — GBT-WP09: Harmonic Cone-Angle and k-Grid Convergence

## Mission

Establish a numerically controlled small-cone regime for frozen spin spirals and separate physical cone-angle dependence from finite-k-grid pure-gauge quadrature error.

This work package must preserve raw energies. Do not report only a processed dispersion.

## Dependency

WP08 modes must be available and energy semantics settled.

## Core quantities

For each q and cone angle theta record at minimum:

\[
E(q,\theta),\qquad E(q,0),\qquad E(q_{ref},\theta),\qquad E(q_{ref},0).
\]

For k-space MFT, the same-q gauge-subtracted quantity

\[
\Delta E_{gauge}(q,\theta)=E(q,\theta)-E(q,0)
\]

is useful because at finite k mesh the pure-gauge `theta=0` energy can acquire q-dependent quadrature error even when the exact band problem is equivalent to shifted k points.

However, also retain

\[
\Delta E_{pure}(q)=E(q,0)-E(0,0)
\]

and verify that it tends to zero with k-mesh refinement.

## Required cone-angle sweep

For representative small q values use a deterministic set such as:

- 2 deg;
- 5 deg;
- 10 deg;
- 15 deg;
- 20 deg as a nonlinearity monitor.

Do not assume all are harmonic.

For each angle evaluate

\[
K(q,\theta)=\frac{\Delta E(q,\theta)}{\sin^2\theta}.
\]

Identify the plateau region from numerical evidence.

## Required k-grid convergence

Use several systematically refined meshes. At each mesh report:

- pure-gauge drift `E(q,0)-E(0,0)`;
- gauge-subtracted spiral energy;
- raw reference-subtracted spiral energy;
- Fermi energy/electron-count consistency;
- harmonic ratio `K(q,theta)`.

Convergence target must be set relative to the small spiral energy being fitted, not a generic total-energy threshold.

## q-grid policy

Check whether q/2 shifts map the chosen k mesh onto itself. When they do, the pure-gauge drift should be exceptionally small. Use at least one such commensurate mesh as a control if practical.

## Modes

At minimum characterize bare MFT and corrected constrained MFT. Add SCF where computationally reasonable after the former are numerically stable.

## Tests/automation

Create a small analysis utility/script that reads machine-readable sweep files and reports:

- plateau spread;
- k-grid convergence of each energy definition;
- angle range admitted to harmonic fitting.

Do not bake system-specific expected energies into generic unit tests. Use functional thresholds on convergence/invariance.

## Deliverables

- convergence fixtures/scripts;
- `docs/dev/GBT_REAUDIT_WP09_HARMONIC_KGRID.md`;
- recommended default diagnostic output for future GBT runs.

## Completion checklist

- [ ] Raw `E(q,theta)` retained.
- [ ] Same-q `E(q,0)` retained.
- [ ] Pure-gauge drift plotted/tabulated versus k mesh.
- [ ] Several cone angles tested.
- [ ] `DeltaE/sin^2(theta)` plateau identified.
- [ ] Nonlinear angles excluded rather than forced into fit.
- [ ] k-grid convergence is tighter than target spiral energy uncertainty.
- [ ] At least one mesh mapping control considered.
- [ ] Bare and corrected MFT compared.
- [ ] Evidence report gives admitted angle/k-grid regime.

## Commit

`test: establish GBT harmonic and k-grid convergence`
