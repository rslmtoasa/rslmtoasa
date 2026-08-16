# VAL-11 — abspinlib spin-integrator validation

## Scope

This validation isolates the standalone `abspinlib` spin solver from the
self-consistent LMTO field generator. The tested API is the production
Depondt-Mertens predictor/corrector in
`source/include_codes/abspinlib/depondt.f90`, exposed as
`depondt_evolve_first` and `depondt_evolve_second`. The `abSpinlib` wrapper's
`asd_pred`/`asd_corr` procedures delegate to the same implementation and are
not a second numerical integrator.

The abspinlib CMake target contains no separate Heun implementation. The
application module also contains `spin_dynamics::asd_pred_euler` and
`asd_corr_euler`, but those are application-level helpers rather than
abspinlib integrators and are not promoted by this report. The existing
`Example_bulk_bccFe_sd_smoke` remains an LMTO orchestration smoke test.

## Oracles and test problems

At zero temperature the test supplies deterministic effective fields directly
to the Depondt API. No Python numerical integrator or duplicate production
solver is used as an oracle.

For a constant field (\mathbf B=B\hat z), damping \(\alpha), and unit spin,
the low-level API follows

\[
 \dot{\mathbf m}=\frac{\gamma}{1+\alpha^2}
 \left[\mathbf B\times\mathbf m + \alpha\left(\mathbf B-
 \mathbf m(\mathbf m\cdot\mathbf B)\right)\right].
\]

The analytic solution used by the test is

\[
 u(t)=\tanh\left[\operatorname{atanh}(u_0)+
 \frac{\gamma\alpha B}{1+\alpha^2}t\right],
 \qquad
 \phi(t)=\frac{\gamma B}{1+\alpha^2}t,
\]

with the transverse component rotated by (+\phi) around \(\hat z). The
zero-damping limits (\lvert\mathbf m\rvert) and
\(\mathbf m\cdot\mathbf B) are exact invariants, and the phase sign checks
the precession direction.

The exchange case uses two equal-magnitude spins with prescribed isotropic
fields (\mathbf B_1=J\mathbf m_2) and
\(\mathbf B_2=J\mathbf m_1\), refreshed from the predictor configuration.
For (E=-J\mathbf m_1\cdot\mathbf m_2), the zero-damping textbook
invariants are the two spin norms, energy, and total spin.

## Measured evidence

The run was the serial GNU Fortran Debug build, using
`ctest --test-dir build-rf-debug -R '^UnitAbspinlibIntegrators$'`.

| Check | Result |
|---|---:|
| Constant-field precession norm error | (3.89\times10^{-15}) |
| Constant-field (\mathbf m\cdot\mathbf B) error | (0) |
| Constant-field analytic final-vector error | (4.80\times10^{-15}) |
| Damped-field minimum projection increment | (7.71\times10^{-8}) (positive) |
| Damped-field analytic final-vector error | (1.59\times10^{-7}) |
| Damped-field norm | (1.0000) to printed precision |
| Timestep errors, (2,1,0.5\times10^{-11}) s | (6.37\times10^{-7}, 1.59\times10^{-7}, 3.98\times10^{-8}) |
| Observed timestep order | (2.0001) |
| Two-spin norm error | (4.06\times10^{-14}) |
| Two-spin energy error | (2.46\times10^{-18}) |
| Two-spin total-spin error | (5.72\times10^{-14}) |

## Maturity and limitations

**Depondt-Mertens deterministic CPU integrator: Validated (scoped).** The
claim covers the direct effective-field API at zero temperature for constant
fields and a two-spin isotropic exchange field, including norm preservation,
precession direction/frequency, damping relaxation, analytic agreement, and
second-order timestep convergence.

This does not validate stochastic thermal noise, spin-transfer or spin-Hall
torques, induced-moment handling, time-dependent fields, or the self-consistent
LMTO field feedback. No UppASD comparison was used because the analytic
constant-field solution and exact exchange invariants provide the independent
oracles for this scope. Full ab-initio spin dynamics therefore remains
**Experimental**.

## Checklist

- [x] Available abspinlib integrators inventoried.
- [x] Constant-field precession validated.
- [x] Spin norm validated.
- [x] Energy invariant checked in the zero-damping exchange case.
- [x] Damping relaxation checked.
- [x] Timestep convergence measured.
- [x] Optional external UppASD comparison assessed as unnecessary for this analytic scope.
- [x] No electronic-structure code changed.
- [x] Solver maturity documented.

## Changed files and test command

- `tests/validation/test_abspinlib_integrators.f90`
- `CMakeLists.txt`
- `docs/dev/PHASE_I_STABILIZATION.md`
- `docs/dev/VAL-11_ABSPINLIB_INTEGRATORS.md`
- `docs/dev/PHASE_II_VALIDATION.md`

The focused command is:

```sh
ctest --test-dir build-rf-debug -R '^UnitAbspinlibIntegrators$' --output-on-failure
```

Commit message: `Validate abspinlib spin integrators`
