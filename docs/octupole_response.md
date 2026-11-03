# Orbital octupole response

These are local orbital-angular-momentum multipoles, not collective magnetic
multipoles of a multi-site spin configuration. L matrices use hbar=1.

## Components and normalization

For a != b, the new pure rank-three component is

    Oabb = {La,Lb^2}/2 - {La,L^2}/10 - La/10
         = {La,Lb^2}/2 - (L^2 + I/2) La/5.

Oxyz is the six-permutation symmetrization of Lx Ly Lz divided by six.
The physical PDF convention (cubic angular momentum divided by hbar^2)
is obtained by multiplying these dimensionless matrices by hbar.

| Preferred pol | Harmonic-label alias | Old input (now rejected) | Explicit mixed observable |
|---|---|---|---|
| Oxyy | Ox5y2r2 | Ox3y2r2 | Mx3y2r2 |
| Oyzz | Oy5z2r2 | Oy3z2r2 | My3z2r2 |
| Ozxx | Oz5x2r2 | Oz3x2r2 | Mz3x2r2 |
| Oxzz | Ox5z2r2 | Ox3z2r2 | Mx3z2r2 |
| Oyxx | Oy5x2r2 | Oy3x2r2 | My3x2r2 |
| Ozyy | Oz5y2r2 | Oz3y2r2 | Mz3y2r2 |
| Oxyz | -- | Oxyz remains valid | -- |

Lowercase aliases are accepted. Bare xyy/yzz/zxx/xzz/yxx/zyy/xyz also work.
Harmonic-label aliases have exactly the same normalization as Oabb: they do
NOT multiply the operator by five. These seven components are independent
but need not be mutually orthogonal or equally normalized.

The explicit M labels reproduce the old formula {La,Lb^2}/2 - {La,L^2}/6.
They mix dipole and octupole and must not be reported as pure octupoles.
Old Oa3b2r2 inputs fail with a migration message rather than silently changing
the physical observable or silently retaining a mixed moment.

Pure components vanish in s and p, have zero dipole projection in every shell,
and obey sum_a [La,[La,O]]=12 O. The distinction applies to both accumulation
and conventional current. Selecting an M label deliberately requests a mixed
observable, including its dipole-current contribution.

## Input slots

| Setting | Role |
|---|---|
| linear_out | Measured response operator, slot a |
| pol_alpha | Multipole component of the output operator |
| v_alpha | Output current direction, ignored by a local accumulation |
| linear_in | Driving operator, slot b |
| pol_beta | Multipole component of the input operator |
| v_beta | Input current direction; electric-field direction for charge input |

Both slots accept `octupole` and its synonym `octupole_current` for J={O,v}/2,
and `octupole_accumulation` for the local O. Operator keywords are lowercase.
Use `charge` input for the electric-field-driven response. Other input
operators require defining the conjugate external perturbation and units.
Directions are normalized; their magnitudes are not field strengths. Both
vectors must be finite and nonzero; setup validates and normalizes both directions.
pol_alpha/pol_beta have length 10; every documented component fits. Defaults
are 'z', which is not an octupole component: set the relevant pol explicitly.

## Hall-current example

Merge `example/transport/octupole_hall.nml` into the matching sections of a
complete converged material input. It is a response-settings fragment, not
a standalone material fixture. It requests J_y^(Oxyy) driven by E_x,
i.e. the transverse coefficient sigma_yx^(Oxyy). A nonzero value depends on
the material symmetry; this input does not guarantee a nonzero Hall effect.
Do not create duplicate &control or &hamiltonian sections.

For accumulation, replace linear_out='octupole' with
linear_out='octupole_accumulation'. The accumulation example provides this.
For an input octupole current instead, use linear_in='octupole', pol_beta,
and v_beta; set linear_out independently with pol_alpha and v_alpha as needed.

## Current scope and limitations

All supported `linear_in`/`linear_out` operators now work with either value of
`hoh`, including local impurity rows (`nmax > 0`). The stochastic Kubo kernel
uses one matrix-free Hamiltonian action throughout:

$$H=h-hoh+E_\nu+H_{\rm SOC},\qquad V=v-voh-hov.$$

Here `o` is the on-site overlap correction and `v = -i[X,h]` in the code's
units. For `hoh=.false.`, the action reduces to `H=h+H_SOC` and `V=v`.
Both product-rule terms are retained without assuming commuting matrices.

| Input keyword | Operator action |
| --- | --- |
| `charge` | Full velocity `V` |
| `spin`, `orbital`, `quadrupole` / `oap`, `octupole` / `octupole_current` | `(A V + V A)/2` for the selected S/L/Q/O |
| Corresponding `*_accumulation` (including `oap_accumulation`) | Local selected S/L/Q/O |
| `spin_torque` | `-i[S,H-H_SOC]`, including spin-dependent HOH and E_nu terms |
| `spin_soc_torque` / `soc_spin_torque` | `-i[S,H_SOC]` |
| `orbital_torque` | `-i[L,H]` |

Polarization selectors and direction conventions are unchanged. The input slot
retains the symmetric bond weight `max(velocity_scale(i),velocity_scale(j))`
on each bare velocity insertion in both HOH terms. The output slot is unweighted.
For ordinary Hall response use `linear_in='charge'` and the desired output
current, with transverse `v_beta` and `v_alpha`; set `hoh=.true.` in `&hamiltonian`.

All intermediate actions traverse the complete cluster so two-hop paths cannot
be truncated by the older recursion support masks. This costs extra sparse
matrix actions and temporary vectors, but does not assemble a dense Hamiltonian.
The multiplication routines live in `recursion.f90`, alongside the existing
`ham_vec_matmul`/`ham_hoh_vec_matmul` routines. They access the existing Hamiltonian
and lattice objects directly and use explicit `zgemm` calls. The standalone
transport module and unused real-space operator block builders have been removed,
including their outdated HOH rejection. The selected S/L/Q/O matrices still come
from `hamiltonian.f90`.

The Kubo caller allocates a six-plane scratch array once per moment calculation.
Planes 1:3 serve H/V products; planes 4:6 hold the operator intermediates.
Input, output, and scratch must not overlap. The multiplication routines perform
no explicit allocations and retain every intermediate site. This cleanup does
not restore the old one-sided HOH formula or the old support masks. It changes
no `input.nml` settings and does not change non-transport recursion routines.

This change adds and tests the operators, not a validated material conductivity.
Previously identified energy-integration bounds, normalization/unit conversion,
and local-projection issues in the common conductivity code remain outstanding.
Quantitative Hall results must wait for those repairs and end-to-end validation.
The current is the conventional anticommutator; no conserved-current or
torque-dipole correction is implied.

Malformed &control and &hamiltonian namelists now fail immediately rather than
continuing with partially read settings. Obsolete cond_type/js_alpha inputs
must be migrated to linear_in/linear_out and pol_alpha/pol_beta.

## Validation

Run from the repository root:

    python tests/regression/test_quadrupole_operators.py
    python tests/regression/test_octupole_operators.py

The octupole checks evaluate the actual selector expressions and actual hcpx
basis transformation extracted from source. They test all aliases, pure rank,
p-shell vanishing, dipole orthogonality, legacy migration, seven-component
independence. Source routing and example namelist checks cover the input
contract. These Python-only checks are complemented by the compiled tests below.

The compiled HOH test requires gfortran, CMake, BLAS/LAPACK, NumPy, f90nml, and pytest:

    FC=gfortran OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/regression/test_hoh_transport.py

It builds and links the production recursion module in Debug mode with bounds
checking, then compares its Fortran actions against dense matrices for all
operator families, all seven pure octupoles, HOH on/off, zero/nonzero overlap,
impurity rows, input weighting, and one/two OpenMP threads. It checks complete
matrix elements and Hermiticity, including an induced two-hop bond. Both response
slots and all supported operator keywords are also exercised through compiled
setup and selector routines, with independently normalized directions. This is
operator-level validation, not a material-level Kubo integration benchmark.

`FC` and `CMAKE` select the compiler and CMake executable. Additional CMake
configuration arguments (for example explicit BLAS/LAPACK library paths) can be
provided through `RS_LMTO_TEST_CMAKE_ARGS`. The fixture runner is an explicit
`transport_operator_test` target and is not built by default.
