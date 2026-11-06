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

The original `ham_vec_matmul` local/bulk loops now serve both conductivity and
`orbital_modern`. `ham_hoh_vec_matmul` retains the original two explicit passes,
each with separate impurity and bulk loops. Its second pass forms h_ij*o_j
using the actual intermediate site before multiplying by h*psi. The routines apply the original
shift and scale internally. Intermediate bare h products and torque products
use a=1, b=0; velocities remain unscaled. The non-HOH action is h+SOC.

All products honor `izero` as input support and return conservative output
support in `idum`. HOH expands support for each hop, preserving two-hop paths.
The velocity first pass builds h*psi and v*psi together, as in the original
routine. They share a conservative one-hop map. The second pass expands that
map for both v*o*h and h*o*v, including the on-site h*o contribution. Local accumulation preserves site support.
The input map is restored on return, including for the output observable, which
must not change the support of the evolving Chebyshev state. Masks include the
original support so the recurrence subtraction and energy shift remain valid.
Zero neighbor entries are checked before indexing the map.

The routines remain in `recursion.f90`, using explicit BLAS `zgemm` calls and
OpenMP loops over output sites. Each thread writes its own output site and map
entry. Local impurity rows use `hall`; bulk rows use `ee` with the type mapping.
The caller allocates six scratch planes once: 1:3 for H/V, 4:6 for observable
products. Within the HOH loops these have the familiar names psi1, psi2 and
hohpsi, assigned to existing storage once at entry. Older HOH Hamiltonian callers
can still use the existing object scratch arrays.
Input, output and scratch must not overlap. No new Fortran file/module or input
setting is introduced. The standalone transport module and obsolete operator
builders remain removed; selectors still come from `hamiltonian.f90`.

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

## Which multiplication routines are used?

The conductivity and conductivity_p2rs workflows call `compute_moments_stochastic`.
Its H dispatcher selects the original multiplication routines and passes a,b;
it contains no independent Hamiltonian implementation. Observable selection
remains in `apply_kubo_operator` / `kubo_operator_vec_matmul`.

| Routine | Current use |
| --- | --- |
| `setup_kubo_operators`, `set_kubo_operator_slot` | Configure slots and selectors |
| `apply_kubo_operator`, `kubo_operator_vec_matmul` | Apply current, accumulation or torque |
| `kubo_ham_vec_matmul` | Select original H routine; unscaled defaults for torques |
| `ham_vec_matmul` | Original masked local/bulk H loops, including shift/scale; bare h when SOC disabled |
| `ham_hoh_vec_matmul` | Two bare h products, overlap, E_nu and optional SOC; shift/scale |
| `velo_vec_matmul` | Masked bare velocity, including local impurity rows |
| `velo_hoh_vec_matmul` | Full v-v*o*h-h*o*v with intermediate support propagation |
| `onsite_vec_matmul` | Masked local SOC products for SOC torque |
| Former `hopping_vec_matmul` | Removed; bare H uses the original H loops |
| Former operator block builders and standalone transport module | Removed |

`chebyshev_orbital_mod` still calls the H routines directly with its original
arguments. Other SCF/Lanczos and Chebyshev routines are unchanged.

## Direct comparison with the old multiplication routines

The initial checks compared the new operators with independent dense matrices.
They did not execute the old kernels. The regression driver now also calls the
test-only archived copies of the pre-refactor H kernels and the old velocity
kernels from commit `e604f4029e59d35b6e53551dd2d95919a9baf8e0`. Active statements
are unchanged apart from procedure names; archived kernels are not linked into
the solver. This keeps the old/new comparison independent of the modified H routines.

The comparison accounts for these different contracts:

- Both H paths are compared with the same internal `(H-b)*psi/a` shift/scale.
- Old routines require initialized scratch and valid support masks. Tests cover
  full support and a state localized on the first site, repeated calls,
  `nmax=0` and `nmax=2`, zero/nonzero overlap, and one/two OpenMP threads.
- For the old routines, neighbor lists are compact: their compound `.and.`
  guard does not guarantee short-circuit evaluation before indexing a neighbor.
  The new routines remain separately tested with zero neighbor entries.
- Bulk non-HOH velocity agrees with the old ordinary (`c_or_n='n'`) action.
- HOH velocity deliberately differs: `V_new*psi - V_old*psi = -h*o*v*psi`.
  The old routines do not implement impurity velocity rows; these are checked
  against the dense reference instead of being compared with omitted rows.

The deterministic three-site fixture checks against independent dense H/V and
observable matrices, including the HOH-generated 1--3 bond absent from bare h.
Additional tests compare masked and all-site actions for every observable kind,
both slots, both Hamiltonians and impurity settings, and assert the output maps.
Six successive shifted/scaled Chebyshev steps check support preservation and
scratch reuse. Empty support must produce zero output and an empty map.

Each run writes `legacy_comparison.csv` in its pytest temporary directory.
These are multiplication-kernel comparisons; they do not validate the complete
`orbital_modern` algorithm or the unresolved conductivity integration/units.

The routine headers document scaling, scratch, maps and callers. Products use
explicit `zgemm` calls. No `ASSOCIATE` blocks are used in the transport routines;
the HOH scratch pointers only select existing caller or object storage. The
explicit H/V loops contain no calls to other multiplication helpers. The
observable dispatcher and six-plane storage interface remain in place; the
routines are therefore not literal copies of the historical source.

Validation of the restored masked implementation: GNU Fortran 13 Debug with
bounds checks, reference BLAS/LAPACK, OpenMP 1 and 2 threads: 10 tests and
15 subtests passed. Both caller-owned and original object HOH scratch paths
are covered. The complete `rslmto.x` executable also builds successfully.
No material-scale performance benchmark was performed.
