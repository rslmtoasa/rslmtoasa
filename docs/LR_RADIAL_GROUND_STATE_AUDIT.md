# LR-01 radial magnetic ground-state audit

## Verdict

**PASS** for the supported scalar-relativistic, collinear, two-channel radial
ground-state contract. The exact converged radial density and multiplicative
spin-dependent XC arrays now survive the `ATOMSC` call in
`symbolic_atom%radial_ground_state`; a machine-readable copy is emitted by
`report`. The live legacy-XC and libXC oracles, independent charge and moment
quadratures, channel/decomposition checks, and the focused unit supplement all
pass.

This is a ground-state data contract only. It does not implement TDDFT,
`chi0`, response reconstruction, Ward/Goldstone repair, a terminator, or a
generalized-overlap response route.

## Audit identity and initial evidence

The audit started on branch `fable_v4` at the exact HEAD
`14f67ca8b92416b60b3704f2cab67f9d71c67079`:

```text
14f67ca Validate LMTO radial augmentation against SCF density
497ff42 Certify one-electron Green-function contracts for linear response
2ba708f td-dft: establish clean-room linear-response baseline
```

The initial worktree was clean. Before changing the radial path, the existing
LR-GF-01 and LR-BASIS-00 evidence was present and the focused baseline tests
were run. The existing evidence files are
[`LR-GF-01_GF_CONTRACT_EVIDENCE.md`](LR-GF-01_GF_CONTRACT_EVIDENCE.md) and
[`LR-BASIS-00_RADIAL_AUGMENTATION_EVIDENCE.md`](LR-BASIS-00_RADIAL_AUGMENTATION_EVIDENCE.md).
The existing `UnitLrBasisRadial`, `UnitLrBasisAngular`, and
`UnitLrBasisAugmentation` tests remained passing after the change.

The validation build used GNU Fortran 13.3, OpenMP, MPI off, and libXC 5.2.3.
The no-libXC configuration also passed the negative selector guard;
`Val22LrRadialGroundStateLibxc` was then run with libXC enabled.

## Scope and supported state

The certified state is:

| item | contract | status |
| --- | --- | --- |
| radial equation | existing scalar-relativistic `RSEQSR`/`PHDFSR` machinery | PASS |
| radial channels | local channel 1 = up, local channel 2 = down | PASS |
| global modes | collinear two-channel state, including the local radial path when `control%nsp` is not 2 | PASS |
| XC | legacy RS-LMTO and supported libXC LDA/GGA routes | PASS |
| SOC, noncollinear radial spinors, Hubbard response, generalized overlap | not silently inferred by this contract | DEFERRED / out of scope |

## Live call graph and capture point

The audited SCF path is:

```text
self%run
  -> run_scf
     -> lmtst
        -> atomsc
           -> POISS0(rho_in)
           -> VXC0SP(xc_obj, rho_in, ..., B_fsm)
              -> legacy XCPOT_hybrid -> XCPOT
              or libXC LDA wrapper / radial GGA helper
           -> NEWRHO(V, PL, QL, radial solutions)
           -> residual and radial mixing
           -> final accepted iteration
              -> VXC0SP(..., snapshot=atom%radial_ground_state)
                 -> radial_ground_state%capture
              -> radial_ground_state%mark_accepted
  -> bands%calculate_magnetic_moments
  -> report
     -> set_reported_moment
     -> radial_ground_state%write_file
```

`ATOMSC` fixes `n_radial_spin_channels = 2`. This is deliberately independent
of `control%nsp`, which selects the global Hamiltonian mode and is not the
radial density-array channel count. `rho_in` is initialized from
`atom%rho0(2)` and is the density actually passed to `POISS0`, `VXC0SP`, and
the accepted final `NEWRHO` comparison.

The last radial iteration is not captured before mixing. `ATOMSC` first tests
the preceding residual, mixes the density, sets `LAST`, and on the next pass
evaluates `POISS0` and `VXC0SP` on that mixed accepted density. The snapshot is
captured there, before `NEWRHO` produces the next trial density; the newly
computed residual is then retained in `accepted_residual_control` by
`mark_accepted`. Thus XC and density arrays are from the same accepted input
state, not from an initial guess or a pre-mixing trial.

## Radial density semantics, mesh, and Jacobian

`POISS0` documents the live meaning of the array as

```text
RHO = spherical charge density = 4*pi*r*r*RHOTRUE
```

For each local spin channel,

\[
  \mathrm{RHO}_\sigma(r)=4\pi r^2 n_\sigma(r),
  \qquad
  n_\sigma(r)=\frac{\mathrm{RHO}_\sigma(r)}{4\pi r^2}
\]

for nonzero mesh points. The origin is regularized using the same `VXC0SP`
linear extrapolation (`RHO0`) used by the production XC evaluation. The
snapshot stores both the original weighted density and the physical density;
it does not replace `RHO` with a differently normalized quantity.

`ATOMSC` and `POISS0` construct the logarithmic mesh with

\[
 r_i=B\left[\exp(A(i-1))-1\right],
 \qquad
 \frac{dr}{di}=A(r_i+B),
 \qquad i=1,\ldots,NR.
\]

The exact code Simpson weight is

\[
 w_1=w_{NR}=\frac13,
 \quad w_i=\frac43\text{ for even }i,
 \quad w_i=\frac23\text{ for odd interior }i.
\]

Therefore the native electron-number integral is

\[
 N_\sigma=\sum_i w_i A(r_i+B)\,\mathrm{RHO}_\sigma(r_i)
          =\int_0^{R_{max}}4\pi r^2 n_\sigma(r)\,dr.
\]

The independent oracle evaluates both forms, `RHO` with the Jacobian and
`4*pi*r**2*n`, with the same explicit weights. `radial_ground_state` exposes
the same quadrature through `charge_integral`, `spin_number_integral`, and
`log_mesh_integral`.

## Channel conventions and mapping

There are two related conventions in the legacy boundary routine. They are
now recorded explicitly rather than left to positional inference.

| location | channel 1 / first | channel 2 / second |
| --- | --- | --- |
| `ATOMSC`/snapshot `RHO(:, :)` | up | down |
| `VXC0SP` output arrays `V(:, :)`, `vxc_*_radial` | up | down |
| historical `XCPOT` arguments `RHO1,V1` and `RHO2,V2` | down | up |
| libXC wrapper input order | up | down |
| global collinear coefficient-space spin blocks (`spin_off`) | block 1, up | block 2, down |

`VXC0SP` passes its local column 2 as historical `RHO1` (down), local column 1
as historical `RHO2` (up), and maps the returned historical `V2` to local
up/output column 1. The libXC route has an explicit up/down interface, so no
historical swap is applied inside `xcpot_libxc_gga_radial`.

The native snapshot convention is written in every file as
`up=RHO(:,1)=Vxc(:,1); down=RHO(:,2)=Vxc(:,2)`.

## `NEWRHO` density construction

For each orbital channel and spin, `NEWRHO` obtains the live radial solution
and energy derivatives from `RSEQSR`/`PHDFSR`, then accumulates

```text
RHO += Q0*(GFAC*G**2 + Gs**2)
     + 2*Q1*(GFAC*G*GP + Gs*GPs)
     + Q2*(GFAC*(GP**2 + G*GPP) + GPs**2 + Gs*GPPs)
```

where `Q0=QL(1)`, `Q1=QL(2)`, and `Q2=QL(3)`. `GFAC` is the scalar-relativistic
large-component metric factor; the small component has unit weight. The
source explicitly identifies `FUN2` as probability density rather than charge
density, while the accumulated `RHO` is the spherical charge-density array
consumed by `POISS0` and `VXC0SP`. No state-wise orbital reconstruction is
used as a surrogate for the accepted live `RHO`.

## Spin moment convention

The native magnetic DOS path forms

```text
dz = -aimag(g_uu - g_dd)/pi
```

and integrates it to `potential%mz`; the x/y projections use the corresponding
off-diagonal spin quantities. In the collinear z state this is the positive
up-minus-down number convention:

\[
  m_s = N_\uparrow-N_\downarrow,
  \qquad
  M=\mu_B(N_\uparrow-N_\downarrow).
\]

The report and snapshot store the numerical value in `mu_B` units, with the
factor `\mu_B` understood by the unit label. No extra electron-magnetic-moment
minus sign is inserted. `integrated_moment_muB` is the radial spin number
under this native convention; `reported_moment_xyz_muB` is the independently
reported DOS moment. The XC field is not derived from this reported moment.

For a future response consumer, `n_up-n_down` is a number-density difference
in `bohr^-3`; multiplication by `mu_B` is required when a magnetic-moment
density is desired. A physical electron-moment sign convention must be an
explicit downstream conversion, not silently mixed into `Bxc`.

## Exact multiplicative Vxc source

The snapshot takes `vxc_up_radial` and `vxc_down_radial` directly from the
active `VXC0SP` evaluation. It does not use `potential%vrmax`, `C`, a boundary
value, or any compressed LMTO potential field as a Vxc surrogate.

* Legacy selectors use the existing `XCPOT_hybrid`/`XCPOT` implementation.
  The historical down/up argument order is mapped back to the canonical local
  up/down arrays at the `VXC0SP` boundary.
* A libXC LDA route uses `xcpot_libxc_wrapper`. It maps RS-LMTO down/up inputs
  to libXC `[up,down]`, accumulates the selected native components, converts
  libXC Hartree outputs to internal Rydberg exactly once with `2*`, and maps
  outputs back to the canonical local channels.
* A libXC GGA route uses `xcpot_libxc_gga_radial`. It evaluates `vrho` and
  `vsigma`, constructs the spherical radial flux and its divergence, and then
  forms the multiplicative radial potential. The explicit final `2*` converts
  the Hartree-valued libXC result to Rydberg after the radial functional is
  assembled.

The snapshot retains the resulting exact pointwise arrays after the local
automatic `VXC0SP` temporaries have expired.

## \(\Delta V_{xc}\), Pauli field, and constraining field

The convention-neutral fields are defined pointwise by

\[
 \Delta V_{xc}=V_{xc}^{\uparrow}-V_{xc}^{\downarrow},
 \qquad
 V_{xc}^{0}=\frac{V_{xc}^{\uparrow}+V_{xc}^{\downarrow}}2,
 \qquad
 B_{xc}^{\rm Pauli}=\frac{\Delta V_{xc}}2.
\]

Because global block 1 is up and block 2 is down,

\[
 \begin{pmatrix}V_{xc}^{\uparrow}&0\\0&V_{xc}^{\downarrow}\end{pmatrix}
 =V_{xc}^{0}I+B_{xc}^{\rm Pauli}\sigma_z,
 \qquad \sigma_z=\operatorname{diag}(1,-1).
\]

The historical reversed quantity
`Bxc_historical_reversed=(Vxc_down-Vxc_up)/2` is not used by this contract.
The snapshot field `bxc_pauli` has the explicit Pauli sign and factor.

`B_fsm` is a separate scalar radial constraining/finite-spin-moment field. The
live insertion is

```text
V(:,1) += Vxc_up   + B_fsm
V(:,2) += Vxc_down - B_fsm
```

so the total Kohn-Sham spin splitting is

\[
 V_{KS}^{(1)}-V_{KS}^{(2)}
 =\Delta V_{xc}+2B_{fsm},
\]

with the common Hartree/nuclear scalar part canceled. The snapshot stores
`total_ks_spin_splitting`, `delta_vxc`, and `constraining_field_ry`
separately. The unit supplement changes only the synthetic constraining field
and verifies that the captured Vxc arrays are invariant while the total
splitting shifts by twice the field.

## Units

| quantity | unit / meaning |
| --- | --- |
| `r`, `rmax`, `B` | bohr |
| `A` | dimensionless logarithmic mesh increment |
| `n_up`, `n_down` | electrons bohr`^-3` |
| `RHO` weighted density | electrons bohr`^-1`; `RHO dr` integrates to electrons |
| `N_up`, `N_down`, `N_up-N_down` | electron-number units |
| reported spin moment | `mu_B` under native positive up-minus-down convention |
| `Vxc`, `delta_vxc`, `bxc_pauli`, `B_fsm`, total splitting | Ry |
| libXC raw `exc`, `vrho`, `vsigma` | Hartree at the API boundary, converted once to internal Ry |

The radial mesh is not an Ångström mesh. Lattice inputs may use their own
documented units; those units are not copied into this radial contract.

## XC provenance

`radial_xc_provenance` is copied from the live `xc` object while the active
XC route is known. It records:

* backend (`legacy RS-LMTO` or `libXC`), `TXCH`, `TXC`, functional name, and
  mapping-quality label;
* internal energy units (`Ry`);
* libXC route family and `nspin`;
* selected native component IDs and their family/kind metadata when libXC is
  active; and
* whether the route is spin-polarized and whether radial GGA evaluation was
  used.

The supported bundle TXC=101 (`XC_LDA_X + XC_LDA_C_VBH`) was run through the
live radial path with libXC enabled. The existing selector mapping document
[`XC_SELECTOR_AND_LIBXC_MAPPING.md`](XC_SELECTOR_AND_LIBXC_MAPPING.md) remains
the authority for the full selector namespace and equivalence labels. A
response consumer must use the provenance from the same XC state as the
arrays; it must not infer a functional from a compressed potential name.

## Snapshot and data lifetime

Before this audit, the exact `RHO` and multiplicative Vxc arrays were automatic
locals of the `ATOMSC`/`VXC0SP` path and were unavailable after that call. The
minimal accepted-state fix is:

* `symbolic_atom` owns one `radial_ground_state` component;
* the optional `VXC0SP(..., snapshot=...)` argument copies the exact accepted
  `RHO`, mesh, physical densities, Vxc arrays, total spin splitting, and field;
* `ATOMSC` marks it accepted with the final inner iteration and residual; and
* `report` adds the native magnetic moment and writes
  `radial_ground_state_<symbol>_<site>.dat` after all atomic temporaries have
  expired.

The in-memory accessor for a future consumer is
`symbolic_atom%radial_ground_state`. The file is an auditable persistence seam,
not the source of a second calculation. `valid` and `accepted` prevent a
consumer from treating an uninitialized or pre-acceptance object as a
ground-state contract.

## Numerical oracle: converged magnetic bcc Fe

The primary fixture is `tests/scf/cases/bulk/bccFe`, with TXC=1, two local
radial channels, and a real outer SCF convergence. The oracle uses the emitted
file after process completion, independently parses all radial rows, and
checks charge, moment, Vxc decomposition, total-field separation, provenance,
and report lifetime.

Legacy TXC=1 result:

| quantity | value |
| --- | ---: |
| accepted inner iteration | 58 |
| accepted residual control | `3.9856590649e-7` |
| independent `N_up` | `14.0520296321884` |
| independent `N_down` | `11.9479778001210` |
| independent total charge | `26.0000074323093` |
| independent spin number | `2.1040518320674` |
| reported `m_z` | `2.104059` `mu_B` |
| maximum `|Vxc_up-Vxc_down|` | `2.528595e-1` Ry |
| `B_fsm` | `0` Ry |

Selected legacy rows show that the snapshot contains the local density and the
pointwise multiplicative spin splitting, rather than only a compressed scalar:

| row | `r` (bohr) | `n_up` | `n_down` | `delta_vxc` (Ry) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | `0` | `8.67890789e3` | `8.67971334e3` | `0` |
| 2 | `2.75304689e-6` | `8.46952391e3` | `8.47030993e3` | `1.54142240e-3` |
| 248 | `1.89116664e-2` | `2.31053742e3` | `2.31075318e3` | `9.99385824e-4` |
| 495 | `2.6622` | `1.71842390e-2` | `1.98018104e-2` | `2.09601468e-2` |

The same live oracle with TXC=101 and libXC 5.2.3 recorded backend `libXC`,
functional `Slater exchange + von Barth & Hedin`, family `1`, `nspin=2`, and
native component IDs `[1,17]`. It returned
`N_up=14.0520296605947`, `N_down=11.9479777717116`, spin number
`2.1040518888832`, and reported `m_z=2.104059` `mu_B`.

The six-decimal report value differs from the full radial integral by less
than `8e-6` in both live runs; the oracle tolerance is `2e-4` because the
report is intentionally printed at six decimals.

## Test matrix and evidence status

| test / evidence | what it proves | status |
| --- | --- | --- |
| `tests/validation/val22_lr_radial_ground_state.py` | live converged magnetic bcc-Fe oracle; independent charge/moment quadrature; accepted-state lifetime; exact pointwise decomposition; no compressed Vxc surrogate | PASS |
| `Val22LrRadialGroundState` | CTest registration and legacy TXC=1 execution | PASS |
| `Val22LrRadialGroundStateLibxc` | same live radial snapshot contract for supported libXC TXC=101, including backend/provenance and native IDs | PASS with libXC enabled |
| `tests/unit/test_lr_radial_ground_state.f90` / `UnitLrRadialGroundState` | independent logarithmic-mesh quadrature, lifecycle, Pauli sign/factor, and nonzero-field separation supplement | PASS; algebraic supplement |
| `UnitLrBasisRadial`, `UnitLrBasisAngular`, `UnitLrBasisAugmentation` | previously certified live radial augmentation and angular/Hamiltonian contracts remain intact | PASS |
| `UnitXcSelectorSemantics`, `UnitLibxcXcBaseline`, `UnitLibxcGgaRadial`, selector negative guards | XC namespace, mapping, libXC metadata, LDA/GGA route behavior, and fail-closed unsupported paths | PASS |
| `source/self.f90`, `source/xc.f90`, `source/potential.f90` audit | exact call path, RHO/Jacobian semantics, channel order, source Vxc, Ha-to-Ry conversion, and distinction from compressed `vrmax` | PASS by source inspection plus oracle |
| SOC/noncollinear spinor radial contract | no false extension of the two-channel collinear contract | DEFERRED |
| TDDFT response / `chi0` | explicitly out of scope | NOT IMPLEMENTED |

The focused libXC run was made possible by the local libXC 5.2.3 installation.
On builds with `ENABLE_LIBXC=OFF`, the libXC CTest is not registered; the
legacy oracle and the existing fail-closed libXC selector tests remain valid.

## Completion checklist

- [x] Exact branch and starting HEAD recorded.
- [x] LR-GF-01 and LR-BASIS-00 evidence/tests checked before the audit.
- [x] Live `ATOMSC -> POISS0 -> VXC0SP -> NEWRHO` path traced.
- [x] `RHO` meaning, origin handling, logarithmic Jacobian, and Simpson rule recorded.
- [x] Up/down channel order and global `spin_off` block mapping recorded.
- [x] Native `M=mu_B(N_up-N_down)` convention recorded separately from XC fields.
- [x] Exact multiplicative Vxc arrays retained from the accepted state.
- [x] Legacy and libXC mapping, GGA radial divergence, provenance, and Ha-to-Ry conversion checked.
- [x] `Delta Vxc`, Pauli `Bxc` sign/factor, and separate `B_fsm` decomposition tested.
- [x] Real converged magnetic oracle plus independent charge and moment quadratures added.
- [x] Snapshot lifetime and post-process persistence tested.
- [x] No TDDFT response implementation added.
