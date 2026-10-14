# Phase 3 — functional evaluator and total-energy accounting

## Functional convention and API

The production fixed-Liechtenstein evaluator is
`hamiltonian%evaluate_hubbard_energy(density,potential,e_int,e_dc,e_hub,tr_vp)`.
`density` and `potential` are explicit snapshots with shape
`(ntype,4,2,7,7)`, where the orbital pair is the conventional spherical
operator projector `P_ab=psi_a psi_b*` and operator matrix `V_ab`; scalar
outputs have shape `(ntype,4)` and are spin-summed per type and `s/p/d/f`
channel. The trace uses `sum_ab V_ab P_ba = Tr(VP)`. Explicit snapshots let a
caller report `EHub(P_out)`, `Tr(V(P_in)P_in)`, and the cross quantity
`Tr(V(P_in)P_out)` without confusing the density used to build the current
Hamiltonian with the post-mixing density.

For `T_abcd=<ac|bd>`, the interaction evaluator uses

```text
Eint = 1/2 sum_{s,s',a,b,c,d} T_abcd P^s_ba P^s'_dc
       -1/2 sum_{s,a,b,c,d} T_abcd P^s_da P^s_bc
EDC  = U N(N-1)/2 - J sum_s N_s(N_s-1)/2
EHub = Eint - EDC
```

The evaluator uses the same Slater-parameter conversion as the fixed-U
potential, and rejects ACBN0. It returns energy components in the code's Ry
units. The method is a diagnostic/evaluation API; it does not change the
SCF total-energy assembly.

## Algebra and covariance evidence

`UnitLdaUHamiltonian` exercises the production potential and energy evaluator.
The p-shell nonzero-J and d-shell J=0/nonzero-J complex Hermitian directional
derivatives compare `Tr(V delta P)` with central differences at steps
`1e-3`, `1e-4`, and `1e-5`. The J=0 d-shell check also compares against the
independent Dudarev trace expression
`U/2 sum_s [Tr(P_s)-Tr(P_s^2)]`; a normalized complex rank-one integer
projector has zero `Eint`, `EDC`, and `EHub`. The d-shell production
potential transforms covariantly under a physical SO(3) rotation, and
`EHub`/`Tr(VP)` remain invariant.

`tests/validation/fixed_ldau/phase3_d_energy_reference.py` is a separate
SymPy Wigner-3j construction using normalized complex Gaunt coefficients. For
the d-shell complex fixture (`U=1.4 Ry`, `J=0.22 Ry`) it gives
`Eint=6.530967462948`, `EDC=5.334218750000`, `EHub=1.196748712948`,
`Tr(VP)=0.402247425896`, and
`V_up(1,2)=-0.028236422844+0.004210279647i`. The Fortran unit test asserts
these independent values. The J=0 analytic fixture (`U=2.1 Ry`) gives
`Eint=4.390487535`, `EDC=2.8709690625`, `EHub=1.5195184725`, and
`Tr(VP)=0.700161945`. A J=0 integer rank-one projector check confirms zero
self-interaction.

Reproduce focused checks with:

```sh
cmake --build /tmp/rslmto-ldau-phase0-build \
  --target UnitLdaUHamiltonian UnitLrBasisAugmentation -j4
ctest --test-dir /tmp/rslmto-ldau-phase0-build \
  -R '^(UnitLdaUHamiltonian|UnitLrBasisAugmentation)$' --output-on-failure
python3.11 tests/validation/fixed_ldau/phase3_d_energy_reference.py
```

At the current Phase 3 checkpoint the focused `UnitLdaUHamiltonian` and
`UnitLrBasisAugmentation` CTests pass. The final full unit command
`ctest --test-dir /tmp/rslmto-ldau-phase0-build -L unit --output-on-failure`
reported 112/112 tests passed (110 tests carry the `unit` label); its captured
stdout and exit marker are in `/tmp/phase3-final-full-unit.log` (SHA256
`648150307f560e4e4df7497a851acea4f66eba0f27a015379a84fe7004eb8d78`). The radial replay
reports zero change in `SUMEV`, `SUMEC`, and density when only U/J metadata
changes with radial potential, QL, and PL held fixed. This test changes
metadata only; it does not inject a Hubbard Hamiltonian or potential into
NEWRHO.

A separate carrier sensitivity test uses production POTPAR-generated
`C/QPAR/SRDEL`, inverts `calculate_pl` to calibrate a synthetic baseline to
`PNU` (maximum baseline PL error `1.51e-13`), and then perturbs the d-shell
gravity center by 0.025, 0.0125, and 0.00625 times `Tr(VP)_s/N_s`. The test invokes
production `calculate_pl` and passes its PL output through unchanged
`NEWRHO`. Although the focused test passes, the response is strongly
nonlinear: at scale 0.025, d-shell PL changes are
`[0.600698,0.636219]`, `SUMEV` changes are `[4.58404,1.83087] Ry`, and the
ratios to injected `Tr(VP)` are `[7.82e3,2.40e3]`; at scale 0.00625, the
corresponding values are `[0.436709,0.596751]`, `[3.00385,1.58337] Ry`, and
`[2.05e4,8.32e3]`. These observations do not identify a physical coefficient
or sign. This is only a parameterization/path sensitivity test on a calibrated
synthetic radial setup, not a self-consistent spectral-to-radial energy
identity. The earlier attempt using `gravity_center=ENU` failed its baseline
calibration and is not used as evidence. The latest focused test stdout is
reproducible in CTest's `Testing/Temporary/LastTest.log`; the tested source
hash is `524198b8653ecc8d27d63ecad4cf206d5b5aa3f9e0ba70539af0d60e883c2647`.

## Gate 3 disposition

The functional energy and potential algebra is independently checked, but
the physical total-energy adjustment remains unresolved. The carrier test's
nonlinear response is insufficient to establish whether production radial
`SUMEV` contains all, part, or none of `Tr(VHub P)` after the fixed-U
spectral moment is mapped through `gravity_center`, `PL`, and `NEWRHO`. No
`EHub` or `EHub-Tr(VP)` term has been added to `physical_total_energy`.
Therefore variational consistency and the Phase 3 release gate are **not
passed**.

The radial projector is the code's ASA/LMT orbital moment expansion
(`QL`/`PL`); this phase has not established equivalence to a full PAW/overlap
metric projector. Only the supported fixed-Liechtenstein scalar-collinear
mode is considered; ACBN0 and unsupported spin/orbital modes remain outside
scope.
