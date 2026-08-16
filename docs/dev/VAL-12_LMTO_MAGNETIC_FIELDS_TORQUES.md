# VAL-12 — LMTO magnetic fields and torques

## Scope and status

This validation establishes the electronic-structure field/torque seam used by
the current production spin-dynamics path. It uses short bcc-Fe SCF/report
cases only; no long trajectory is integrated. The production campaign is:

```text
python3 tests/validation/val12_lmto_fields_torques.py \
  --binary build-rf-debug/bin/rslmto.x \
  --scratch-root <scratch>
```

The result is **Validated (scoped)** for the ordinary one-site bcc-Fe LMTO
field/torque feedback route, including the tested SOC-on/off limits. This is
not a claim that the complete ab-initio spin-dynamics trajectory, thermal
field, STT/SHE, induced-moment, or constraint-plus-SOC dynamics is validated.

## Formulation map

### 1. LMTO local electronic field and torque — active SD route

`bands%calculate_magnetic_torques()` calls the projected-DOS path and, for each
site, constructs the LMTO local vector

\[
 I_i = p_{0,i}\,m_{0,i} - p_{1,i}\,m_{1,i},
\]

where `p_0` and `p_1` are formed from the spin-resolved LMTO `c` and `srdel`
parameters for the active `d` channel (or the highest available channel).
The code stores

\[
 B_i^{\rm LMTO} = I_i\,(\mathrm{Ry}\rightarrow\mathrm{T}),
 \qquad
 \tau_i = \hat m_i\times I_i\,(\mathrm{Ry}\rightarrow\mathrm{T}).
\]

This is the active Antropov-style local electronic-field/torque formulation
in this code. It is a local LMTO force expression, not the LKAG intersite
`J/D/A` exchange tensor and not the SOC torque-correlation operator used for
damping. The output array is named `bands%mag_for` and is logged as “Magnetic
field”; the same routine logs the derived “Magnetic torque”.

The internal `I_i` is an energy coefficient in Ry. `mag_for` and the logged
torque are converted with `ry2tesla = 2.35051754997e5 T/Ry`. The moment in the
cross product is the unit site direction `potential%mom`; the report's raw
spin-moment projections are not substituted for this direction.

### 2. Spin-dynamics field hand-off — active SD route

`processing='sd'` performs the following exact sign hand-off:

```text
bands%mag_for = B_LMTO
field_in     = -bands%mag_for
spin_dynamics%beff = -field_in = B_LMTO
```

The Depondt library receives `beff` in Tesla. The explicit Euler helper uses
`-gamma * m x field_in`, which is therefore `+gamma * m x B_LMTO` under the
current wrapper convention. The Depondt rotation has the same effective
sign for this seam. This records the implementation convention; it does not
reinterpret `mag_for` as a separate energy-gradient definition.

The lower-level Depondt API also has `btorque` and `she_btorque` arguments, but
the current LMTO `processing='sd'` path passes zero arrays. `asd_jij` is read
as a control value but is not consumed by the current SD driver, so no active
LMTO-to-`J_ij` field route is claimed here.

### 3. Constraining field — active SCF contribution, indirect SD input

VAL-10's constraining-field controller returns `bfield` in Ry as a coefficient
of the onsite Pauli block. It is a Hamiltonian field, not a Tesla-valued
`B_eff`. Its convention is

\[
 H_{\rm con}=\mathbf B_{\rm con}\cdot\boldsymbol\sigma,
 \qquad
 \mathbf B_{\rm con}=-\nabla_{m}E_{\rm con}
\]

for the transverse penalty direction; a positive `B_z` raises the up block.
The field is inserted once into the onsite `m=1` Hamiltonian block shared by
the RS and reciprocal builders. `self%apply_constraints()` updates and stores
it during SCF, and the next `self%run()` in SD consumes the resulting
electronic state. `mag_for` has no separate constraint-field component: any
effect after self-consistency is folded into the electronic field route.
The penalty energy is reported separately.

The clean controlled-canting and finite-difference oracle is the existing
`UnitConstrainingField` test from VAL-10. It checks the restoring direction,
the quadratic penalty, global SOC-free rotation, seed retention, and the
finite-difference penalty derivative. A converged fixed-potential electronic
energy-angle sweep for `mag_for` is not exposed by the current production
output and is not claimed here.

### 4. Rotating-frame density residual torque — active GBT diagnostic, not SD

The `spin_density`/reciprocal-spin-density route stores the density in the
rotating frame. Under `density_policy='constrained_spiral'` it reports

\[
 \tau_{\rm residual}=\hat m_{\rm ref}\times m,
\]

in that local rotating frame. Under `relaxed_reference` this residual is zero
by construction. This is a constraint residual/diagnostic, not the
`bands%mag_for` field consumed by `processing='sd'`. Ordinary RS/reciprocal
site directions are Cartesian; GBT's corresponding constrained quantities are
local rotating-frame quantities. GBT with SOC remains outside the supported
scope.

### 5. SOC torque-correlation operator — active damping/inertia, not SD

`hamiltonian%torque_operator_collinear()` builds the commutator-style

\[
 T^a=[\sigma^a,H_{\rm SOC}],
\]

operator consumed by the Gilbert-damping and magnetic-inertia post-processing
routes. It is an energy/operator derivative of SOC, not the Tesla-valued
`mag_for` field and not a torque array passed to `processing='sd'`. Its
SOC-off limit is zero and is validated separately by VAL-08. SOC does affect
the SD LMTO route indirectly through the self-consistent Hamiltonian and
moments; the code must not impose a global-rotation zero-torque expectation
there.

## Validation evidence

The VAL-12 campaign used `nsp=2`/`soc_scale=0` for the collinear case,
`nsp=3`/`soc_scale=0` for the SOC-free rotation pair, and
`nsp=4`/`soc_scale=1` for the SOC pair. Every case checked the production
`m x B` relation using the saved normalized `Fe_out.nml` direction and the
reported field/torque records.

| case | result |
|---|---:|
| collinear SOC-off transverse torque | `0.0 T` |
| SOC-free global-rotation energy difference | `7.614e-7 Ry` |
| largest production `m x B` reconstruction error | `9.08e-7 T` |
| SOC-on rotated torque norm | `3.548e-3 T` |
| SOC-off rotated torque norm | `3.562e-4 T` |

The SOC-on finite torque is evidence that SOC is treated separately; it is not
an assertion that this one-site fixture supplies a converged magnetic
anisotropy benchmark. The small SOC-free residual is within the one-step,
printed-output envelope and does not alter the zero-torque invariant for the
collinear state.

## Checklist

- [x] Field/torque routes inventoried
- [x] Units documented
- [x] Frame conventions documented
- [x] Sign conventions documented
- [x] Collinear zero-torque limit checked
- [x] Global rotation checked without SOC
- [x] Controlled canting checked through the VAL-10 constrained-field oracle
- [x] Finite-difference energy relation checked where valid (VAL-10 penalty)
- [x] SOC behavior checked separately
- [x] No unrelated spin-dynamics changes made

## Tests and files

Passed before and during this task:

- `UnitAbspinlibIntegrators`
- `UnitConstrainingField`
- `UnitConstrainingFieldSource`
- `Example_bulk_bccFe_sd_smoke`
- `Val07ExchangeFormulations`
- `Val12LmtoFieldsTorques`

Changed files are `tests/validation/val12_lmto_fields_torques.py`, this report,
the VAL-12 CTest registration, the developer map, and the Phase-II evidence
ledger.
