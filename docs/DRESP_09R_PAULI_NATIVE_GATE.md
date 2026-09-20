# DRESP-09R Pauli projected native-rotation gate

DRESP-09R adds an independent Pauli projected direct radial/angular source
operator.  It uses only the accepted `phi_large`, `phidot_large`, and
`enu_work` arrays, emits the explicit `B00`, `B01`, `B10`, and `B11` endpoint
branches, and covers the complete `L=0..4` complex response space.  The
oracle does not call `evaluate_pauli_transition_vertex` or the DRESP-09
compact metric-adjoint routine.

The accepted bcc-Fe material gate uses the production radial snapshot and
64-point `4x4x4` k-space handoff.  The field is the normalized `L=0` field
`sqrt(4*pi)*radial_ground_state%bxc_pauli`, with the DRESP-07/08 transverse
rotation sign and factor.  It compares the combined sigma-plus/minus direct
operator with the DRESP-08 native tangent `-i[G,H2(k)]`, and reports per-k
relative Frobenius, weighted RMS, maximum element, occupied action, and
near-Fermi action.

| Capability | Result | Evidence |
|---|---|---|
| Pauli direct radial/angular operator | Certified | `UnitDresp09PauliDirect` |
| Explicit B00/B01/B10/B11 endpoint orientation | Certified | `UnitDresp09PauliDirect` |
| Complex mixed `L/M`, all `L=0..4` | Certified | `UnitDresp09PauliDirect` |
| Pauli direct vs compact adjoint | Closed; `8.56e-16` max relative | `/tmp/dresp09r_fe_4k.dat` |
| Pauli-projected `bxc_pauli` vs native tangent | Blocked; `1.637e-1` max relative, `1.191e-1` weighted RMS | `/tmp/dresp09r_fe_4k.dat` |
| Full scalar-relativistic mixed-spin route | `BLOCKED - MIXED-SPIN RADIAL METRIC NOT CERTIFIED` | Existing DRESP-09 capability gate |
| Ward / explicit ALSDA follow-up | Not run | Native gate failed as required |

The material verdict is
`PAULI_PROJECTION_INSUFFICIENT_FOR_NATIVE_TANGENT`.  The mismatch is a gate
result, not a license to alter `chi0`, `Kxc`, POTPAR differentiation, radial
SCF, BES/Halle, or the full-SR blocker.
