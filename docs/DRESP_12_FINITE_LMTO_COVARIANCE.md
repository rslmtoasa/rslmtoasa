# DRESP-12 — finite-LMTO covariance decomposition of the Ward defect

Starting HEAD: `cd2fd3fd69c42b7b992cfc138683f391e89af39d`

DRESP-12 is a diagnostic-only continuation of DRESP-11. It tests whether the
distributed fixed-basis Ward residual is the representation covariance seam
between direct scalar-relativistic field insertion and the production LMTO
commutator. It does not alter the six-branch basis, the 348-dimensional
Pauli space, the Fréchet response, `Kxc = Bxc/P3`, the raw SR source vertex,
or the P3 target.

## Architecture

```text
                   fixed basis
P3 -> Bxc -> F_SR -------------> deltaH_B
                                  |
                                  v
                              Frechet
                                  |
                                  v
                              delta m_B
                                  |
                                  +----------> r_fixed
                                  |
                                  |
              finite-LMTO connection
                    deltaH_conn
                                  |
                                  v
                              Frechet
                                  |
                                  v
                            delta m_conn
                                  |
                                  |
observable-frame connection      |
      delta m_O -----------------+
                                  |
                                  v
                              exact P3
```

The field-level definition is

```text
deltaH_conn(k) = deltaH_cov(k) - deltaH_B(k)
deltaH_cov(k) = -i [G, H(k)]
deltaH_B(k)   = F_SR[Bxc](k)
```

The identity `deltaH_B + deltaH_conn = deltaH_cov` is checked independently
at every accepted k point. The connection is never fitted or projected into
the compact field basis.

## Response and accounting

`deltaH_conn` is measured with the exact finite-temperature static Fréchet
derivative and the frozen DRESP-10F Pauli measurement. The same perturbation
is independently evaluated with the certified 348-dimensional transition
machinery. The observable term is the frozen DRESP-09Y augmentation-frame
tangent, retaining upper, lower-small, lower-angular, total, and the
production Pauli-observable component separately.

The primary physical-space gate is

```text
r_fixed + delta_m_conn + delta_m_O = 0
r_fixed = delta_m_B - m_P3
```

The compact reconstruction uses only the independently derived covariance
terms:

```text
D m_G = delta_m_conn + delta_m_O
```

No BES/Halle correction, SVD zeroing, rank-one repair, kernel rescaling, or
dynamics is enabled.

## Representation capability boundary

The production endpoint tangent is linear in independent endpoint rotations;
the unit regression therefore closes a two-endpoint sitewise superposition
fixture. This supports sitewise rigid rotations from the existing production
map. An arbitrary `L>0` transverse field inside an ASA sphere is not treated
as a rotation of spherical radial functions: it requires nonspherical
radial/basis response. A general local field is therefore not promoted to a
348×348 connection operator from the one rigid vector.

```text
sitewise rigid rotation       DERIVABLE_FROM_EXISTING_PRODUCTION_MAP
arbitrary L,M within sphere   ARBITRARY_L_COVARIANCE_REQUIRES_NEW_BASIS_RESPONSE
multi-site nonuniform field   DERIVABLE only for sitewise rigid rotations
```

## Artifacts and tests

The accepted-state backend is
`tests/integration/tddft_driver_smoke/input_dresp12_fe.nml` and writes:

```text
/tmp/dresp12_fe_4k.dat
/tmp/dresp12_fe_4k.dat.kpoints.csv
```

The integration artifact checks the field identity, native commutator seam,
Fréchet/compact closure, DRESP-09Y observable term, physical master identity,
compact `D m_G` reconstruction, and the explicit no-correction/no-dynamics
boundary. DRESP-11 remains the authoritative owner of the expensive 348-column
denominator/SVD assembly; DRESP-12 reads its frozen singular spectrum and
target overlaps rather than rebuilding that matrix. The DRESP-09U
representation unit regression includes the two-site endpoint superposition
fixture. DRESP-10F and DRESP-11 remain frozen prerequisite gates.

## Result

The accepted 64-k Fe run is `PASS-B`:

- field identity maximum: `1.0431e-17`;
- independent Fréchet/compact connection residual: `1.8802e-13`;
- frozen DRESP-11 denominator reconstruction: `4.0238e-16`;
- master identity relative residual: `1.4619`;
- compact `D m_G` reconstruction residual: `1.4619`.

Thus the finite-LMTO connection is numerically well-defined, but the currently
available observable/accounting terms do not close the full Ward defect. The
next step remains representation-consistent Goldstone restoration; no
correction or dynamics was run. The machine-readable final report, including
the k-resolved table and frozen DRESP-11 SVD spectrum/target overlaps, is
`/tmp/dresp12_fe_4k.dat` with sidecars `.kpoints.csv`, `.DRESP09U`, and
`.DRESP09Y`.
