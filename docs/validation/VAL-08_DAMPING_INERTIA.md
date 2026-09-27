# VAL-08 — Damping and magnetic moment of inertia

## Scope and status

This validation uses the small metallic bcc-Fe one-site-pair fixture from
`tests/regression/triad_bccFe_damping`, with collinear `nsp=2`, `hoh=.true.`,
and the supported on-site p/d SOC terms. The existing recursion/Lehmann/Dyson
damping triad is retained. Campaign runs are outside the quick gate:

```text
python3 tests/validation/val08_damping_inertia.py \
  --binary build-rf-debug/bin/rslmto.x \
  --scratch-root <scratch>
```

Damping is **Validated (scoped)** for this route/material/operator envelope;
this does not define an intrinsic metallic alpha at one arbitrary broadening.
The moment-of-inertia quantity is **Experimental**.

## Formula and conventions

The damping routine builds the SOC-derivative torque operators

$$
T^a=[\sigma^a,H_{SOC}],
\qquad A_{ij}=G_{ij}-G_{ji}^{\dagger},
$$

and contracts the resulting torque products at the energy-grid point nearest
the Fermi level. Its scalar prefactor is

$$
-\frac{1}{2\pi m_i},
$$

where `m_i` is the spin moment accumulated from the spin-resolved `ql`
charges. Thus damping is a Fermi-point estimator here, not an energy-integral
whose convergence can be summarized by a single integrated value. The relevant
energy-grid check is whether the nearest-(E_F) point is stable.

The inertia diagnostic uses the same torque operators and defines

$$
B_{ij}=G_{ij}+G_{ji}^{\dagger},
\qquad
sB_{ij}=\frac{d^2 B_{ij}}{dE^2},
$$

with the central second difference on the uniform production energy mesh. The
raw component is the real or imaginary trace of

$$
T_i^a A_{ij}T_j^{b\dagger}sB_{ji}
+T_i^a sB_{ij}T_j^{b\dagger}A_{ji}.
$$

`allinertias.out` reports this raw tensor at the Fermi-point mesh index;
`inertia-energy.out` reports the summed real tensor over the energy mesh.
There is no independently implemented production relation, published
normalization wired into this code, or legitimate equality to J/D/A or
Gilbert damping. The common damping-related limit is the SOC-off limit: both
quantities vanish because the torque operator vanishes.

The production exchange output labelled `Iij` is the anisotropic exchange
tensor, not this dynamical moment of inertia, so it is not an independent
inertia oracle.

## Damping evidence

The values below are alpha = \(\frac12(\alpha_{xx}+\alpha_{yy})\). All tested
on-site damping tensors had zero antisymmetric part to the printed precision;
the smallest principal symmetric minor was positive. The Lehmann/Dyson
reference values agree to the existing \(10^{-5}\) absolute route pin.

| sweep | alpha |
|---|---:|
| recursion, (N_k=8^3, \eta=0.02) | 0.001341155 |
| Lehmann, (N_k=4^3, \eta=0.02) | 0.004156500 |
| Lehmann, (N_k=8^3, \eta=0.02) | 0.002682953 |
| Lehmann, (N_k=12^3, \eta=0.02) | 0.002527619 |
| Lehmann, (\eta=0.01, N_k=8^3) | 0.001866412 |
| Lehmann, (\eta=0.02, N_k=8^3) | 0.002682953 |
| Lehmann, (\eta=0.04, N_k=8^3) | 0.003491698 |
| Dyson, (\eta=0.02, N_k=8^3) | 0.002682953 |
| Lehmann, (N_E=120,240,360), (N_k=8^3,\eta=0.02) | 0.002682953, 0.002682953, 0.002682953 |
| SOC off, Lehmann, (N_k=8^3,\eta=0.02) | 0.000000000 |

The positive symmetric-part minimum over the nonzero-SOC damping sweep was
\(1.619\times10^{-6}\); the largest printed antisymmetric component was zero.
The \(N_E\) sweep is stable because this implementation samples the same
nearest-\(E_F\) point, not because an energy-integrated damping has converged.
The eta and k sweeps demonstrate why no single-eta alpha is promoted as an
intrinsic material value.

## Inertia evidence

The Fermi-point raw tensor was finite and symmetric in every run, with the
largest real-tensor antisymmetric component below \(4\times10^{-15}\). This is
a symmetry/implementation check, not a physical normalization check.

| sweep | max absolute raw tensor component |
|---|---:|
| Lehmann, (N_E=120, \eta=0.02) | 18.4738058 |
| Lehmann, (N_E=240, \eta=0.02) | 41.8893563 |
| Lehmann, (N_E=360, \eta=0.02) | 48.9303019 |
| Lehmann, (\eta=0.01, N_E=240) | 108.353528 |
| Lehmann, (\eta=0.02, N_E=240) | 41.8893563 |
| Lehmann, (\eta=0.04, N_E=240) | 11.090465 |
| recursion, (N_E=240, \eta\approx0) | 14.5634034 |
| SOC off, Lehmann, (N_E=240) | 0.0000000 |

The \(N_E=120\to240\to360\) changes in the max absolute component are
23.4155505 and 7.0409456, while the broadening sweep spans 11.090465 to
108.353528. The raw inertia is therefore not converged in this campaign and
is not promoted. In particular, the fact that the code now writes a finite
number is not treated as evidence of a production magnetic moment of inertia.

## Checklist

- [x] Damping triad retained.
- [x] Eta dependence measured.
- [x] K-mesh dependence measured.
- [x] Damping tensor conventions, sign, and symmetry checked.
- [x] Inertia formula and assumptions documented.
- [x] Inertia energy and broadening dependence measured.
- [x] SOC-zero and tensor symmetry behavior checked.
- [x] No single-eta result overinterpreted.
- [x] Maturity updated separately: damping **Validated (scoped)**; inertia **Experimental**.
