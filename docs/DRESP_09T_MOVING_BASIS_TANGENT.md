# DRESP-09T — LMTO moving-basis / representation tangent

Status: **BLOCKED at the live orthogonalization response**  
Primary classification: **`ORTHOGONALIZATION_RESPONSE_REQUIRED`**

This milestone is restricted to a global `L=0` rigid spin rotation. It does
not implement arbitrary `L>0` ALSDA fields, BES/Halle, a Dyson solve, a loss
spectrum, or a magnon calculation.

## Result

The independent radial connection is implemented and closes against the
finite-angle cross-overlap oracle for all four energy-linearization branches.
The connection is not yet a valid tangent in the accepted orthogonal LMTO
coefficient representation. The missing map is the derivative of the live
orthogonalization/representation transformation, not a fitted residual field.

The accepted 64-k Fe run reports:

| diagnostic | result |
|---|---:|
| fixed-basis SR field norm | `2.3569853e-1` |
| radial basis/representation candidate norm | `2.3750277e-1` |
| native total tangent norm | `2.5566670e-1` |
| reconstructed total norm | `1.0973712e-1` |
| field/basis cross term | `-9.9919129e-2` |
| first-order `h` residual | `4.7051065e-1` |
| first-order `o` residual | `2.0034791` |
| first-order `E_nu` residual | `2.1180467` |
| full second-order residual | `1.0056789` weighted; `1.0757936` max-k |
| endpoint metric-covariant residual | `1.7646379` |
| fixed-basis density mismatch | `4.0820980e-2` |
| candidate corrected density mismatch | `9.4416871e-1` |

The corrected density number is diagnostic only. It is not accepted as a
production correction because the coefficient-space connection has not passed
the orthogonalization gate.

The independent native finite-rotation oracle remains healthy:

| theta | native central-difference residual |
|---:|---:|
| `1.0e-2` | `1.6666583e-5` |
| `5.0e-3` | `4.1666615e-6` |
| `2.5e-3` | `1.0416663e-6` |

This is the expected `theta^2` convergence. The predicted field plus current
moving-basis candidate does not pass the corresponding total-matrix oracle;
the three predicted residuals are approximately `11.63` and do not converge
to zero.

## Radial connection

For a site and orbital `(l,m)`, define the linearized endpoint functions

\[
 X^{(0)}_{s,l}=\phi_{s,l}-\varepsilon_{s,l}^{\rm work}\dot\phi_{s,l},
 \qquad X^{(1)}_{s,l}=\dot\phi_{s,l}.
\]

The augmented scalar-relativistic radial overlap between endpoints `p,q` is

\[
 S^{pq}_{s t}(r)=
 \frac{1}{r^2}\left[
 G_{s,p}G_{t,q}\left(1+
 \frac{l(l+1)}{TMC_sTMC_t r^2}\right)+
 \Phi_{s,p}\Phi_{t,q}\right].
\]

This is the mixed-spin augmentation metric obtained from the accepted
DRESP-09S large, lower, and angular radial components. The `-1/3` factor in
the DRESP-09S transverse *field* bilinear is not inserted into a basis
overlap: in the moving-basis derivative the spin generator acts on the
spin-angular ket before the scalar-relativistic lower metric is contracted.
Using the field insertion factor here would fail the equal-channel limit.

For a positive `y` rotation, the spin-angular derivative is

\[
 \partial_\theta U(0)=
 \begin{pmatrix}0&-1/2\\1/2&0\end{pmatrix}.
\]

Therefore the raw branch connection is

\[
 \Gamma_{\uparrow\downarrow}^{pq}=-\frac12 S_{\uparrow\downarrow}^{pq},
 \qquad
 \Gamma_{\downarrow\uparrow}^{pq}=+\frac12 S_{\downarrow\uparrow}^{pq},
\]

with all same-spin blocks zero. The implementation is in
`source/lr_lmto_basis_connection.f90`.

The finite-angle oracle constructs

\[
 M^{pq}(\theta)=\langle\chi^{(p)}(0)|\chi^{(q)}(\theta)\rangle
\]

from the same radial endpoint integrals and the exact two-by-two spin
rotation, then evaluates

\[
 \Gamma_{\rm FD}^{pq}=
 \frac{M^{pq}(+\theta)-M^{pq}(-\theta)}{2\theta}.
\]

The unit fixture checks `00`, `01`, `10`, and `11`, gives a factor-of-four
error reduction under halving `theta`, and checks the equal up/down radial
limit. No native Hamiltonian result enters this calculation.

## Metric convention and gauge

`Gamma_raw=<chi|d chi/dtheta>` is not the coefficient connection in a
nonorthogonal basis. The covariant coefficient connection is defined by

\[
 S K=\Gamma_{\rm raw}.
\]

Its metric identity is

\[
 \delta S=K^\dagger S+S K.
\]

Only when the relevant representation has `S=I` may one identify `K` with
an anti-Hermitian connection. The fixture verifies this distinction and does
not anti-Hermitian-project the raw matrix.

There is also a basis-connection gauge: a differentiable transformation
within the retained local basis changes the intermediate `K`. The physical
matrix tangent must be reconstructed after the same representation map is
applied on both sides. DRESP-09T does not use gauge freedom to alter the
material residual.

## Live LMTO chain

The source-level chain is:

```text
orthogonal C, Delta, ENU, Q
        |
        v
predls: center_band, width_band, shifted_band, obar
        |
        v
cx/wx/cex/obx and cx0/cx1, wx0/wx1, cex0/cex1, obx0/obx1
        |
        v
ee, obarm, enim, eeo
        |
        v
h(k), o, E_nu
        |
        v
H2(k) = E_nu + h - h o h
```

The native first-order and second-order tangents are independently supplied
by DRESP-08. The fixed field is independently supplied by the complete
DRESP-09S SR `L=0` source. The DRESP-09T candidate then does:

```text
radial Gamma_raw[pq], S[pq]
        |
        v
endpoint contraction H^q X[pq] H^p
        |
        v
K = S^-1 Gamma_raw
        |
        v
delta X_basis = K^H X + X K
        |
        v
delta H_basis by the H2 product rule
```

The endpoint-contracted `S`/`Gamma` pair is not the final accepted
orthogonal coefficient pair: its metric-covariant residual is `1.7646379`.
The live code currently stores the transformed LMTO parameters (`cx`, `wx`,
`cex`, `obx`, `obar`) but does not expose the orthogonalization matrix `X` or
its derivative `delta X` as a representation object. Consequently the radial
connection cannot yet be transported into the production coefficient basis.

This is the specific reason for the classification. A residual such as
`delta H_native-delta H_field` is reported only as a diagnostic and is not
used to construct `delta H_basis`.

## Density side

The fixed-basis SR density observable retains the certified DRESP-09S radial
metric. The candidate uses the same `K` to form the represented density
tangent,

\[
 \delta\rho_{\rm basis}=K^\dagger\rho+\rho K,
\]

and contracts it with the same physical radial observable. It is therefore
not a separate fitted density correction. Since `K` has not passed the live
orthogonalization identity, the resulting `9.4416871e-1` density mismatch is
rejected and the accepted fixed-basis reference remains `4.0820980e-2`.

## Tests

Implemented checks:

1. independent mixed-spin radial connection fixture;
2. finite-angle cross-overlap derivative;
3. all `00/01/10/11` branches;
4. equal up/down radial limit;
5. generalized metric identity and no anti-Hermiticity projection;
6. finite-angle matrix connection fixture;
7. 64-k Fe native `h`/`H2` and density bridge;
8. DRESP-08 and DRESP-09S remain separate regression gates.

The material artifact is generated by
`tests/integration/tddft_driver_smoke/input_dresp09t_fe.nml` and checked by
`tests/validation/dresp09t_fe_artifact.py`.

## Required follow-up

Expose the live LMTO orthogonalization/representation map at the `predls` →
`cx/wx/cex/obx/obar` boundary, including its first derivative under a global
spin rotation. Then transport the independently closed radial `Gamma_raw`
through that map and rerun the same first-order, second-order, density, and
finite-angle gates. No `L>0` ALSDA generalization should start before that
closes.

Primary classification: **`ORTHOGONALIZATION_RESPONSE_REQUIRED`**
