# Full-spinor Liechtenstein/FLL functional: independent algebra oracle

This is an algebraic R08 contract and oracle result, not evidence that the
production SOC+U path is enabled or physically validated. The orbital basis is
the complex spherical d shell `m=-2,...,+2`; the ten-component index is
spin-major (`up m=-2..2`, then `down m=-2..2`). The density is the standard
operator matrix

\[
P_{a\sigma,b\tau}=\langle c^\dagger_{b\tau}c_{a\sigma}\rangle,
\qquad P=P^\dagger,
\]

so all four spin blocks, including orbital coherences in the spin-offdiagonal
blocks, are retained. `T_abcd` follows the repository's independent complex
spherical tensor reference convention (`<a c|b d>`).

## Hartree–Fock interaction and derivative

For a spin-independent Coulomb tensor, the full spinor Hartree–Fock energy is

\[
E_{\rm int}=\frac12\sum_{abcd,\sigma\tau}T_{abcd}
\left(P_{b\sigma,a\sigma}P_{d\tau,c\tau}
-P_{d\tau,a\sigma}P_{b\sigma,c\tau}\right).
\]

The first term is direct and the second is exchange. Spin indices are summed
independently; exchange therefore includes terms involving `P_up,down` and
`P_down,up`. The derivative is defined by
\(\delta E=\operatorname{Re}\operatorname{Tr}(V^\dagger\delta P)\), for
Hermitian variations. This convention avoids an implicit transpose ambiguity
between the density and Hamiltonian matrices.

The FLL correction audited here is the spin-density (sFLL) extension of the
existing collinear expression. Define

\[
D_{\sigma\tau}=\sum_a P_{a\sigma,a\tau},\quad N=\operatorname{Tr}D,
\qquad E_{\rm DC}=\frac U2N(N-1)-\frac J2[\operatorname{Tr}(D^2)-N].
\]

It reduces for collinear `D=diag(N_up,N_down)` to
\(UN(N-1)/2-J\sum_\sigma N_\sigma(N_\sigma-1)/2\), the current
spin-dependent FLL convention. Its derivative has orbital-diagonal spin
blocks
\[
(V_{\rm DC})_{a\sigma,a\tau}
=[U(N-\tfrac12)+\tfrac J2]I_{\sigma\tau}-J D_{\sigma\tau}.
\]
The Hubbard correction is `E_int-E_DC`, with potential `V_int-V_DC`.

This spin-density extension is an explicit convention. Collinear reduction
and global spin-rotation covariance alone do not select a unique
noncollinear/SOC double-counting functional. The literature distinguishes
charge-only cFLL from spin-density sFLL and notes that the choice depends on
the exchange-correlation functional; see [Han et al., *Scientific Reports*
8, 9559 (2018)](https://doi.org/10.1038/s41598-018-27731-4). Full spinor
occupations and offdiagonal spin potentials in rotationally invariant DFT+U
are discussed by [Kozub et al., *Physical Review B* 94, 125113
(2016)](https://doi.org/10.1103/PhysRevB.94.125113).

## Required identities and tested evidence

At `J=0`, the monopole tensor must give the full ten-dimensional Dudarev
identity for an arbitrary Hermitian spinor density:

\[
E_{\rm Hub}=\frac U2\operatorname{Tr}(P-P^2),\qquad
V_{\rm Hub}=U(\tfrac12 I_{10}-P).
\]
This is stronger than applying a scalar formula to two separate 5×5 blocks:
the identity acts on the full 10×10 matrix and includes its spin coherences.

For a rank-one projector, direct and exchange cancel exactly in `E_int`; at
`J=0`, `E_DC` also vanishes, so the Hubbard correction has no one-electron
self-interaction. For `J>0`, the chosen sFLL term can leave a nonzero remainder
for an orbital-spin-entangled rank-one state: the reduced 2×2 `D` can have
fractional eigenoccupations even when full `P` is a projector. The oracle
reports this remainder separately. It must not be “fixed” by altering the
double counting without a separately justified functional choice.

The runnable oracle is
[`spinor_functional_oracle.py`](../../../tests/validation/fixed_ldau/spinor_functional_oracle.py).
It constructs the d-shell tensor from the independent SymPy complex-harmonic
reference, evaluates the full spinor second-quantized interaction, and checks
the derivative by central differences along a Hermitian direction containing
complex spin-offdiagonal elements at three step sizes. It also checks the
`J=0` Dudarev energy and matrix potential, rank-one direct/exchange and `J=0`
cancellation, collinear sFLL reduction, global SU(2) covariance with SOC off,
and orbital basis relabeling with all tensor indices transformed together.

Run with the validation Python environment:

```text
/opt/homebrew/bin/python3.11 tests/validation/fixed_ldau/spinor_functional_oracle.py
```

The recorded deterministic run (seed `20261014`, `U=2.3`, `J=0.41`,
`F2=14J/1.625`, `F4=0.625F2`) gives `Eint=26.5992101356`,
`EDC=25.1061687628`, and `EHub=1.49304137276` for a Hermitian 10×10 matrix
with spin-coherence Frobenius norm `0.453068319074`. The directional derivative
is `-0.0822676964559`; absolute central-difference errors for steps
`1e-3, 1e-4, 1e-5` are `1.59e-11, 1.21e-10, 3.34e-10`. The J=0 Dudarev
energy and potential residuals are `1.42e-14` and `1.09e-14`. A rank-one
spinor gives `Eint=4.87e-17`, `EDC(J=0)=0`, `EHub(J=0)=5.29e-17`; for this
entangled state `EDC(J=0.41)=0.0722164179803`. SU(2) and orbital relabeling
energy/potential covariance checks pass. These are independent functional
algebra results only (evidence class A for the oracle), not production-SCF
evidence.

## SOC basis convention note

The source basis uses spin-major up/down blocks separated by `spin_off=norb`;
the complex spherical orbital order is `m=-l,...,+l`. In `build_lsham`, the
spin blocks have the standard `xi/2` pattern
`H_up,up=Lz`, `H_up,down=Lx-iLy`, `H_down,up=Lx+iLy`,
`H_down,down=-Lz`, assuming the usual angular-momentum ladder convention.
This source inspection is a convention audit only; it does not validate the
SOC eigenvalues, radial `xi`, spinor occupations, SCF energy, or an enabled
SOC+U calculation. Any production claim still requires SOC-off covariance
and `xi -> 0` limits, direct `j=l±1/2` checks, and actual SOC+U RS/KS tests
after the radial/XC and convergence contracts are closed.

There is also a separate spin-operator sign discrepancy to resolve before
using transverse spin observables as SOC validation. `math.f90` documents
standard `sigma_y=[[0,-i],[i,0]]` in up/down order, but assigns the opposite
offdiagonal signs. These matrices are consumed directly as operators in
spin-current anticommutators and spin-accumulation setup
(`hamiltonian_build.f90`, `build_realspace_spin_operators`; and
`recursion_transport.f90`). Therefore the source-level `S_y` has the opposite
expectation for `(up+i down)/sqrt(2)` under the documented basis convention.
`build_lsham` constructs its SOC blocks from `Lx,Ly` directly and does not use
`S_y`, so that path needs its own angular-momentum spectrum test. No sign was
changed as part of this algebra-only oracle.
