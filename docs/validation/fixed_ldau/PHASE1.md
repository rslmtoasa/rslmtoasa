# Phase 1 — orbital basis and Coulomb tensor

## Scope and convention

The new `Coulomb_mat_spherical(l,m1,m2,m3,m4,F0,F2,F4,F6)` returns

\[
T_{abcd}=\int Y_a^*(1)Y_b(1)\,r_{12}^{-1}\,Y_c^*(2)Y_d(2)\,d1\,d2
       =\langle a c|b d\rangle.
\]

With `G_ab(kq)=integral Y_a* Y_kq Y_b`, its implementation sums
`F^k 4π/(2k+1) G_ab(kq) conj(G_dc(kq))`. Thus the selection rule is
`m_a+m_c=m_b+m_d`. In the current stored reciprocal density convention
`nS_ab=psi*_a psi_b`, the fixed Liechtenstein contraction is

`V_ab = sum_cd [T_abcd n_other,cd + (T_abcd - T_adcb)n_same,cd]`.

The energy derivative in this stored orientation is componentwise
`dE=sum_ab V_ab d(nS_ab)`; the conventional projector is `P=nS^T` and
transforms as an operator. Under a unitary change of orbital coordinates
`C`, `nS'=C^T nS C*`, `T'_{abcd}=C*_{ia}C_{jb}C*_{kc}C_{ld}T_{ijkl}`,
and `V'=C†VC`. The FLL scalar correction is diagonal in this basis.

`hcpx('cart2sph')` uses spherical columns ordered `m=-l..l`; its d rows
are `[xy,yz,zx,x2-y2,3z2-r2]`. The reference script explicitly maps the
legacy real-Gaunt orbital order through that Cartesian order before applying
the four-index basis transform. `chbar_nc` applies the same `hcpx` map before
Hamiltonian assembly; the reciprocal LDM producer uses the matching onsite
`m` slices directly. The current LDM is real, so Phase 1 validates its
real-symmetric subset. Complex Hermitian tensor mathematics is tested, but
end-to-end complex occupation storage, mixing, restart, MPI, and Hamiltonian
insertion remain Phase 2 work. The legacy real Gaunt routines remain intact;
the ACBN0 proxy path intentionally remains on its legacy contraction.

## Independent checks

`tests/validation/ldau_tensor_reference.py` uses SymPy's Wigner-3j symbols
and normalized complex harmonics; it neither imports nor calls production
math routines. It also transcribes the old real Gaunt construction and the
actual hcpx basis map so the old/new same-index discrepancy and corrected
basis transform can be measured independently. For `F0=1,F2=.37,F4=.19`,
direct same-index old-real versus spherical maxima are 0.0444 (p) and
0.0636507936508 (d); after the explicit hcpx orbital map the old tensor
matches the independent spherical reference within `4.45e-16`. The allowed
d element `T[-2,-1,-1,-2]=0.0474603174603` and forbidden
`T[-2,-1,-2,-1]=0` discriminate the physical selection rule.

The Fortran unit test covers p/d F0/F2/F4 reference anchors, tensor
symmetries and selection, actual `hcpx` output and random-Hermitian roundtrip,
arbitrary complex U(5) tensor/density/potential covariance, an SO(3) d-shell
rotation, and a nonzero-J finite-difference energy derivative with complex
Hermitian densities and perturbations. FLL terms are included in the finite
difference. A test-helper loop-placement error initially omitted the FLL
diagonal contribution; moving it to an explicit diagonal loop resolved the
test, with the production contraction unchanged.

Reproduction in the Phase 0 build directory:

```sh
cmake --build /tmp/rslmto-ldau-phase0-build --target UnitLdaUTensor UnitLdaUHamiltonian rslmto.x -j4
ctest --test-dir /tmp/rslmto-ldau-phase0-build -R '^(UnitLdaUTensor|UnitLdaUHamiltonian)$' --output-on-failure
/opt/homebrew/bin/python3.11 tests/validation/ldau_tensor_reference.py \
  --fortran-executable /tmp/rslmto-ldau-phase0-build/bin/UnitLdaUTensor \
  --report /tmp/rslmto-ldau-phase1-tensor-reference-final.txt
```

The independent comparison covers 2,118 p/d elements across isolated F0, F2,
and F4 grids. Maximum Fortran/reference residual is `2.09e-17` (tolerance
`2e-12`). The dedicated tensor and existing Hubbard Hamiltonian unit tests
pass. The full unit CTest run passed 111/111 (Phase 0 had 110; Phase 1 adds
the tensor test). The rebuilt executable also passed `tests/validation/val04_lda_u.py`:
U=0/2/4 eV produced the expected no-U identity and U/J response; the eight-step
U=2 eV residual decreased from `2.300011e-2` to a minimum `4.693231e-3`.
This is a finite-iteration smoke check, not a convergence claim.

## Gate 1 disposition

Gate 1 validates the tensor definition, hcpx orbital mapping, corrected
fixed-Liechtenstein tensor use, symmetries/covariance, and the current
real-symmetric active density subset. It does not claim end-to-end support
for imaginary occupations, nor total-energy self-consistency (Phase 3).
