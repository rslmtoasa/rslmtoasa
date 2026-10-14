# Phase 2 — complex occupations and Hamiltonian potential

## Representation and data flow

The fixed-Liechtenstein path now stores the standard operator projector
`P_ab = psi_a psi_b*` in `potential%ldm(l,spin,a,b)`. The orbital indices are
the spherical `m=-l..l` basis used by the onsite Hamiltonian. It is Hermitian
and transforms as `P' = C† P C`; the potential is added to the same Hamiltonian
matrix orientation. `hamiltonian_hubbard.f90` contracts the Phase 1 tensor
`T_abcd=<a c|b d>` with `P_dc`, so the complex matrix derivative is
`dE = Tr(V dP)`. A nonzero-J finite-difference test uses general complex
Hermitian occupations and perturbations and checks both this derivative and
`V=V†`.

The reciprocal eigenvector path forms `psi_a psi_b*`; the real-space path
integrates the Hermitian spectral matrix `i(G_ab-conjg(G_ba))/(2π)`. Both
paths use the shared scalar kernels in `source/math.f90`. Fixed-Liechtenstein
occupations keep their imaginary off-diagonal terms through density mixing,
flatten/expand, and MPI reduction. Existing real restart data remain readable:
`ldm_flatten` remains the real field, `ldm_flatten_imag` is an optional paired
field, and an absent imaginary field is initialized to zero. The optional
field is written only when nonzero, which preserves legacy no-U output bytes.

Both fixed-U producers report per-shell/spin trace, minimum and maximum
eigenvalue of the Hermitian part, and anti-Hermitian residual. They warn when
`max|P-P†| > 1e-6`, an eigenvalue is outside `[-1e-3,1+1e-3]`, or the trace is
outside `[0,2l+1]` by more than `1e-3`. The `1e-6` Hermiticity bound is a
conservative warning threshold rather than a measured integration error
bound; the supported two-step production checks below report zero
anti-Hermitian residual. The `1e-3` bounds are diagnostic tolerances only.
LAPACK eigenvalues are computed from `(P+P†)/2` for diagnosis, while the
original matrix is preserved. No clipping or positivity projection is done.

The ACBN0 proxy remains on its legacy real-Gaunt contraction and rejects a
nonzero imaginary correlated occupation explicitly. SC-U and Hubbard-V
occupation producers retain their previous real behavior. The guards validate
fixed-U for scalar-collinear (`nsp=1`) d-shell calculations with an spd-or-
larger basis. SOC (`nsp=2`), noncollinear (`nsp=3`), p-shell-only, invalid
channel/type indices, and `hoh=.true.` are rejected with diagnostics. This
explains why the original nsp=2 VAL-04 U sweep cannot serve as a fixed-U
validation after the scope guard; no fixture was edited to conceal that fact.

## Tests and runtime evidence

The latest serial Debug build completed successfully in
`/tmp/rslmto-ldau-phase0-build`, executable SHA-256
`1710c0f8325071ed54dfabc55056cadb20656d5125865f7213f453e7da56ad1b`.
The focused LDA+U unit set passed 3/3:
`UnitLdaUTensor`, `UnitLdaUHamiltonian`, and
`UnitLdaUOccupationProducers`. The latter compares KS occupations with an
analytic projector and checks RS quadrature against the finite-window Lorentz
integral, reporting broadening/window error separately from quadrature error.
The Hamiltonian test also reads a legacy namelist without `ldm_flatten_imag`,
checks zero-imaginary recovery, complex flatten/expand round-trip, a Hermitian
complex potential, and the nonzero-J energy derivative. Full unit selection
passed 112/112 tests (130.78 s total; the CTest summary labels 110 tests exactly
`unit`). Reproduction:

```sh
cmake --build /tmp/rslmto-ldau-phase0-build -j4
ctest --test-dir /tmp/rslmto-ldau-phase0-build \
  -R '^(UnitLdaUTensor|UnitLdaUHamiltonian|UnitLdaUOccupationProducers)$' \
  --output-on-failure
ctest --test-dir /tmp/rslmto-ldau-phase0-build -L unit --output-on-failure
```

`tests/validation/val04_lda_u_inputs.py` passed against the latest binary. It
checks no-U and explicit U=J=0 inactivity, U-only J=0 defaulting, U/J eV-to-Ry
conversion, and the unsupported-mode/index guards. The retained scalar U-only
deck and log are in
`/tmp/rslmto-ldau-phase2-inputs-final-diagnostics/u_only/`. The same latest
binary's U-only RS diagnostics report trace and eigenvalue bounds within
physical ranges and `max_antiherm=0` for s, p, and d shells. For d, the logged
spin-up trace/eigenvalue range is `4.260847/[0.789420,0.946245]`; spin-down is
`2.072287/[0.286728,0.499614]`. Under the
same supported scalar fixture, no-U and U=J=0 `Fe_out.nml` files are byte
identical (SHA-256
`e862a18c431862d6336fe151ee0c54ee370d11853d1f5aa88c216c2f7cdba880`).

An additional no-U replay used the original Phase 0 nsp=2, four-iteration
fixture. It exited 0 and reproduced its saved Phase 0 `Fe_out.nml` byte for
byte (SHA-256
`6ee007e1d156ac890937096ec4d159343c8e239191bbb76e3f6f92827bb3f06e`). This
demonstrates that the fixed-U guard leaves the inactive no-U route intact even
for that input mode.

The production mixer is exercised by `UnitLdaUHamiltonian`: a complex
off-diagonal entry passes through `mix_ldm_linear` and `save_ldm_to('current')`
into the current potential. A separate serialized restart test writes the
imaginary companion with `print_state_formatted(unit=...)`, reloads it with
`build_from_file`, and verifies the complex entry. The legacy restart test
omits the imaginary field and verifies zero initialization.

Six two-iteration scalar runtime smokes are recorded in
`/tmp/rslmto-phase2-smoke/REPORT.md` using binary SHA-256
`759f4586068f3d225f92fe4a4a9513a7dc307832b0d62d5d45e109a1bf4cbc89`:
real-space and reciprocal routes, each with no U, U=J=0, and U=2 eV/J=0. All
six exited 0. Both routes preserve byte identity between no-U and U=J=0;
U=2 eV becomes `0.1469971714457357 Ry`, produces a nonzero onsite split, and
retains small imaginary d-LDM entries (~1.2–1.3e-8). The two-step runs are
finite smoke checks, not convergence claims. A synthetic producer test checks
a much larger imaginary coherence and sign.

The original Phase 0 source fixtures and the user-owned `cleanup/` contents
remain unchanged. Runtime inputs, logs, and outputs are outside the repository.

## MPI gate and limits

The isolated Debug MPI build completed with GNU Fortran 16.2/mpifort, MPI
Fortran 3.1, and OpenMP disabled. The committed
`tests/validation/fixed_ldau/mpi_ldau_projector.f90` harness uses the full
8-component sp basis and tests an imaginary p-shell projector through the
production reciprocal occupation call, real/imaginary MPI reductions, and
flatten/expand. One- and two-rank outputs agree to `1.11e-16`; both retain
`P12≈-0.5i`. The separately run scalar real-space U=2 fixture also gives
byte-identical `Fe_out.nml` at one and two ranks, including all 150 imaginary
entries; SHA-256 is
`63ec65b8f32c84160036c34f59f8fdb1a9aca6139cd1fe5879a76956d44e4017`, with zero
anti-Hermitian residual. Reproduction commands and rank logs are described in
[`tests/validation/fixed_ldau/README.md`](../../../tests/validation/fixed_ldau/README.md);
the harness SHA-256 is
`cf92ce8719db25984213a49780109d56e93bb337e00452962b3875b800513aff`.

These scalar runtime smokes stop after two iterations; total-energy
variational consistency and converged physical validation remain later phases.
