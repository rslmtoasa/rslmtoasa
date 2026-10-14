# Fixed-U MPI checks

`mpi_ldau_projector.f90` is a small integration harness for the production
reciprocal fixed-Liechtenstein occupation path. It splits two equally weighted
k points across ranks; each carries an occupied p-shell projector with
`P(1,2) = -i/2`. The call to
`reciprocal%calculate_ldm_from_projected_dos(...,.true.)` executes the
production `sync_lattice_ldm` real and imaginary MPI reductions. The harness
then exercises the potential LDM flatten/expand storage round trip.

In this macOS/GNU Fortran environment, configure and build the isolated MPI
tree with:

```sh
cmake -S . -B /tmp/rslmto-ldau-phase2-mpi-build \
  -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON -DRUN_UNIT_TESTS=ON \
  -DENABLE_MPI=ON -DENABLE_OPENMP=OFF \
  -DCMAKE_Fortran_COMPILER="$(command -v mpifort)" \
  -DCMAKE_C_COMPILER=/usr/bin/cc
cmake --build /tmp/rslmto-ldau-phase2-mpi-build --target rslmto.x
```

Compile the harness against that build's modules and libraries (the remaining
link inputs are the same as `UnitLdaUHamiltonian.dir/link.txt`):

```sh
BUILD=/tmp/rslmto-ldau-phase2-mpi-build
SPGLIB=/opt/homebrew/Cellar/spglib/2.8.0/lib
mpifort -ffree-line-length-0 -cpp -mtune=native -DCOLOR -g -O0 -DUSE_MPI \
  -I"$BUILD/modules" tests/validation/fixed_ldau/mpi_ldau_projector.f90 \
  -o "$BUILD/bin/mpi_ldau_projector" -L"$SPGLIB" -Wl,-rpath,"$SPGLIB" \
  "$BUILD/source/librslmto.a" -framework Accelerate -lm -ldl \
  "$SPGLIB/libsymspg.dylib" \
  "$BUILD/source/include_codes/abspinlib/libabspinlib.a" \
  "$BUILD/source/include_codes/strux_lib/libstrux_lib.a" \
  -framework Accelerate -lm -ldl
"$BUILD/bin/mpi_ldau_projector"
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  mpirun -n 2 "$BUILD/bin/mpi_ldau_projector"
```

Both rank counts should report `P12=(0,-0.5)`, `restart_expand_error` below
`2e-16`, and `PASS`. The rank-1 and rank-2 captured logs and the MPI build are
preserved under `/tmp/rslmto-ldau-phase2-mpi-build/harness/`.

The real-space Green-function route was also exercised with the supported
scalar nsp=1, two-step U=2 eV/J=0 deck prepared from the repository's
`tests/scf/cases/bulk/bccFe_hubbard` fixture. The original fixture input had
SHA256 `57ce6bef82ad66ae0449a6a2a5b108a9487a2f62d598e6306a964e6c505bf7bd`;
the original nsp=2 deck was preserved and not run because fixed-U currently
rejects SOC/noncollinear input. The derived deck and provenance are in
`/tmp/rslmto-phase2-smoke/REPORT.md` and
`/tmp/rslmto-phase2-smoke/rs_u2/input.nml`.

Using the same MPI-enabled executable, separate copies of that deck were run
at one and two ranks for exactly two SCF steps:

```sh
mkdir -p "$BUILD/rs_serial" "$BUILD/rs_mpi2"
cp -R /tmp/rslmto-phase2-smoke/rs_u2/. "$BUILD/rs_serial/"
cp -R /tmp/rslmto-phase2-smoke/rs_u2/. "$BUILD/rs_mpi2/"
(cd "$BUILD/rs_serial" && RSLMTO_OMP_THREADS_SERIAL=1 \
  /bin/bash /Users/andersb/Jobb/rslmto_fable_v4b/tests/run_binary.sh "$BUILD/bin/rslmto.x")
(cd "$BUILD/rs_mpi2" && OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  RSLMTO_OMP_THREADS_SERIAL=1 mpirun -n 2 \
  /bin/bash /Users/andersb/Jobb/rslmto_fable_v4b/tests/run_binary.sh "$BUILD/bin/rslmto.x")
```

The exact run directories are `/tmp/rslmto-ldau-phase2-mpi-build/rs_serial`
and `rs_mpi2`. Both returned success. Their `Fe_out.nml` files are byte
identical (SHA256
`63ec65b8f32c84160036c34f59f8fdb1a9aca6139cd1fe5879a76956d44e4017`), including
all 150 serialized imaginary LDM entries; the largest magnitude is
`1.838751264926971e-8`. `HUBBARD_LDM_CHECK` reports zero anti-Hermitian
residual on both ranks, and the logged Hubbard band-comparison traces agree
at both iterations. These are finite two-step smoke runs, not converged SCF
results. Full rank-specific logs and output namelists remain in the two run
directories.
