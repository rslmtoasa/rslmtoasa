# C6.6 RS identity diagnostic: LibXC build setting

After the current-tip RS identity gate failed, the user authorized continuing the
investigation with LibXC disabled. The current branch was rebuilt in
`/private/tmp/c66-nolibxc-build` with `ENABLE_LIBXC=OFF`, GNU Fortran 16.2.0,
Release mode, and the same recorded Fortran flags as the earlier C5.7 build. The
resulting executable SHA-256 is
`8b4273aaae07df7e02b93a98664da0dd0449ae8684fe8a1034e838911ec381b2`; `otool -L`
confirms it has no LibXC library dependencies.

The archived RS deck and restart were copied unchanged: their SHA-256 values are
`45e28a6c2f2853949d4bf947e9916cd763a71e350d394bbd9e3810b729e0054d` and
`6ded298a3512b920211a8761b05a191a76402d3c985080380034cdff42596673`. The run was
launched through `tools/run_rslmto_checked.py`; its `preflight.json` records
`PREFLIGHT_PASSED` and `RUN_COMPLETED` (return code 0). Artifacts are in
`/private/tmp/c66-nolibxc-rs/`.

The no-LibXC run still fails RS identity: 12/27 artifacts differ from the
registered C5.2 baseline, with the same mismatch set as the LibXC-enabled run.
Both current-branch runs report free Fermi energy `-0.069588`, moment `2.105425`,
and `Diff=0.0001112571`. All 27 candidate output artifacts from the no-LibXC and
LibXC-enabled runs compare byte-identically after the registered `str.out`
`cpu_s` normalization. This rules out the LibXC build setting as the cause of
the observed mismatch on this fixture.

The source diff since the last passing C6.2 gate is limited to the surface and
impurity local-moment paths in `bands.f90` and `hamiltonian_build.f90`; the
identity deck is bulk (`calctype='B'`). No direct cause has been identified.
The RS gate remains failed, so surface campaign runs remain stopped. The
no-LibXC diagnostic adds one solver invocation; no surface probe was launched.

At the user's request, an additional checked RS run used the retained executable
at `/private/tmp/c511-release-build/bin/rslmto.x` (current file SHA-256
`46524061a637a1c9ec7f25df74d26d1d0a74cec65222a2241a2842c990366082`). It also
produced the same 12/27 baseline mismatch and byte-identical output to the current
no-LibXC build. This is the executable currently at that path; its hash is not
bound to the earlier C6.2 passing run. The separate KS comparison passes 24/24;
RA2 passes all 13 fixtures; `UnitAutoImpurityLocalMode` and
`UnitLdaUProjectorMultisiteFrozenH` both pass. No reference oracle was changed.
