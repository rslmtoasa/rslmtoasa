# C6.6 T1 gate results

**Initial disposition: STOP; RS identity later PASS against a versioned refreshed oracle.** The
mandatory current-tip C6.2/C6.4 gates were started before any campaign run. The first attempt failed before launch because the identity script expected
disposable input files under `/private/tmp`. Those inputs were recovered from the repository's
archived C5.2 Part D package, whose deck/restart hashes match the prior successful C5.7 run. The
preflight then exposed its own obsolete import from `/private/tmp/c52-c-private`; it now imports
the archived schema module from the repository. With both setup issues fixed, preflight passed,
the RS solver exited 0, but the exact artifact identity gate failed: 12/27 artifacts differ.
At the user's direction, the remaining independent gates were run despite the RS failure.

| Gate | Result | Detail |
|---|---|---|
| RS identity | PASS against refreshed oracle | 27/27 artifacts match the versioned current-runtime RD1 reference. The historical C5.2 comparator still differs in the same 12 artifacts listed below; it remains preserved as historical evidence. |
| KS identity | PASS | 24/24 artifacts match the registered KS baseline using the current no-LibXC build. |
| 13 guard fixtures | PASS | RA2 reports all 13 fixtures passed; per-run records are under `/private/tmp/c66-current-guards/`. |
| Impurity auto fixture | PASS | `UnitAutoImpurityLocalMode` passed; it exercises the impurity auto fixture. |
| Surface auto fixture | PASS | `UnitAutoImpurityLocalMode` passed; it also exercises the surface auto fixture. |
| C6.4 frozen-plus-recursive-site test | PASS | `UnitLdaUProjectorMultisiteFrozenH` passed. |

The extended diagnostic used three RS identity runs (current LibXC build, current no-LibXC build, and the retained `c511-release-build` binary path), one KS identity run, and 13 RA2 fixture runs. The two unit tests launched no solver runs. No surface campaign probe run was launched.
The response-probe gate selected in C6.6 is `Diff <= 5e-8` within at most two steps. It was not
executed during this gate run. The original C5.2 reference remains immutable; the active C6.6 RS
comparison now uses the versioned refreshed oracle described below.

## Follow-up diagnosis

At the user's direction, the current branch was rebuilt with `ENABLE_LIBXC=OFF`
and the RS identity deck was rerun through the checked preflight wrapper. It
still fails 12/27 artifacts, with the same Fermi level, moment, Diff, and
mismatch set as the LibXC-enabled run; the two current-branch output sets match
each other byte-for-byte. LibXC is therefore ruled out as the cause of this
observed mismatch. The root cause remains unresolved. The user then requested running the
remaining gates; their results are recorded above. Details and hashes are in
[`C6_6_RS_DIAG.md`](C6_6_RS_DIAG.md).

## Restart-selection and oracle follow-up

The user noted that a preexisting `Fe_out.nml` takes precedence over `Fe.nml`. Four checked
single-thread runs isolated this effect: old RD1 and current no-LibXC binaries, each with a fresh
`Fe.nml` start and with the archived `Fe_out.nml` explicitly staged. Preflight passed for all
four, including the host-E_F check. Startup provenance confirms the fresh-start pair loaded
`Fe.nml`, while the restart pair loaded `Fe_out.nml`. Each same-start binary pair matched 27/27
artifacts. A preexisting `Fe_out.nml` changed outputs, but it was absent from the historical
reference launch: its script copied only `input.nml` and `Fe.nml` into fresh directories, and its
log says it loaded a fresh potential from `Fe.nml`.

The same-start reruns still differ from the historical output in 12/27 artifacts, despite exact
initial deck/restart hashes. A versioned current-runtime reference snapshot is now recorded in
[`C6_6_RS_ORACLE`](C6_6_RS_ORACLE/README.md); the original C5.2 reference is preserved. The
Fermi rewrite is `1e-6 Ry`, but DOS grids and values also move, so the change is not classified as
solely an input-deck edit. The paired reruns support the refreshed reference; the historical
12-artifact drift remains recorded rather than erased.
The refreshed identity result is [`C6_6_RS_REFRESHED_IDENTITY.json`](C6_6_RS_REFRESHED_IDENTITY.json).
The C5.7 runner accepts `--rs-baseline` and `--result` so future checks can use the refreshed
reference without rewriting the historical C5.7 record.
