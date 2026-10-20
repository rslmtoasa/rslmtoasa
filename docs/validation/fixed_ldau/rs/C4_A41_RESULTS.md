# C4 A4.1 — Fe embedding comparison and §2.4 triage

The matched Fe impurity +U run has no accepted endpoint, so no bulk/impurity energy ratio is reported. The preregistered §2.4 result is **PRE-EXISTING ISSUE**: the same accepted-seed deck is non-finite on clean `8b8f7ef`, and the same-deck U=J=0 control is also non-finite. A4.2 independently found coefficient drift but no evaluator inconsistency; the A4.3 decision is therefore **NEITHER**, closing the proposed Fe embedding known issue without a production formula change.

| Run | Outcome | Evidence |
|---|---|---|
| Production U=3.8/J=1.1, rc=100, explicit spectral | Fatal non-finite at attempted step 103 after 102 complete outer iterations; pre-LU Q and mixed qia become NaN, then invalid PL. No accepted C1 endpoint. | `/private/tmp/c4-a41-fe-u38j11-rc100` |
| Same accepted-seed U+ deck, Debug bounds/FPE rerun | SIGILL after 27 complete iterations. Crash report symbolicates to Accelerate `zheev`, called at `source/haydock_fast.f90:72`; this does not establish a bounds violation. | `/private/tmp/c4-a41-triage-debug` |
| Same accepted-seed U+ deck, clean `8b8f7ef` | After stripping only two unsupported potential-output metadata assignments, the restart fields and explicit spectral input were preserved. Fatal non-finite at attempted step 103 after 102 complete rows, before a usable endpoint. | `/private/tmp/c4-a41-triage-clean8b-compatible` |
| Same deck and accepted U+ seed, U=J=0 triage | Fatal non-finite at attempted step 35 after 34 complete rows; no converged U=0 endpoint. | `/private/tmp/c4-a41-triage-u0-from-uplus-seed` |

The clean-8b direct invocation with uncaptured output and its captured parser-incompatibility retry are retained as separate attempts. The retry stopped before loading the seed (potential namelist `iostat=5010`); it is setup evidence only. The subsequent compatible clean-8b run removed exactly the unsupported `hubbard_radial_moments` and `hubbard_backend` assignments from `Fe.nml` and `Fe_out.nml`; `input.nml` remained byte-identical and retained `hubbard_radial_moments='spectral'`. The accepted restart's physical fields were independently verified before launch.

The matched U=0 bulk reference is still completing its prescribed beta continuation. Its current stage is recorded as in progress in the provenance JSON; it does not change the A4.1 result because there is no accepted matched U+ impurity endpoint from which to form the requested ratio. The stage must finish under the registered iteration cap/continuation rules before Part A is closed.

All raw attempt directories, manifests, logs, traces, accepted seed, static inputs, and executable hashes are enumerated with SHA256 in [C4_A41_PROVENANCE.json](C4_A41_PROVENANCE.json). No value is presented as an accepted A4.1 energy ratio.
