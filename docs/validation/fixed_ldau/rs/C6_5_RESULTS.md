# C6.5 — Mechanical preflight enforcement

C6.5 adds `tools/run_rslmto_checked.py`. Run it with the repository `.venv` interpreter,
the target run directory, and the intended `rslmto.x` binary. The wrapper invalidates any old
`preflight.json`, invokes `C5_4_PREFLIGHT.py` with the selected deck, guard, restart, and optional
host E_F trace, then atomically writes its decision and input hashes to `RUN_DIR/preflight.json`.
The solver is invoked through the same thread-capped bash guard only after a successful preflight;
any preflight error exits non-zero and records `launch.status=REFUSED`. On success, the final
record binds the pass to the deck, restart, host-E_F trace, guard, and binary hashes and records
the launch result. Campaign auditors must reject a run if the record is missing, failed, or does
not match the run's inputs.

The unit test `tests/validation/fixed_ldau/test_run_rslmto_checked.py` passed all three cases
using `.venv/bin/python3`: an unclosed namelist was refused; a host-E_F mismatch was refused; and
a valid surface deck passed preflight and invoked a no-op `rslmto.x` stub. No SCF solver was run.

The surface stage is PAUSED pending a user-provided converged no-U Ni(001) or Fe(001) reference
with `nsp=1`. No merge was performed.
