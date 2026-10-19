# RA4 retained evidence

This directory preserves the small preregistered input/seed files needed to reproduce the registered RA4 attempts, plus the read-only analyzer and compact extraction snapshots. The corrected matched Fe U=0 impurity inputs/host/seeds are under `prerun_inputs/fe_u0_impurity_matched/`; their exact hashes are recorded in `RA4_PREREGISTRATION_ADDENDUM.md`.

Full run trees remain in `/tmp/rslmto-ra4-frozen/`, as required by the protocol. No full logs are copied into the repository; `logs/README.md` records raw paths, hashes and short decision-relevant excerpts. Paths, executable/input/seed/output hashes, commands and measured endpoints are documented in `RA4_RESULTS.md` and `RA4_PREREGISTRATION_ADDENDUM.md`. The retained files are compact reproducibility inputs and result artifacts, not substitutes for the raw runs.

`ra4_analyze.py` is read-only: it reports extracted values and hashes, labels LDOS units/axis, and does not apply acceptance thresholds. The Fe U=0 comparison JSON records the identical repeat pair, including full-precision potential ql traces and common-grid LDOS differences.

The JSON files are compact, reproducible extraction snapshots: the original Fe U=0 impurity output is exploratory because its energy bounds mismatched the bulk; the corrected matched Fe U=0 pair and bulk comparison are measured and retained separately; the Co cap-100 snapshots are provisional, and the registered cap-300 extensions also failed to converge; their compact incomplete snapshots are retained separately. Their hashes are listed in `RA4_RESULTS.md`.
