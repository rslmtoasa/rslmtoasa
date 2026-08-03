# WP0 finite-q GBT diagnostics

`cases.json` launches the current legacy reciprocal GBT implementation through
the frozen-magnon path. `legacy_finite_q_gbt` is a full-BZ smoke test; the
`*_guard` cases are expected failures. The q=0 reference point is allowed, and
the first nonzero q point must stop before eigensolution when an unsupported
term is present.

Run it manually and expect a nonzero status:

```bash
/Users/andersb/envs/p311/bin/python tests/run_test.py --binary build_13/bin/rslmto.x \
  --cases-json tests/gbt_wp0/cases.json --case-name finite_q_gbt_soc_guard \
  --scratch-root /tmp/rslmto-gbt-wp0
```

The bare-case log must contain `nonzero-q GBT rebuilding the full chemical BZ
mesh`. A guard-case log must contain its named unsupported-combination fatal.
These fixtures have no golden numerical output and must not be promoted to one
before the WP1+ operator migration.
