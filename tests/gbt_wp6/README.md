# GBT WP6a HOH/overlap probes

`wp6_hoh_k4` and `wp6_hoh_rs12` are finite-q production smoke tests for the
same primitive-linked `ee` plus derived `eeo=ee*obarm` operator. The proxy and
Kanpur cases are negative probes: both must stop with the explicit incomplete
formal-overlap diagnostic before an eigensolve.

The independent dense operator proof is `tests/unit/test_gbt_wp6_hoh.py`.
