# RA2 result: potential moment-mode provenance

**Measured: pass.** A KS-local fixed-U run wrote both `hubbard_radial_moments='local'`
and `hubbard_backend='ks'` to its U-carrying potential. The potential restarted in
an RS run with local mode and completed. Reusing that same host in an RS run
whose mode resolved to spectral produced the actionable fatal mismatch. Removing
both new fields from the U-carrying restart produced the legacy warning and the
spectral run completed. A local-mode restart with its mode retained but the
informational backend field omitted also completed. A no-U output contained
neither new field.

Preregistered decisions and error rationale are in
[`RA_PREREGISTRATION.md`](RA_PREREGISTRATION.md). The driver is
[`rs_embedding_ra2.py`](../../../../tests/validation/fixed_ldau/rs_embedding_ra2.py).

Command:

```text
python3 tests/validation/fixed_ldau/rs_embedding_ra2.py --binary /tmp/rslmto-ra-current-build/bin/rslmto.x --scratch-root /tmp/rslmto-ra2-backend-verified
```

Final-source validation used source commit `8b8f7ef7bfc026fd34a56b37a3ad3281fd0a6832`
and executable SHA256
`6ec0e5820c711115570577ddf3b2ffe9fe8c9e222e6a50061fadf773cffc0041` in
`/tmp/rslmto-ra-current-build/bin/rslmto.x`.
Inputs, potential outputs, and run-log hashes for each scenario are recorded in
[`RA2_PROVENANCE.json`](RA2_PROVENANCE.json); raw files remain in
`/tmp/rslmto-ra2-backend-verified/`.

No scientific tolerance is involved: the mode mismatch and missing legacy field
have categorical fatal/warning outcomes. Exact no-U output identity is measured
separately against a clean executable built from the pinned commit and passed
for the one-step RS and KS bcc-Fe output sets recorded in
[`RA_NO_U_IDENTITY.md`](RA_NO_U_IDENTITY.md).
