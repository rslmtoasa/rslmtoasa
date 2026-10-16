# A5 common U=J=0 baseline continuation, 240-step attempt

This is the full residual log, lightweight trace and final restart from the
uninterrupted continuation started from the archived `Fe_out_input.nml`
(copied from `/private/tmp/a5-baseline-finalrun/Fe_scf.nml`, SHA256 below)
with the archived `input.nml`, executable and physical settings. The later
common campaign restart is a separate mid-continuation snapshot at 00:10:18;
it was copied from this run's `Fe_scf.nml` and is archived in
`../a5_common_seed/`. This baseline segment reached the
240-step cap without a `Converged!` marker; the last residual was `4.8e-9`
under the existing `1e-9` threshold. This is preserved as a non-converged
attempt, not an accepted endpoint.

Input SHA256: `ca65d9695b049f10beb3ee9461c05572b48f3973885de99ea7f8fee3dfc5699c`

Starting seed SHA256: `873c0fb7205319511b16c2c8f088383bc4a18577ba9b33012316f688e9115328`

Executable SHA256: `52eb94c6e6f6ae12fca629a4eea1b30a6800e29736c9fdfdfd836538077f2225`

Final restart `Fe_scf.nml` SHA256:
`9b7ddeace03fab3cf58c56987b65133509cefbe4bee40db7859f99099e108ec2`

Residual log SHA256:
`5c379ede903a6a5792d95f34ce4afa17f3188b5701c5d39f1552f4e7fbbdda55`

Trace SHA256:
`c59718c77b08edaa5b94ae3990308c3a403701da2b8929de332c02a4b9e0595a`
