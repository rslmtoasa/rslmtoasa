# A5 non-KS setup failure (2026-10-17)

The first A5 campaign launches did not select the preregistered k-space SCF
backend. Their inputs set `calctype='B'` and a reciprocal mesh, but omitted
`&self use_kspace=.true.`. The executable initializes `use_kspace` to false,
so `calctype='B'` alone selected real-space recursion. Runtime logs confirm
this: they contain recursion-step messages and the Hubbard total reports
`entropy=UNAVAILABLE_KS_ONLY` and `TS=UNAVAILABLE_KS_ONLY`.

These are setup failures, not scientific outcomes. Their observables,
convergence histories and energy comparisons are not used for A5 decisions.
Raw files are retained in place; no historical receipts have been rewritten.
No W2 endpoint or candidate comparison was accepted from these attempts.

The common restart `25c1c826…fec5` came from the real-space U=J=0
continuation, not a KS endpoint. It may be used uniformly as an initial state
for corrected KS cases, and is labelled as an RS-origin seed. The earlier A4
KS seed recorded as `aab3f56d…234df` is absent from the current scratch area;
no replacement is described as that historical seed. Corrected case decks
must explicitly set `use_kspace=.true.` and verify the reciprocal-SCF runtime
marker plus finite entropy and TS before campaign relaunch.

All attempts below used executable SHA256
`52eb94c6e6f6ae12fca629a4eea1b30a6800e29736c9fdfdfd836538077f2225`; the
instrumented source hashes are recorded in the adjacent A5 preregistrations
and case receipts.

| Attempt | Pre-run input SHA256 | Starting seed SHA256 | Retained outputs and status |
|---|---|---|---|
| Main `u1.6_spectral_standard` | `26e60bb2090c239bfe4539bd4cffb211928e51867313ac1821a7a7e8d77089c9` | `25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5` | 240 recursion steps; no positive `Converged!`; stdout `f1985784df6e6d737a8da8a5334c7a441acf213dfd0865a7fc01e75225762691`, trace `b8ce7fc88f99e818b2b45c17f6a68de49ece919cfae15cc26fc1a2518678d67b`, final restart `db03a7bb8a430ad6a626fecfdab6e488724d0250a7bf8017ba4ee3ea34701158`. |
| Main `u1.8_spectral_standard` | `c3fa51631e82030038601ae3129e2919383d63bd016a2b07ef9488574ac27f72` | `25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5` | Interrupted at 32 steps; no positive marker; stdout `8a5c37f4a3298a7f8b9f9f69ab4b7cd3b700c0b5f7fa199fdf74cab13831f4d4`, trace `95d00ac9c6ea2866e29e1db67890b3cf97a80d3c113542c2ace5777e14ce288b`, restart `f6a80572dae0f1057ffb6dece04a7fc736a4a0ac2b75c41bc1c41e6386f3838f`. |
| Main `u1.6_spectral_beta005` fallback | `66acdd307e86c52ba96ebb4d775e69a5eff29a110d230e33acb95b97bd0be5bf` | `db03a7bb8a430ad6a626fecfdab6e488724d0250a7bf8017ba4ee3ea34701158` | Interrupted at 39 steps; no positive marker; stdout `fc58f7541bba8664c685dbda4e1698bc54a48ce43cdb0330546b984ab6f71e4c`, trace `e8e631787b974495c61feab94308b1bf5cd79826a840c59fbbad6a3b55e42ecf`, restart `90920a0ec2d2aa2de4e6fbfbb84f4bf270bae3344688fc68f101001646e78c1c`. |
| Priority `u2_local_standard` | `f0653bef1e3daa7f8206f3c5f8a80b084ec1aed7a3265be0f6236f39abc521d5` | `25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5` | Interrupted at 144 steps; no positive marker; stdout `0e63d5842cb966deebae8396b0b3b5e42074d4fe17111cd7a9d5bff98237657b`, trace `5da0dce5187342b5754bec4feac0341aaaed28a85d2194845f1e7a3d3da52dd1`, restart `d15fe854036825a26d79f5246d274838279714cfbfb3406f8757dbbe0dbb1503`. |
| A5.2 local `lambda=-0.005` | `575c8194ccf744f4b390e9d89facfcfd9690c9f49f862b5d0aaf856d1491f4e0` | `25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5` | Reached 240-step cap; stdout `6a29baac338298e4c18dc924c8102b1e1a922fa141fd1872f596e9a1de27ae56`, trace `08ef54fe528bd248cb3ddc5ea051be6738ea6cf8010bfe6e9d5f26da7c05e8db`, restart `4b85ffa187ff252ec8aa92b8d7ef8e99a6480f47a8c40cc97a1cf15f84c030c2`. |
| A5.2 lambda fallback | `a841dba251b6f0ee4d8f79f8df0098a3c7219812e473b495f42e5b48e1c2c09d` | `4b85ffa187ff252ec8aa92b8d7ef8e99a6480f47a8c40cc97a1cf15f84c030c2` | Interrupted at about 67 steps; streamed stdout was not persisted; trace `25d7979b75870a3517fb9980db5568d65226f3a9ed1f189c53e7a6bb068989f9`, restart `dd424d7cba7fd51b9fb68305240cd11f5364be4cc00ac20d11cb2eefa850e896`. |
| W1 `base_h0` | `b193dada822b6d45a1439acbf15f373b7a6841f42e8d51d59c63282554c8ee93` | `25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5` | Reached 240-step cap; stdout `de59f7d4bbb7af3f207f5af65162b717fffa695e83dfbb91e8b243a6fef0e81f`, trace `7553acab0405c0f9e84a436f44a6642a9a55995b952ede69a069b9db0ba16c52`, restart `acca0bba35b3df651b5f1acc9b25d85c51278b96593d1b62210e138ae70c16c7`. Reviewer setup classification SHA `2228373d36dfedb2e494aee07b7c3f613dc57bd4d46e7ee300f98d25baefdb7d`. |
| W1 beta fallback | `8ac5db3c7288434735ea25be674882e35f8111879268a8eb25cafc9d555eca96` | `acca0bba35b3df651b5f1acc9b25d85c51278b96593d1b62210e138ae70c16c7` | Stopped by reviewer; stdout `74aa3b0235b95fe5441f2a526f0848757d739ecc8ca22888a63398067f4149ac`, trace `3a91c46c29f5905653946468798c29e415a360ff11148a6b06cf9ddde263488f`, restart `aa83781d731cf3d5b5973158df9a44a4dfcc4f5e087cbf9e9b4f1a783f5c5a07`. Reviewer setup classification SHA `320d873f287bbfe3f8b76b19abc5878abba82015505486190f07b049ee02ff20`. |

The archived A5 U=J=0 baseline attempt is also real-space: its tracked input
SHA256 is `ca65d9695b049f10beb3ee9461c05572b48f3973885de99ea7f8fee3dfc5699c`
and it omitted `use_kspace`. Its 240-step history remains useful only as the
provenance of the RS-origin restart, not as an accepted KS baseline.
