# RA4 raw-log manifest and excerpts

Full raw run trees and logs remain in `/tmp/rslmto-ra4-frozen/` as required by the protocol. This compact manifest preserves original log hashes, paths, and only the decision-relevant terminal SCF marker; no full raw logs are copied into the repository. The detailed extracted observables and their source/output hashes are in `../RA4_RESULTS.md`.

| Run | Raw log path | Raw log SHA256 | Terminal excerpt |
|---|---|---|---|
| Fe U=0 bulk repeat1 | `/tmp/rslmto-ra4-frozen/fe_bulk_corrected/u0/repeat1/run.log` | `fe3abf4bd00af0532b628f2bc345284cfcc7767f02688bff123952c88d18aad1` | `[Oct 18 2026 21:28:42.393] [self.f90:4011]             [1m[32m<info> [0m[0mConverged!0.0000000046` |
| Fe U=0 bulk repeat2 | `/tmp/rslmto-ra4-frozen/fe_bulk_corrected/u0/repeat2/run.log` | `bc398d67dcef5f83856cb4a430ff193ea0437bbe4aed318b3484204d40f23c77` | `[Oct 18 2026 21:37:45.818] [self.f90:4011]             [1m[32m<info> [0m[0mConverged!0.0000000046` |
| Fe U=0 impurity repeat1 | `/tmp/rslmto-ra4-frozen/fe_impurity/u0/repeat1/run.log` | `50ffb7236ad959d31ae3806b716c34b99c6a251af8dfe08661df283df337e982` | `[Oct 18 2026 22:00:11.314] [self.f90:4011]             [1m[32m<info> [0m[0mConverged!0.0000000043` |
| Fe U+ bulk initial 100-step run | `/tmp/rslmto-ra4-frozen/fe_bulk_corrected/u38j11/repeat1/run.log` | `8111030bfeeb21521a9801bc6082c3843a07a98cd26f3113971d402513417317` | `[Oct 18 2026 21:40:40.823] [self.f90:4019]             [1m[32m<info> [0m[0mNot converged! Diff= 0.0000000205` |
| Fe U+ beta .01 continuation stopped at step 105 | `/tmp/rslmto-ra4-frozen/fe_bulk_corrected/u38j11/repeat1_extend/run.log` | `7bb89e41fa01102526d4975339819a79acdb5f929b1fed609eafb828d51789fd` | `[Oct 18 2026 22:03:27.390] [self.f90:4019]             [1m[32m<info> [0m[0mNot converged! Diff= 0.0000000144` |
| Co U=0 bulk initial 100-step run | `/tmp/rslmto-ra4-frozen/co_bulk/u0/repeat1_full/run.log` | `47eb9d2746bed444f53faa78e1956fb351fea5d30164c99e28478d3f89936974` | `[Oct 18 2026 22:06:05.323] [self.f90:4019]             [1m[32m<info> [0m[0mNot converged! Diff= 0.0000000148` |
| Co U+ bulk initial 100-step run | `/tmp/rslmto-ra4-frozen/co_bulk/u3j07/repeat1_full/run.log` | `422d36f49bf5051a807a726dcdbf68f414bb439563dd4e47ba8a5c28c335d0c4` | `[Oct 18 2026 22:07:40.913] [self.f90:4019]             [1m[32m<info> [0m[0mNot converged! Diff= 0.0000000063` |

| Fe U=0 impurity matched one-step feasibility | `/tmp/rslmto-ra4-frozen/fe_u0_impurity_matched/feasibility/run.log` | `d64755f2a2be6a794f17d3f0c6e51f1881983d75fa4dad61eadfbead33e1cb28` | `[Oct 18 2026 22:36:08.671] Not converged! Diff=0.0798779711` (expected nstep=1 cap; normal exit) |
| Fe U=0 impurity mismatched repeat2 interrupted after step 39 | `/tmp/rslmto-ra4-frozen/fe_impurity/u0/repeat2_threadcap/run.log` | `28685b2c3864dd2332ca92305024a8e51e275bc9fffc7ad953f38214c112f628` | `Perform recursion at step 40` (partial step; excluded because bounds mismatch) |
| Fe U=0 matched impurity repeat1 | `/tmp/rslmto-ra4-frozen/fe_u0_impurity_matched/repeat1/run.log` | `d94308d8bdd21f7c7b5966fcc36aad4c39141ef25b8fb9c4bbb3fe0ec8ce7f34` | `[Oct 18 2026 23:02:40.658] Converged!0.0000000046` at iteration 86 |
| Fe U=0 matched impurity repeat2 | `/tmp/rslmto-ra4-frozen/fe_u0_impurity_matched/repeat2/run.log` | `087b7202cde54ceb6189808b83c9e1377c7b3ab1f868c4bfe0e75e8acde3bb29` | `[Oct 18 2026 23:02:32.493] Converged!0.0000000046` at iteration 86 |
| Fe U+ beta=0.03 continuation at nstep=300 | `/tmp/rslmto-ra4-frozen/fe_bulk_corrected/u38j11/repeat1_beta003/run.log` | `989055468d5f77e165910b03d8922d55a306304fa5f54a20be83376e4f4a6d91` | `[Oct 18 2026 23:07:26.172] Not converged! Diff=0.0000000201`; same-step fixed-U `rp_abs=1.0982e-7`, `rv_abs=2.0164e-8 Ry`, streak 0 |

| Co U=0 bulk extension, 300-step cap | `/tmp/rslmto-ra4-frozen/co_bulk/u0/repeat1_extend300/run.log` | `0af9f23741635af71978391dc6d2025b748a814aad9ec1167eee4fb2267824d1` | 300 iterations; `Not converged! Diff=0.0000000179`; normal timer report |
| Co U=3/J=0.7 bulk extension, 300-step cap | `/tmp/rslmto-ra4-frozen/co_bulk/u3j07/repeat1_extend300/run.log` | `5b87ddea01dc44ed2d6a27a09370a000b2fd2efe5abc34344784f649e3289490` | 300 iterations; `Not converged! Diff=0.0000000114`; same-step fixed-U `rp_abs=1.2498e-7`, `rv_abs=2.4256e-8 Ry`, streak 0; normal timer report |
