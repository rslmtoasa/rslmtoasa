# S0b results

Base `fable_v4b_u` at `876063c19a4f2821dc55a25127244f0838b2bb56`; work branch `fix/ldau-s0b`. The approved scope is preserved verbatim in [`S0B_AMENDMENT.md`](S0B_AMENDMENT.md). All solver launches used `.venv/bin/python tools/run_rslmto_checked.py` and the checked single-thread launcher.

## T0: S0 checks

The ten-short-run allowance was fully used. Run records and outputs are in `/private/tmp/rslmto_s0b_t0`.

| Run | Case | Outcome |
|---:|---|---|
| 1 | RS no-U identity case | Completed; identity comparison uncovered |
| 2 | KS no-U identity case | Completed; identity comparison uncovered |
| 3 | RS +U identity case | Completed; identity comparison uncovered |
| 4 | KS +U identity case | Completed; identity comparison uncovered |
| 5 | +U impurity with the default host | Completed; identity comparison uncovered |
| 6 | +U impurity with a fresh RS +U density-matrix host | Completed; identity comparison uncovered |
| 7 | Referenced density-matrix file missing | Fatal as expected: `Referenced density-matrix file is missing: ./Fe_dm.nml` |
| 8 | Initial SHA mismatch probe | Invalid fixture: the edit did not change the sidecar; not scored |
| 9 | Corrected SHA mismatch probe | Fatal as expected: `Density-matrix SHA-256 mismatch for ./Fe_dm.nml` |
| 10 | Legacy inline `ldm_flatten` host | Completed; input had inline LDM and no `dm_file`/`dm_sha256` |

The identity gate is **uncovered**, not passed or failed. The copied S0 output directories contained solver-mutated decks. For example, the fresh RS no-U deck hash was `2a33bb6c…`, versus the saved S0 preflight hash `0ee4025d…`; the Fermi value was `-0.069377` rather than the pre-run `-0.069378`. That changed the scientific outputs, so the fresh outputs are not comparable to the baseline. The +U sidecar hash did match its potential's `dm_sha256` (`d367ef68…`), but old/new +U values were not compared. The request says “all 6 decks” while enumerating five; run 6 was an additional DM-backed-host impurity case. The 13 RA2 guard fixtures each launch the solver and remain uncovered under the ten-run limit.

## T1: magnetic Ni host

After the first T1 pair, the user supplied bulk Ni references and directed that they be used. Their sphere moments from `Ni_out.nml` are:

| Supplied reference | Stored U/J | Moment | Use |
|---|---:|---:|---|
| `/Users/andersb/Jobb/RS-utils/tmp/Cheb/fccNi/Ni_out.nml` | 0/0 | `0.6707045323 μB` | Magnetic no-U cross-check |
| `/Users/andersb/Jobb/RS-utils/tmp/Cheb/fccNi/zeroU/Ni_out.nml` | 0/0 | `0.6707047347 μB` | T1(a) seed |
| `/Users/andersb/Jobb/RS-utils/tmp/Cheb/fccNi/U/Ni_out.nml` | 4/1 eV | `0.8086561055 μB` | Independent finite-U magnetic reference; T1(b) must restart from T1(a) |

The user references were generated with Chebyshev recursion at `rc=60`; the checked T1 runs used the requested RS fcc block recursion at `rc=80`, `nsp=1`, and fixed Fermi energy `-0.061141 Ry`.

| State | Setup and C1 acceptance | Endpoint moment | Result |
|---|---|---:|---|
| T1(a), U=0 | Seeded from the explicit-zero-U reference. Earliest registered C1 window: iterations 53–72; Diff `1.15e-8`, spin-q0 spreads `3.7788e-7` / `4.6913e-7`, moment spread `1.6313e-7` | P_OUT q0-derived moment `7.5853e-8 μB`; saved post-mix state ql moment `1.6378e-9 μB` | RS-accepted; collapses below the `>0.4 μB` magnetic gate |
| T1(b), U=3/J=0.7 | Restarted from the accepted T1(a) state; `auto` resolved to `local` (`radial_source=default`, clean-bulk eligibility reason). Earliest registered C1 window: iterations 19–38; Diff `1.87e-8`, spin-q0 spreads `4.9098e-7` / `4.5887e-7`, moment spread `1.3731e-7` | P_OUT q0-derived moment `1.3419e-8 μB`; saved post-mix state ql moment `4.8535e-9 μB` | RS-accepted; non-magnetic collapse below the `>0.4 μB` gate |

The accepted T1(b) endpoint is the specified stop condition. T2 and T3 were therefore not run; their ratios, controls, predicted spin-swap separation, and moment/occupation differences remain unevaluated. No workaround or production change was made. Corrected T1 artifacts and scans are under `/private/tmp/s0b_ni_t1a_userseed` and `/private/tmp/s0b_ni_t1b_userseed`; scan SHA-256 values are `917b676a1234863a899d518fae7dedaff5be03131f061f9979c51ca6ed745abf` and `06b754f1b11622387d89b990539f59796dfc9b7c6fc447561e294297279e02ce`.

Run accounting: five T1 solver starts total. One initial deck-preparation launch was discarded and corrected. The first T1 pair used the earlier Ni(001) seed; the supplied magnetic bulk references then superseded that pair. The corrected user-seed T1(a)/(b) pair used the additional start explicitly authorized by the user. The formal results above refer to the corrected pair.

## T2/T3 and merge gate

T2 RD2 and T3 RD3(a) were not run because corrected T1(b) collapsed to a non-magnetic state, which is the explicit stop condition. Their ratios, controls, predicted spin-swap separation, and moment/occupation differences remain unevaluated.

The 13 guard fixtures, a scientifically comparable T0 identity result, and the post-merge S0 identity gate remain uncovered. The ten T0 solver starts were exhausted, so no additional short identity/fixture runs were launched. The final merge therefore carries no fresh post-merge identity verdict.
