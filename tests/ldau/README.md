# LDA+U tests

| Path | Content |
|---|---|
| `test_guards.py` | `LdaUGuards`: every LDA+U guard scenario in one test (nsp=1 host at nsp=2, U-less host, tilted axis with block recursion, mixed axes, host-axis mismatch); prints PASS/FAIL per scenario |
| `reference_rows.py` | stored-reference rows: `reference_rows.py --list`, `reference_rows.py ROW --binary rslmto.x` |
| `references/` | JSON scalars of every row; five rows are registered with ctest, the others are archived data |
| `data/` | fixtures read by the Fortran oracles (`UnitHubbardSpinor`, `UnitLdaULocalMoments`) |
| `oracles/` | the independent numpy/SymPy generators of those fixtures (not run by ctest) |
| `decks/ni001_soc_u3/` | reduced Ni(001) surface deck and frozen hosts used by the surface row |
| `rotate_spin_frame.py` | rotates a host and seeds to a tilted axis (see `docs/validation/fixed_ldau/TESTING_POLICY.md`) |

References are rewritten only with `reference_rows.py ROW --binary BIN --update "reason"`; see `docs/validation/fixed_ldau/TESTING_POLICY.md`
and `tests/README.md` (reference policy). Rows with a `<row>.<platform>.json` use that file on that platform.
