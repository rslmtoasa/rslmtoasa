# VAL-07 — Advanced exchange formulations

## Scope

This validation uses the existing two-pair bcc-Fe exchange fixture. The
fixture is sufficient for the compact algebraic and indexing contracts that
VAL-06 identified; it is not a noncentrosymmetric SOC material benchmark.

The campaign is registered as `Val07ExchangeFormulations` under the
`validation` CTest label. It checks production output from
`calculate_exchange` and `calculate_exchange_twoindex` without reimplementing
the Green-function or exchange integrals.

## Oracles and evidence

| Contract | VAL-06 justification | VAL-07 check |
|---|---|---|
| Full tensor assembly | (T=J I+D^{\mathrm{skew}}+A) is an exact source identity | Reconstruct both printed 3×3 tensors from `jij.out`, `dij.out`, and `aij.out` |
| Two-index scalar exchange | (J_{SO}=CD-SD+CC-SC), (J_{FO}=CD+SD-CC-SC) | Recombine `jijparts.out` and compare with `jijso.out`/`jijfo.out` |
| Two-index DMI | (D_{SO}=2(SC+CC)), (D_{FO}=2(SC-CC)) | Recombine `dijparts.out` and compare component-wise |
| Two-index anisotropic exchange | (A_{SO}=SD+SC), (A_{FO}=-SD+SC) | Recombine `aijparts.out` and compare all nine components |
| Site-index structure | VAL-06 maps pair-local arrays and MPI reduction by pair slot | Check both distinct pair vectors and both rows independently |
| DMI limit | bcc-Fe fixture is scalar-relativistic and has no physical nonzero-DMI term | Require the reported D vector to remain at numerical zero; this does not validate nonzero DMI |

The campaign observes recombination errors at the precision of the printed
parts files (six decimal places). It does not compare SO/FO results with the
full result: VAL-06 identifies those as intentionally inequivalent
perturbative orders.

On the current GNU debug build, the maximum absolute recombination residuals
were (1.05\times10^{-7}) for J, (2.94\times10^{-6}) for D, and
(5.62\times10^{-7}) for A. The largest D vector in the scalar-relativistic
bcc-Fe limit was (2.36\times10^{-6}), below the campaign's (10^{-5})
symmetry-limit envelope.

## Explicitly uncovered formulations

The following remain **Experimental** in the Phase-II ledger:

- Gauss-Legendre exchange: the current PAOFLOW fixture does not share a
  Hamiltonian, Fermi level, contour, and normalization with the Simpson route;
  no equality is asserted.
- Auxiliary-GF exchange: no active caller or matched orthogonal collinear
  representation fixture exists.
- `calculate_jijk`: no active caller; its legitimate oracle is a controlled
  finite-difference displacement derivative of pair J.
- `calculate_exchange_rs2pao`: no native export/import round-trip campaign;
  its A-file I/O contract is also outside this validation.
- `calculate_jij_auxgreen`, `calculate_jijk`, and the `calculate_dij`/
  `calculate_aij` stubs are not promoted by the bcc-Fe result.

No exchange formula was changed to obtain agreement. In particular,
`jtens.out` remains unwritten, as documented by VAL-06, and is not treated as
an output oracle.

## Checklist

- [x] Every tested equality justified by VAL-06
- [x] scalar/tensor relation tested where applicable
- [x] DMI/tensor relation tested where applicable
- [x] anisotropic part tested where applicable
- [x] integration-route equivalence assessed where applicable — no legitimate
  matched Simpson/Gauss comparison exists at current HEAD, so that route
  remains Experimental
- [x] multi-site indexing exercised where required
- [x] no formula changed merely to force agreement
- [x] uncovered formulations remain explicitly Experimental
- [x] maturity ledger updated by formulation

## Files and command

- `tests/validation/val07_exchange_formulations.py`
- `CMakeLists.txt`
- `docs/dev/PHASE_II_VALIDATION.md`
- `docs/dev/VAL-07_ADVANCED_EXCHANGE.md`

Run with:

```text
ctest --test-dir <build> -R Val07ExchangeFormulations --output-on-failure
```
