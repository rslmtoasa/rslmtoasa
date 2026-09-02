# TDDFT-12 AF/ferri chirality validation

## Result

The collinear SOC-free transverse path now treats the ordered circular
responses independently:

\[
\chi^{+-}: (O_+,O_-), \qquad \chi^{-+}: (O_-,O_+).
\]

The new `&tddft` key `circular_channel` accepts `plus_minus`, `minus_plus`,
or `both`; `both` is the default.  A `both` run retains the historical
primary filenames for `+-` and emits separate `_minus_plus_...` files for the
reverse channel.  Every chi0/Dyson/Goldstone record carries the ordered
channel in its metadata.

## Signed multi-sublattice handling

The response-projector population remains a positive magnitude for legacy
radial ALSDA normalization.  The signed `P_site sigma_z` population is stored
separately and is used explicitly for the transverse Ward/Goldstone vector.
The LMTO pair-potential source selects independently assembled `Q^-` or
reverse `Q^+` tiles; the reverse route is not obtained by reusing the
historical single-channel output.

## Reproducible evidence

`tests/unit/test_tddft_af_ferri_chirality.f90` is a two-sublattice analytic
AF/ferri fixture with one up-polarized and one down-polarized site.  It checks:

- the positive-energy AF poles occur in opposite `+-` and `-+` channels;
- the two ordered channels remain numerically distinct;
- the eigenpair and K-space Green-function chi0 backends agree for both
  channels;
- a signed `[+2,-1]` multi-sublattice Goldstone vector satisfies the static
  Ward identity and is retained in the diagnostic record.

Commands used:

```text
cmake --build build --target UnitTddftAfFerriChirality -j2
ctest --test-dir build -R UnitTddftAfFerriChirality --output-on-failure
```

The named FeRh and CrMnSb production campaigns remain optional higher-cost
material cases: this checkout contains no corresponding converged material
input/evidence set.  The compact fixture is therefore the committed
regression gate for the implementation-level TDDFT-12 behavior.
