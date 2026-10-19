# C1 Stage RA4 results

**Status: RA4 complete as a RECORD; all six matched impurity states are C1 RS-accepted.** The table reports saved-potential central-site values (`ql(1,:,spin)`, `l=0..2`), with sphere moment `|Σ_l ql↑ − Σ_l ql↓|`. Ratios are absolute discrepancy magnitudes, `|Δ(U>0)|/|Δ(U=0)|`. C1 thresholds classify state status; RA4 values themselves are records.

| Case | Impurity (d↑, d↓, |M| μB) | Bulk (d↑, d↓, |M| μB) | Impurity − bulk (d↑, d↓, |M| μB) | Matching U=0 difference | Ratio to U=0 (d↑, d↓, |M|) | Max d-LDOS difference (↑, ↓, total; states/Ry/site) | Status |
|---|---:|---:|---:|---:|---:|---:|---|
| Fe U=0, no-U seed | (4.3492017571, 2.1499334725, 2.1038373286) | (4.3491996549, 2.1499410532, 2.1038282828) | (+2.10218e−6, −7.58075e−6, +9.04578e−6) | control | — | (4.19950e−4, 6.59995e−4, 6.59995e−4) | ACCEPTED |
| Fe U=3.8/J=1.1, no-U seed | (4.5541590116, 2.1026148367, 2.3834221902) | (4.5576066748, 2.0990288901, 2.3900696691) | (−3.44766e−3, +3.58595e−3, −6.64748e−3) | (+2.10218e−6, −7.58075e−6, +9.04578e−6) | (1640.046, 473.033, 734.871) | (0.118152, 0.529063, 0.528684) | ACCEPTED |
| Fe U=3.8/J=1.1, +U seed | (4.5541590334, 2.1026148963, 2.3834221455) | (4.5576066748, 2.0990288901, 2.3900696691) | (−3.44764e−3, +3.58601e−3, −6.64752e−3) | (+2.10218e−6, −7.58075e−6, +9.04578e−6) | (1640.035, 473.041, 734.876) | (0.118122, 0.529013, 0.528633) | ACCEPTED |
| fcc Ni U=0, no-U seed | (4.7015212620, 4.7015214254, 2.6115967788e−7) | (4.7015195582, 4.7015196124, 7.9129303110e−8) | (+1.70379e−6, +1.81298e−6, +1.82030e−7) | control | — | (3.39887e−4, 2.89931e−4, 6.09931e−4) | ACCEPTED |
| fcc Ni U=3/J=0.7, no-U seed | (4.7710985382, 4.7710985680, 5.6703342466e−8) | (4.7710981197, 4.7710981488, 2.5491526223e−8) | (+4.18514e−7, +4.19153e−7, +3.12118e−8) | (+1.70379e−6, +1.81298e−6, +1.82030e−7) | (0.245638, 0.231196, 0.171465) | (7.60000e−4, 7.39959e−4, 1.49995e−3) | ACCEPTED |
| fcc Ni U=3/J=0.7, +U seed | (4.7710985264, 4.7710985409, 4.5142321881e−8) | (4.7710981197, 4.7710981488, 2.5491526223e−8) | (+4.06629e−7, +3.92086e−7, +1.96508e−8) | (+1.70379e−6, +1.81298e−6, +1.82030e−7) | (0.238662, 0.216266, 0.107953) | (7.60000e−4, 6.69959e−4, 1.42996e−3) | ACCEPTED |

The d-LDOS maxima compare impurity and bulk on the same fixed `E−EF` grid separately for each material/U setting. Fe U=0 uses `[−0.93075, 0.57489]` Ry; Fe U>0 uses `[−0.83656, 0.66913]` Ry; Ni uses `[−0.93886, 0.56631]` Ry. Each grid has 2510 points. The “total” column is the maximum absolute difference of the spin-summed d-LDOS, so it need not equal the sum of the separate spin maxima.

## Interpretation

1. All six RA4 RECORD states meet C1 acceptance using the complete last-20 P_OUT history; per-run hashes, paths, commands, and thresholds are in the JSON sidecar.
2. Ni U>0 discrepancies are smaller than their matched U=0 control in all three observables for both seeds.
3. Fe U>0 discrepancies remain much larger than the U=0 control; the C2 D3 requested ratio below 10 is **not met** (absolute moment ratio ≈734.9).
4. D1 found matched first-iteration Hubbard blocks, values, and insertion counts; D2 isolated an impurity-type indexing inconsistency, which was fixed without changing a physics formula and covered by regression.
5. The frozen Green comparison shows the first differences in computed Green spectra after matching operators; weighted Green integration reproduces P_OUT d occupations. It does not explain the later multistep Fe moment discrepancy, which remains unresolved.
6. The supported hcp-Co setup failure was triaged as pre-existing; C2’s selected substitute, generated-lattice fcc Ni, completed the matched controls and both U>0 seed records.
7. D3’s below-10 target is reported as unmet within this RECORD; it is not a new gate or a C1/C2 stop condition. Continue with the authorized RB stage.

## C1 close (5 lines)

Passed: C1 RA4 acceptance for Fe and fcc Ni records; D1/D2 operator-consistency regression; final mapped Debug no-U identity.
Recorded: six impurity-vs-bulk cases, matched U=0 controls, d-LDOS maxima, and C2 triage/D3 Green localization.
Provisional: none among the six RA4 records; Fe multistep discrepancy remains unresolved.
Unresolved: Fe U>0 discrepancy ratio target below 10; RB0–RB6 gates and records remain.
Next: execute C1 Stage RB0–RB6 under C2, stopping only on an enumerated C1/C2 condition.
