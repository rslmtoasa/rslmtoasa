# RS-LMTO-ASA embedding: RA preregistration

These registrations precede RA execution. The source baseline is commit
`8b8f7ef7bfc026fd34a56b37a3ad3281fd0a6832`. Raw execution products belong in
`/tmp`; compact results and provenance belong beside this file.

## RA1: impurity diagnostic identity

Prediction: for a one-step scalar `calctype='I'` fixed-U run, the diagnostic
site type must be the type of recursion atom `ia`, obtained through
`lattice%ib(ia)` while the density comes from `symbolic_atom(nbulk+ia)`. The
diagnostic `E_int`, `E_DC`, and `E_Hub` must match an independent evaluation of
the printed recursion-atom occupation matrix and the same input U/J.

Expected difference: floating-point roundoff only (target `<= 1e-12 Ry` per
energy component), since the checker will call the same production energy
functional on values parsed from the emitted namelist. A larger difference
would indicate wrong atom/type mapping or a serialization mismatch and blocks
RA1.

## RA2: potential-file mode provenance

Prediction: nonzero-U/J output stores the resolved radial mode and backend;
fixed-U reads reject a known mode mismatch, accept a matching mode, and treat
missing provenance as legacy spectral with a warning. No-U output retains its
legacy byte representation because the two new fields are omitted unless a
potential has nonzero U/J.

Expected separation: mismatch and legacy checks are categorical, so expected
failure/warning text is the decision observable. The no-U protection is exact
SHA256 equality against a baseline run with identical build and input. This
experiment gates acceptance of embedded host potentials in the selected
moment convention; it does not test physical correctness of either mode.

## RA3: fixture input guards

Prediction: `nsp=1` fixed-U copies pass the existing supported scalar guard;
the retained `nsp=2` decks fail with the documented fixed-U scalar-collinear
message. This protects the supported fixture path and the unsupported spin
layout boundary. No numerical physics tolerance is needed; process status and
the named guard message decide the result.

## RA4: spectral impurity baseline

Prediction: in the current spectral mode, an RS impurity embedded in its
matching RS bulk host reproduces the host's d occupations, sphere moment, and
LDOS within the corresponding U=0 impurity-minus-bulk discrepancy plus the
measured repeat scatter. This is expected for the report's same-material,
same-U setup, while the exact embedded-site energy is outside this test.
Separate predictions apply to bcc Fe (`U=3.8 eV, J=1.1 eV`) and hcp Co
(`U=3 eV, J=0.7 eV`). A no-U seed and a +U seed should converge to the same
impurity observables within repeat scatter. U=0 controls must run for both
materials.

Before interpreting U>0 differences, compute each control discrepancy and
repeat scatter from the same observable definitions. If a control discrepancy
is below repeat scatter, report the control as null and use the independent
bulk/impurity occupation-matrix oracle for that observable. Any mismatch above
control plus repeat scatter is a measured failure of the frozen spectral
baseline, not grounds to widen the criterion. If required host/impurity seed
structures or supported inputs are unavailable, record the exact missing
asset and do not claim a material baseline.

The RA1 mapping assumption was amended after an exploratory oracle mismatch;
the original wording and the evidence are preserved in
[`RA1_MAPPING_AMENDMENT.md`](RA1_MAPPING_AMENDMENT.md).
