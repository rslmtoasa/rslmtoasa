# WP06 nonzero constraining-field fixture

`UnitGbtConstraintCovariance` is the deterministic WP06 regression fixture.
It supplies a fixed target `e_target = (0,0,1)` and a sequence of canted,
unit-length measured moments.  The zero-field snapshots have a measurable
angular/transverse residual; eleven prescribed updates reduce it below `1e-6`
rad while each controller update produces a nonzero onsite field and reports a transverse field for the
orthogonalized controller (`constraints_i_cons = 3`).

The same local target/moment/field triplet is transformed through the
period-four `gbt_frame_for_cell` contract for four explicit supercell sites.
The test requires equal field magnitude, equal target/actual angle, and
cell-by-cell equality with the rotated primitive field.
