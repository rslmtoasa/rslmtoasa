# Deterministic local response

In an existing frozen conductivity deck set cond_calctype='per_type',
recur='chebyshev', transport_backend='legacy', linear_in='charge',
linear_out='spin_torque_onsite' and response_projection='local'.
Set pol_alpha='x' and v_beta=1,0,0 for the x spin component and x drive.

The historical default is seed projection. For onsite torque the measured block
is the electronic torque a_i = −i[Sα,Hii(nonSOC)], with S = σ/2 in the global
transport basis. The unscaled native action on a site identity includes hall
overrides and HOH return paths. Assembled Hubbard and constraint fields enter
through that action. Spin-dependent hoppings contribute to the historical full
spin torque. Strict radial XC torque and the reaction on a magnetic order
parameter require separate definitions.

The local calculation propagates X_n = T_n Z_i once per representative.
Propagation from the alternative seed Z_i a_i would produce X_n a_i by
linearity. The stored block is a_i X_n† B X_m, computed as the adjoint of
X_m† (B X_n) a_i. The right multiplication mixes reference columns. It leaves
the velocity acting on the full propagated cluster. Compact reconstruction
uses the complex conjugate of the contracted diagonal. The two-recurrence
formula and the direct dense trace provide independent checks.

In seed mode the onsite output operator acts at every site, even near an
artificial boundary. Its selected spin component costs 16 nb² bytes per
cluster site in FP64 complex storage. Local mode extracts only the currently
traced site's block. Both use one outer trace loop. Seed uses two Chebyshev
families and one input velocity per trace. Local uses one family and M velocity
applications for M polynomial orders. Gamma, energy integration, normalization
and the meaning of existing outputs are inherited from conductivity.

Physical site contributions add with the existing cell-area or cell-volume
denominator. The raw total averages the trace seeds. The feature-mode CSV
identifies every type and representative independent of repeated chemical
labels. Sum physical contributions over the atoms in a layer. Seed and local
totals agree for the same onsite operator on a complete finite trace or
translation-equivalent cell. In an open cluster, test the lateral extent at
fixed spectral window and order while keeping the physical surfaces.

Current implementation requires reciprocal directed Hamiltonian blocks and
Hermitian onsite terms. Overlap and Eν are checked for HOH. A failed check
reports the sites and residual. Legacy FP64 moment generation supports the
feature. A fast_sp request uses the logged CPU fallback. Reciprocal moment
generation, CCOR, imported effective matrices, local spin frames, textures
and spdf formatted transport output are guarded. See the validation record
for measured coverage and outstanding limitations.
