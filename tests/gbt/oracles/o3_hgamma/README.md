# O3 (draft for developer review): exact-ASA H^gamma reference, bcc Fe

Not committed. Not a test. Reference values and two internal checks for one-site bcc Fe on the fixed `tests/spin_response/fe_bcc` potential,
built so that it does not depend on the code's first- or second-order Hamiltonian.

| File | Content |
|---|---|
| `make_hgamma_reference.py` | numpy-only generator; prints every check and writes `o3_fe_bcc.nml` |
| `o3_fe_bcc.nml` | values: checks, 5-k-point eigenvalues, J(0) - J(q) on four meshes, O3b moments on three meshes (`&o3b_moments_mesh{4,8,16}`, appended at the end) |

## Inputs (exported by the code, none of its results is used except in P0)
Run the built binary once per deck from `tests/gbt/gbt_app_common.py` (`make_deck`, `run`), kind `bccFe_prim`, mesh 4x4x4:

```
gbt_q0   : make_deck(.., "bccFe_prim", "gbt_single_q", (4,4,4), q_dir=(0,0,0), order="first", theta_deg=0.0)   -> gbt_bonds.out, Fe.nml, input.nml, lattice.nml
nc_first : make_deck(.., "bccFe_prim", "periodic_nc",  (4,4,4), order="first",  spin_resolved=True)           -> fermi_surface.{bin,meta,spin.bin}   (P0 only)
nc_second: make_deck(.., "bccFe_prim", "periodic_nc",  (4,4,4), order="second", spin_resolved=True)           (P0 only)
python3 make_hgamma_reference.py <export dir containing gbt_q0, nc_first, nc_second> o3_fe_bcc.nml
```
- S^alpha: the 15 raw blocks `raw_orb(i,j)` of `gbt_bonds.out` (onsite, 8 NN, 6 NNN), d Cartesian in Angstrom (full precision). Used as written;
  phase exp(+2 pi i k_direct . d_direct), d_direct from `lattice.nml` and alat. S^alpha(k) is Hermitian to 1.8e-15.
- Parameters per l and spin from the `par` block of `Fe.nml`: `c, enu, srdel, qpar, vmad, ws_r`; `wsm = wav * ang2au` (`input.nml`, `math.f90:92`),
  `wow = wsm/ws_r`; `dele = srdel*wow^(0.5-l)`, `gamma = qpar*wow^(1-2l)`, `C = c + vmad`, `E_nu = enu + vmad` (`symbolic_atom.f90:297-306`);
  alpha = (0.3485, 0.0530, 0.0107) (`lattice.f90:48`, printed in `str.out`). The stored `center_band`, `width_band`, `obar` equal the recomputation to 2.2e-16.
- Orbital index -> l, from S only: S_ij(-d) = (-1)^(l_i + l_j) S_ij(d) over the bond pairs (+d, -d). The index with a unique onsite S value (3.0921)
  is s (even); parity relative to it separates p (odd) from d (even). Result `[0,1,1,1,2,2,2,2,2]` (counts 1, 3, 5). Onsite S diagonal:
  3.0921 (s), 2.7864 x3 (p), 2.7092 x3, 1.2983 x2 (d). Confirmed by P0 below.

## Construction (R1, local frames)
- Collinear: S^gamma = S^alpha [1 - (gamma - alpha) S^alpha]^-1, H^gamma_sigma = C_sigma + dele_sigma S^gamma dele_sigma.
- Spiral (rigid rotation of the ASA potential): S_rot(k) = U^dag diag(S(k - q/2), S(k + q/2)) U, 18x18, spin-blocked, U = [[cos t/2, -sin t/2], [sin t/2, cos t/2]]
  (`gbt_reference_frame`, phi = 0). Local-frame diagonals Gamma = gamma - alpha, C, dele. H_rot = C + dele S_rot [1 - Gamma S_rot]^-1 dele.
- Omega: Mermin grand potential (R11) of the eigenvalues of H_rot at fixed mu; mu solved once for q = 0, theta = 0 on the same mesh, 8 electrons, T = 300 K
  (k_B = 6.3336814e-6 Ry/K as in `reciprocal_mod`). Gamma-centred mesh, full BZ, q/2 on the mesh.
- J(0) - J(q) = Richardson (Neville in sin^2 theta) of [Omega(q,theta) - Omega(0,theta)]/sin^2 theta over theta = 0.25, 0.5, 1 deg.

## Self-checks (measured by `make_hgamma_reference.py`)
| check | what | measured | target |
|---|---|---|---|
| P0 | my H_bar = center_band + width S width (first order) and h - h obar h + E_nu (second order) against the code's spin-resolved `periodic_nc` eigenvalues, 5 k, both spins | 1.0e-15 and 2.9e-15 Ry | 1e-10 Ry |
| C1 (Andersen & Jepsen 1984) | every eigenvalue eps of H^gamma_sigma(k) and of H_rot(k) makes P^alpha(eps) - S^alpha(k) singular; (P^alpha)^-1 = (P^gamma)^-1 + gamma - alpha, P^gamma = (eps - C)/dele^2; max smallest/largest singular value | collinear 9.0e-15, rotated 5.6e-15 | — |
| C2 (R2) | theta = 0: eig(H_rot) = eig(H^gamma_up(k - q/2)) U eig(H^gamma_dn(k + q/2)), both q | 3.1e-15 Ry | — |
| C3 (R8) | Richardson of Omega from H_rot against the Matsubara LKAG sum below | 16^3: 1.3e-5, 6.2e-6; 32^3: 9.6e-6, 7.3e-6 (relative, q1, q2) | — |

C1 and C2 hold for any correct implementation of the construction; they check the code of the script, not the physics. C3 compares two different routes
(eigenvalue sum of Omega versus a Green-function expression) and is the physics identity.

### C3: LKAG in the alpha representation, normalisation and sign
Not fitted. With Omega = -kT sum_n e^{i w_n eta} ln det(z_n - H), z_n = mu + i w_n, and det(z - H^gamma) = det(dele^2) det(1 + P^gamma Gamma) det(P^alpha(z) - S_rot) / det(1 - Gamma S_rot),
the factors other than det(P^alpha - S_rot) are rotation-independent or z-independent. In the lab frame P^alpha - S = Pbar + (Delta/2) sigma.e(theta) - diag(S(k - q/2), S(k + q/2)),
Delta = P^alpha_up - P^alpha_dn. Expanding to second order in theta (the sigma_x term at second order; the sigma_z term is q-independent after the k sum):

    Omega(q,theta) - Omega(0,theta) = sin^2(theta) [J(0) - J(q)],   J(0) - J(q) = (kT/4) sum_n [ <T(q,z_n)> - T_inf(q) - (<T(0,z_n)> - T_inf(0)) ],
    T(q,z) = Tr[ Delta(z) g_up(k - q/2, z) Delta(z) g_dn(k + q/2, z) ],  g = (P^alpha - S^alpha)^-1,

which is the convention H = -sum_{i/=j} J_ij e_i.e_j, E = -J(0) + sin^2(theta) [J(0) - J(q)] of R8 (so J(q) ~ -(kT/4) sum_n <T>). T_inf is T with P^alpha -> (gamma - alpha)^-1
(removes the z-independent part). Two details found while building this check:
- The summand is not 1/omega^2 until |z| >> a few Ry (P^alpha reaches its asymptote only at |z - C| >> Delta/|gamma - alpha|, about 1 to 3 Ry per l channel).
  The generator sums to 10, 20 and 40 Ry and adds the tail from a least-squares fit of sum_p a_p z^-p (p = 1..5, real a_p, window cutoff/3 to cutoff), summed to n = 4e6.
- The symmetric real-part sum omits the convergence-factor term: kT sum_n e^{i w_n eta} a_1/z_n = a_1 f(-mu), the symmetric sum gives a_1 (f(-mu) - 1/2). The generator adds (1/4)(a_1/2) to J(0) - J(q), a_1 being the fitted
  1/z coefficient of the q-dependent summand. Without it the LKAG values are 15.8% and 8.8% above the H_rot values at 16^3 (3.277, 14.327 mRy against 2.831, 13.170); with it they agree to 1e-5.
  This is the xi/2 term of ln(1 + e^{-beta xi}) = ln cosh(beta xi/2) + xi/2 - ln 2. Cutoff convergence (mRy, tail + convergence-factor term, q1 / q2):
  16^3: 10 Ry 2.834936 / 13.178322; 20 Ry 2.830643 / 13.169846; 40 Ry 2.830693 / 13.169970. 32^3: 3.894416 / 11.192946; 3.890126 / 11.184475; 3.890176 / 11.184599.
- 24^3 cannot be used: q1/2 = (0.0625, -0.0625, 0.0625) is 1.5 mesh points. C3 uses 16^3 and 32^3 (q/2 on the mesh needs N a multiple of 16).
- T(z_bar) = conj T(z) holds to 1e-15 (bcc; used to sum n >= 0 only).

## Values (`o3_fe_bcc.nml`), J(0) - J(q) in mRy from H_rot Richardson, 300 K
| mesh | mu [Ry] | q = (0.125,-0.125,0.125) | q = (0.25,-0.25,0.25) |
|---|---|---|---|
| 32^3 | -0.083673 | 3.890214 | 11.184681 |
| 48^3 | -0.083794 | 3.859043 | 10.980890 |
| 64^3 | -0.083823 | 3.830076 | 10.975906 |
| 96^3 | -0.083819 | 3.820189 | 11.009264 |

C3 LKAG values (cutoff 40 Ry with fitted tail): 16^3 2.830693 / 13.169970; 32^3 3.890176 / 11.184599 mRy.
Also in the file: for 5 k-points of the 4x4x4 mesh (indices 0, 7, 19, 36, 53) and both spins, the eigenvalues of H_bar (first order), the code's second-order H and H^gamma,
and for the eigenvalues within 0.1 Ry of mu = -0.069305 Ry the deviations from H^gamma: first order 0.0006 to 0.022 Ry (H_bar - H^gamma, positive), second order -0.0020 to +0.0022 Ry.
Note mu from H^gamma (-0.0838 Ry) differs from the code's first-order mu (-0.0693 Ry) by 0.0145 Ry.

Context only, used in no check: code GBT first order at 96^3 = 3.3514 and 9.6438 mRy; fe_lswt oracle 2.8615 and 7.7902 mRy.

## O3b: l-resolved moments of H^gamma (collinear bcc Fe, per spin)
Same H^gamma_sigma(k) as above (eigenvalues include vmad), Gamma-centred full meshes 4^3, 8^3, 16^3 (weight 1/N^3 per k), 300 K, mu solved by bisection for 8 electrons
(both spins, no spin degeneracy factor). For eigenvectors u_j of H^gamma (Andersen, Jepsen, Sob; Kanpur notes eqs. 2.51-2.55):

    n^q_{l,sigma} = sum_k w_k sum_j f(e_j - mu) e_j^q sum_{m in l} |u_{m,j}|^2,  q = 0, 1, 2;  centre = n1/n0;  central q2 = n2 - 2 centre n1 + centre^2 n0

Orbital index to l is the S-parity assignment above. Written in `&o3b_moments_mesh{4,8,16}` of `o3_fe_bcc.nml`; arrays have shape (l = 0, 1, 2; sigma = 1 majority, 2 minority), Fortran order.

| name | unit |
|---|---|
| `mu_ry` | Ry |
| `n0` | electrons per cell |
| `n1` | Ry electrons (absolute energy scale, vmad included) |
| `n2` | Ry^2 electrons |
| `centre_ry` | Ry |
| `central_q2_ry2` | Ry^2 electrons |

### Self-checks (printed `ok` or `DEVIATION`; threshold 1e-10 set by the generator)
| check | what | measured (3 meshes) |
|---|---|---|
| M1 | sum over l and sigma of n0 = 8 (tests the bisection and the weights) | 5.0e-14 |
| M2 | all moments after H -> V^dag H V, V a random unitary within each l shell, re-diagonalised | 1.8e-15 |
| M3 | 4^3 eigenvalues at the 5 k-points of the table above against the `h_gamma` values of `&o3_kpoints` as read from the previous `o3_fe_bcc.nml` on disk (n/a if none) | 4.0e-15 Ry |
| M4 | eigenvalues and mu shifted by 0.37 Ry: n0 and central q2 unchanged, centre shifted by 0.37 | 6.9e-15 |

M3 compares the same construction (`h_gamma`) with its own earlier output; it checks the mesh bookkeeping (k ordering, eigh against eigvalsh), not the physics. M2 and M4 hold analytically
(unitary invariance of the shell trace; algebra of the central moment); they check the code of the script. None of M1 to M4 is validation of the moments against an independent route.

### Values (mu [Ry]; moments are in the `.nml`)
| mesh | mu [Ry] |
|---|---|
| 4^3 | -0.067451 |
| 8^3 | -0.083497 |
| 16^3 | -0.082936 |

n0 per (l, sigma) at 16^3: s 0.330816 / 0.337962, p 0.378182 / 0.432434, d 4.401480 / 2.119126 (majority / minority). Central q2 at 16^3 [Ry^2 electrons]:
s 0.005294 / 0.004916, p 0.005354 / 0.006414, d 0.046979 / 0.010940. From 4^3 to 16^3 d n0 changes by +0.0647 (majority) and -0.1519 (minority).
Reproducible as for O3: `python3 make_hgamma_reference.py export o3_fe_bcc.nml` from this directory, numpy only, deterministic apart from round-off (existing values changed between the committed file and the rerun with O3b added by at most 1.9e-14 (LKAG, absolute) and 9.2e-11 Ry (J(0) - J(q), 2.4e-8 relative); the O3b block has no random input except the unitaries of M2, seed 2). Wall time of the full run with O3b: 277 s.

## Not covered
Orbital layout of the strux library itself (the code's own S layout is used as exported); more than one site; SOC; any comparison with the real-space
recursion route. Wall time of the full run: 235 s (6 worker processes, Apple silicon): checks 16^3 17 s and 32^3 124 s, values 3, 8, 20, 64 s for 32^3 to 96^3.
