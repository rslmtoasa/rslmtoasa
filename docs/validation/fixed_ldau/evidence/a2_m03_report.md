# Accepted KS M0.3 frozen-state replay

This is the readable companion to the hash-bound preregistration, machine report, and replay CSV in this directory. It reuses the accepted iteration-32 U=2 eV, J=0 capture; it does not run SCF or alter the locked A2.2 fit.

Capture CSV SHA256: `41e8bec19e50ba9931ccab2077686647300e1142774eea45567d1f65602472bc`  
Manifest SHA256: `4f6f82fa0d6779e07fd6babe21aaa0a70b511c9b8cf65594cf6cb0f17aedc54d`  
Production source SHA256: `d91a468a4493da21710fc7f60e745b5f1bef8f20d5d667ae939c31d5ba3061ef`  
KS executable: `b7ca883ec736211552ca1b8df7e2d519c7806459df9a6238dd597c58f3c49f7f`  
Replay driver SHA256: `0c1c3c12e8501f4642af87168c611d261fb5a34cc4747aaef9c8420c2c6449c4`  
Adapter SHA256: `0743222971e90060ce50e56cbaa5aeb58d6e37b9cdf8f828ea8cb82778ad8a7a`  
A2.2 model SHA256: `0ad6a7703aa26eb132602e7a43773d7619f70651c6b90817e4f945ff259cb993`

Accepted state: iteration 32, charge residual 7.69e-08 (limit 1e-07), final R05 `Rp=4.9022e-08`, `Rv=7.2061e-09 Ry`, streak 3.

## Carrier values

| Spin | q0 | t=Tr(Vin Pout) (Ry) | target shift t/q0 (Ry) | q2 local (Ry²) | SUMEV RAW (Ry) | SUMEV LOCAL (Ry) | ΔSUMEV raw−local (Ry) | signed ratio Δ/t |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 4.52761514798 | -0.271387881185 | -0.0599405807065 | 0.0511128784368 | -1.75501085222 | -1.4852104518 | -0.269800400418 | 0.994150509741 |
| 2 | 2.08410054863 | 0.0167468671703 | 0.00803553704803 | 0.0095449698681 | -0.749738468841 | -0.766483776542 | 0.0167453077002 | 0.999906879888 |

Total `Σt = -0.2546410140151576 Ry`; `SUMEV_RAW−SUMEV_LOCAL = -0.25305509271796467 Ry`; signed total ratio `0.9937719329962356`; absolute carrier error `0.001585921297192927 Ry` ≤ `A=0.001702238347231377 Ry`.

## Five replay cases

| Case | Integrated output l=2 d charge up | down | total | Sphere spin moment | ΔSUMEV from RAW (Ry) |
|---|---:|---:|---:|---:|---:|
| RAW_PL | 4.52762075857 | 2.08410120823 | 6.61172196681 | 2.35630744738 | 0 |
| LOCAL_PL | 4.52762302901 | 2.08410126154 | 6.61172429055 | 2.3563096645 | 0.253055092718 |
| SCALAR_V_PLUS_001 | 4.52762046558 | 2.08410117606 | 6.61172164164 | 2.35630718656 | -0.0661067097743 |
| VOUT_CONTROL | 4.52762302901 | 2.08410126154 | 6.61172429055 | 2.3563096645 | 0.253055083627 |
| P_TRANSPOSE_CONTROL | 4.52762302901 | 2.08410126154 | 6.61172429055 | 2.3563096645 | 0.253055092718 |

The integrated charge values above are from the actual NEWRHO output-density l=2 contribution using the ATOMSC logarithmic Simpson metric. Frozen-input-density case energies are the reported ASA energies; output-density components remain a separate diagnostic in the machine report.

| Case | q2 local up (Ry²) | down (Ry²) |
|---|---:|---:|
| RAW_PL | 0.050473635394 | 0.00847541394772 |
| LOCAL_PL | 0.0511128784368 | 0.0095449698681 |
| SCALAR_V_PLUS_001 | 0.050473635394 | 0.00847541394772 |
| VOUT_CONTROL | 0.0511128786714 | 0.00954496981847 |
| P_TRANSPOSE_CONTROL | 0.0511128784368 | 0.0095449698681 |

## ASA energy components

All values are frozen-input-rho case energies in Ry. Output-rho recomputed components are separately recorded in `a2_m03_cases.json`.

| Component | RAW | LOCAL | Scalar +0.01 | VOUT | P transpose |
|---|---:|---:|---:|---:|---:|
| EKIN | 2573.98262976 | 2574.23568485 | 2573.91652305 | 2574.23568484 | 2574.23568485 |
| ETOT | -2542.22381426 | -2541.97075917 | -2542.28992097 | -2541.97075917 | -2541.97075917 |
| RHOEPS | -108.155618094 | -108.155618094 | -108.155618094 | -108.155618094 | -108.155618094 |
| RHOMU | -142.848619974 | -142.848619974 | -142.848619974 | -142.848619974 | -142.848619974 |
| RHOVH | -3920.7481263 | -3920.7481263 | -3920.7481263 | -3920.7481263 | -3920.7481263 |
| SUMEC | -1487.1093672 | -1487.1093672 | -1487.1093672 | -1487.1093672 | -1487.1093672 |
| SUMEV | -2.50474932106 | -2.25169422835 | -2.57085603084 | -2.25169423744 | -2.25169422835 |
| UTOT | -5008.05082592 | -5008.05082592 | -5008.05082592 | -5008.05082592 | -5008.05082592 |
| VNUCL | 234.436674059 | 234.436674059 | 234.436674059 | 234.436674059 | 234.436674059 |
| ZVNUCL | -6095.35352555 | -6095.35352555 | -6095.35352555 | -6095.35352555 | -6095.35352555 |

## Control results and gate

Scalar +0.01 Ry gives `ΔSUMEV=0.06610670977426558 Ry`, expected `v·Σq0=0.06611715696608197 Ry`, absolute error `1.0447191816392865e-5 Ry` ≤ its locked `6.022982854018227e-5 Ry` allowance. Both local q2 values are nonnegative.

`V(Pout)` is reconstructed independently from the J=0 fixed-Liechtenstein formula `U_Ry(½I−Pσ)`, with manifest U=2 eV, J=0 and source conversion `ry2ev=13.605703976`. It matches production VOUT_CONTROL t, center, q0 and q2 to at most `4.44e-16` in t and `1.12e-16 Ry` in center; q0/q2 match exactly at serialized precision.

| Wrong control | Absolute carrier deviation from correct target (Ry) | A (Ry) | Discriminating under registered test? |
|---|---:|---:|---|
| VOUT_CONTROL | 0.00158593038781 | 0.00170223834723 | No |
| P_TRANSPOSE_CONTROL | 0.00158592129719 | 0.00170223834723 | No |

Both deviations are below A. Thus the mandatory negative controls are nondiscriminating and do not pass; `m03_all_gates_pass=false`. Overall KS M0 is not PASS. Stop before Stage M1.

## A1.3 / A2.6 measurement-based review proposal

The locked error model gives `A_true=0.001702238347 Ry`, while the measured wrong-control deviations are about `0.00158593 Ry`, short by approximately `0.00011631 Ry`. The J=0 VOUT operator differs only slightly from the accepted Hamiltonian input: signed d-shell `t` changes are `+1.06977e-8 Ry` (up) and `−1.41968e-9 Ry` (down); captured `Pout` is symmetric to the precision relevant here, so `Tr(Vin Poutᵀ)=Tr(Vin Pout)` exactly in the recorded calculation. The production radial replay error therefore masks these converged-state control differences. These results show that requiring each wrong control to exceed A is not a useful control gate when their independently measured operator separation is smaller than the radial response allowance.

Proposed future-only A1.3/A2.6 correction for user review: before any future control outcomes are inspected, compute each wrong-control separation `δ_k=Σ_s(t_{k,s}−t_s)` from the frozen matrices. Let `A_true` be the locked allowance at the correct shifts, `A_k` the same locked fit’s allowance at the wrong-control shifts, and `E_num` the predeclared replay roundoff guard. Require an energy-discrimination verdict only when `|δ_k| > A_true + A_k + E_num`; an eligible control must then violate `A_true` in production replay. If neither control is eligible, require independent operator-convention and provenance oracles; keep the gate incomplete if those oracles fail. The criterion protects whether a wrong pairing can be resolved above combined model and numerical error, and its error model is the sum of the locked true/control fit-plus-floor allowances and measured numerical replay error. This proposal introduces no synthetic radial setup; any energy control must use an accepted physical state and production `calculate_pl`/`newrho`. This proposal has not been applied and cannot retroactively change the present M0 result; user approval is required before future measurements.

Files: [`a2_m03_preregistration.json`](a2_m03_preregistration.json), [`a2_m03_cases.json`](a2_m03_cases.json), and [`a2_m03_cases.csv`](a2_m03_cases.csv).
