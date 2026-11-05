# A5.3 U=2 derivative comparison (interim)

The U comparison was opened after recording the fine-grid numerical allowance
in `a5_m22_numeric_allowance_preregistration.md`. This is the U=2 portion of
the test only; the J centers, local tightened U repeat, and full observables
table remain in progress, so it makes no production adoption decision.

Using the same-state candidate free energies at U=1.9, 1.95, 2.05 and 2.1 eV,
the central h=0.05 and 0.10 eV derivatives and their two-step Richardson
estimate are:

| radial mode / candidate | D(0.05) | D(0.10) | Richardson D | tightened D(0.05) |
|---|---:|---:|---:|---:|
| local b | +0.0550977903 | +0.0550980794 | +0.0550976939 | +0.0550977935 |
| spectral a (report only) | +0.0649282487 | +0.0647713781 | +0.0649805390 | +0.0649282324 |
| spectral c | −0.0836987238 | −0.0828785973 | −0.0839720992 | −0.0836987371 |

All derivatives are in Ry/eV. The local b derivative differs from its
production fixed-P `evaluate_hubbard_energy` target, +0.0552396891 Ry/eV, by
−0.0001419952 Ry/eV. This lies within the preregistered local ASA allowance
[0.00102115, 0.00155303] Ry/eV. The local h=0.05 versus h=0.10 derivative
difference is 2.89×10⁻⁷ Ry/eV. The local tightened h=0.05 slope differs from
the standard slope by 3.27×10⁻⁹ Ry/eV.

The wider central-difference grid h=0.05, 0.10, 0.20, 0.40 eV is retained in
the JSON. For local b, adjacent Richardson estimates at (.05,.10), (.10,.20),
and (.20,.40) are +0.0550976939, +0.0551019165, and +0.0550497136 Ry/eV.
For spectral c they are −0.0839720992, −0.0836614654, and −0.0814958600
Ry/eV. These wider steps diagnose curvature; the registered decision uses the
central h=.05/.10 estimate and the measured step sensitivity.

The spectral c derivative differs from its production fixed-P target,
+0.0578276954 Ry/eV, by −0.1417997946 Ry/eV. This miss is about 114 times the
upper spectral ASA allowance, 0.00124622 Ry/eV. The independently measured
spectral `d t/dU` is −0.1489526382 Ry/eV by the same two-step Richardson
formula. The c miss has the same sign and is about 5% smaller in magnitude;
the 0.00715 Ry/eV offset exceeds both the measured finite-difference
sensitivity and the ASA allowance, so no closer agreement is claimed or
attributed. The c-fail comparison remains decisive because its miss is far
outside the allowance.
Its h=0.05 versus h=0.10 derivative difference is 0.00082013 Ry/eV, below the
spectral ASA allowance; the standard-versus-tight h=0.05 result differs by
1.34×10⁻⁸ Ry/eV. Candidate a remains report-only, as preregistered.

The fixed-P target used direct final `P_OUT`: the local U=2 target P file
matches all 50 complex spin×d×d entries in the accepted U=2 records exactly.
The target output and P hashes, endpoint verification, and per-case source,
executable, input, records, stdout and restart hashes are in
`a5_m22_u2_derivative_result.json`. The corrected U=2 local endpoint has the
supplemental raw-log verification for the original `.true.` parser defect;
the original receipt remains unchanged. All compared candidate endpoints are
native reciprocal-KS runs with finite entropy/TS, literal convergence and
their applicable R05 status satisfied.

This supports the preregistered U-only result: b agrees within the measured
ASA base-functional accuracy and c fails by much more than that allowance.
The final M2 decision remains open pending the required J comparisons and
observables reporting.
