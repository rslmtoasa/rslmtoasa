# C6 plan

T0 selects the sole existing fixture, `fccCu001`, and records its host, recursion,
and Fermi setup. T1 adds `calctype='S'` to the existing local-mode eligibility
path, tests ownership-based frozen-host selection, and runs the requested gates.
T2/T3 use accepted local bulk hosts of the same material, with at most four SCF
runs total and a preflight before each. T4 records the workflow and open issues,
then merges locally only if the required gates pass. If no matching accepted host
or valid smoke deck exists, preserve the evidence and apply C5.3 item 2.
