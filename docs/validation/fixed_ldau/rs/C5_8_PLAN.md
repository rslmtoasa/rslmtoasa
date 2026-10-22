# C5.8 plan

One Luna worker executes sequential cards; one batched audit precedes merge. Keep separate ledgers for at most twelve short reproductions and three campaigns; required regression/merge verification remains separate. Preflight every invocation, cap identity threads, monitor disk, preserve references and unrelated jobs, and do not push.

1. Reproduce the exact T2a failure with current Release, the prescribed checked RelWithDebInfo flags, Debug, and both historical Release revisions. Stop each at its first completed iteration or crash. Record only necessary old-parser line removals and exact build/source/input hashes.
2. If 908074f passes, bisect source-touching commits within the remaining reproduction budget. Otherwise retain honest comparability limitations; never infer a pre-existing crash from an input rejection.
3. Fix only a demonstrated programming defect. Verify the T2a and RD4 zero-lambda Release cases, no-U RS/KS byte identity and thirteen guard fixtures. If both historical revisions reproduce the crash, document it as pre-existing and use Debug as authorized.
4. Regenerate accepted KS Fe local/KS provenance, run the Fe RS impurity auto smoke, then complete density-only auto adoption checks and documentation. Impurity energy routing stays off.
5. Audit all gates once; if supported, perform C5.2 §7 local merge and merged-tree checks. On the specified failure or an uncovered outcome, preserve criteria/defaults, commit records and report five lines instead.
