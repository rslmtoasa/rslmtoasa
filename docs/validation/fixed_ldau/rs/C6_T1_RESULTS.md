# C6 T1 — surface auto-local resolver

The resolver now includes ``calctype='S'`` in the existing fixed-Liechtenstein,
scalar-collinear d-shell auto-local eligibility scope. It still requires every
U-carrying frozen host potential to carry local provenance. The fallback reason
is also logged for surfaces. ``is_eligible_local_hubbard`` retains its explicit
bulk-only condition, so Hubbard energy routing remains off for surfaces.

``is_frozen_hubbard_host`` already accepts ``I`` and ``S``. It determines the
active recursion types from ``irec -> iz`` and the frozen host types from
``ib -> iz``; the surface host classification follows the same ownership maps.

The focused ``UnitAutoImpurityLocalMode`` target rebuilt and passed (1/1). It
checks the impurity fixture as before and the new ``calctype='S'`` deck with a
local U-carrying host through the shared eligibility resolver. No SCF run was
performed. The 13 RA2 guard scenarios and single-thread no-U RS/KS identity have
prior PASS evidence on the parent tree (C5.13); they were not rerun for C6.
Those checks require solver launches, and the four-run cap is retained.

T1 code review and the focused unit gate pass. T2/T3 remain unevaluated because
the only surface fixture has no matching accepted local U-carrying Cu host; see
``C6_T0_RESULTS.md``. Under C5.3 item 2, no surface smoke result or adoption is
claimed and no merge is made.
