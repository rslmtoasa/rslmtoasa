#!/usr/bin/env python3
"""Static WP10 guard: only gbt_bond_phase (via gbt_endpoint_link) may
compute a q-dotted-into-a-bond phase anywhere in source/.

The one deliberate exception is prepare_explicit_texture_moments in
hamiltonian_build.f90, which computes a different physical quantity -- the
absolute site-position phase that gives explicit_texture its per-site
moment direction, not the bond gauge link used by gbt_single_q. That
routine is named explicitly below and checked for both the presence of the
pattern (so this test fails loudly if the exemption becomes stale) and its
absence from everywhere else.
"""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source"
EXEMPT_FILE = "hamiltonian_build.f90"
EXEMPT_ROUTINE = "prepare_explicit_texture_moments"

# q_ss(1)/q_ss(2)/q_ss(3) used as a component-wise multiplier is the
# hand-rolled q.bond pattern gbt_bond_phase exists to centralize.
QDOT_PATTERN = re.compile(r"q_ss\(\s*[123]\s*\)")


def main() -> None:
    failures: list[str] = []

    for path in sorted(SOURCE.rglob("*.f90")):
        rel = path.relative_to(ROOT).as_posix()
        if path.name == "gbt_structure.f90":
            continue

        text = path.read_text(encoding="utf-8")

        if path.name == EXEMPT_FILE:
            start = text.find(f"module subroutine {EXEMPT_ROUTINE}")
            end = text.find(f"end subroutine {EXEMPT_ROUTINE}")
            if start == -1 or end == -1:
                failures.append(f"{rel}: {EXEMPT_ROUTINE} not found; update this test's exemption")
                continue
            exempt_body = text[start:end]
            if not QDOT_PATTERN.search(exempt_body):
                failures.append(
                    f"{rel}: {EXEMPT_ROUTINE} no longer contains the expected site-phase "
                    "pattern; narrow or remove this test's exemption"
                )
            text = text[:start] + text[end:]

        if QDOT_PATTERN.search(text):
            failures.append(f"{rel}: computes q.bond outside gbt_bond_phase/gbt_endpoint_link")

    if failures:
        raise SystemExit("\n".join(f"FAIL: {item}" for item in failures))
    print("WP10 source contract PASS: only gbt_bond_phase computes q.bond")


if __name__ == "__main__":
    main()
