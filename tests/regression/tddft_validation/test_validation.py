#!/usr/bin/env python3
"""CI entry point for the deterministic transverse TDDFT campaign fixture."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary:
        report = Path(temporary) / "evidence.json"
        result = subprocess.run(
            [sys.executable, str(ROOT / "tddft_validation.py"), str(ROOT / "fixtures" / "ci" / "campaign.json"),
             "--report", str(report)],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        if result.returncode:
            print(result.stdout)
            return result.returncode
        summary = json.loads(report.read_text())
    assert summary["goldstone"]["sum_rule_applied"] is True
    assert summary["goldstone"]["raw_residual"] > 0.0
    assert summary["stiffness"]["points"] == 2  # Gamma is excluded from D q^2 fitting.
    assert set(summary["routes"]) == {"LR-TDDFT", "GBT", "Jij"}
    assert "numerical-eta dependent" in summary["linewidth"]["interpretation"]
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
