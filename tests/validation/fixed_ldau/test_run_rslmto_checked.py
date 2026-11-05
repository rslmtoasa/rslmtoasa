#!/usr/bin/env python3
"""Safety tests for checked SCF launches; the successful solver is a no-op stub."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[3]
WRAPPER = ROOT / "tools/run_rslmto_checked.py"
PYTHON = ROOT / ".venv/bin/python3"
DECK_FIXTURE = ROOT / "tests/scf/cases/surface/fccCu001/input.nml"
GUARD = ROOT / "tests/validation/fixed_ldau/c5_7_launch.sh"
class CheckedLauncherTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(prefix="checked-rslmto-")
        self.root = Path(self.temp.name)
        self.binary = self.root / "rslmto.x"
        self.binary.write_text(
            "#!/bin/sh\n"
            "printf 'stub launch\\n' >> \"$CHECKED_LAUNCH_MARKER\"\n"
        )
        self.binary.chmod(0o755)
        self.trace = self.root / "host_ef.csv"
        self.trace.write_text("kind,scalar\nfermi_level,-0.088376\n")
        self.marker = self.root / "solver_launches.txt"

    def tearDown(self) -> None:
        self.temp.cleanup()

    def run_deck(self, name: str, deck_text: str) -> tuple[subprocess.CompletedProcess[str], Path]:
        run_dir = self.root / name
        run_dir.mkdir()
        deck = run_dir / "input.nml"
        deck.write_text(deck_text)
        env = os.environ.copy()
        env["CHECKED_LAUNCH_MARKER"] = str(self.marker)
        result = subprocess.run(
            [str(PYTHON), str(WRAPPER), str(run_dir), "--binary", str(self.binary),
             "--host-ef-trace", str(self.trace), "--guard", str(GUARD)],
            cwd=ROOT, env=env, text=True, capture_output=True, check=False,
        )
        return result, run_dir

    def test_unclosed_namelist_is_refused_without_launch(self) -> None:
        deck_text = DECK_FIXTURE.read_text().rsplit("/\n", 1)[0] + "\n"
        result, run_dir = self.run_deck("unclosed", deck_text)
        self.assertNotEqual(result.returncode, 0)
        record = json.loads((run_dir / "preflight.json").read_text())
        self.assertEqual(record["status"], "PREFLIGHT_FAILED")
        self.assertEqual(record["launch"]["status"], "REFUSED")
        self.assertFalse(self.marker.exists())

    def test_host_fermi_mismatch_is_refused_without_launch(self) -> None:
        deck_text = DECK_FIXTURE.read_text().replace("fermi = -0.088376", "fermi = -0.080000")
        result, run_dir = self.run_deck("fermi-mismatch", deck_text)
        self.assertNotEqual(result.returncode, 0)
        record = json.loads((run_dir / "preflight.json").read_text())
        self.assertEqual(record["status"], "PREFLIGHT_FAILED")
        self.assertIn("differs from host EF", record["preflight_stderr"])
        self.assertEqual(record["launch"]["status"], "REFUSED")
        self.assertFalse(self.marker.exists())

    def test_valid_deck_launches_noop_solver_stub(self) -> None:
        result, run_dir = self.run_deck("valid", DECK_FIXTURE.read_text())
        self.assertEqual(result.returncode, 0, result.stderr)
        record = json.loads((run_dir / "preflight.json").read_text())
        self.assertEqual(record["status"], "RUN_COMPLETED")
        self.assertEqual(record["preflight_result"]["status"], "PASS")
        self.assertEqual(record["preflight_result"]["host_ef_check"]["delta_ry"], 0.0)
        self.assertEqual(record["launch"]["returncode"], 0)
        self.assertEqual(self.marker.read_text(), "stub launch\n")

    def test_runtime_chebyshev_keyword_passes_preflight(self) -> None:
        deck_text = DECK_FIXTURE.read_text().replace("recur = 'block'", "recur = 'chebyshev'")
        result, run_dir = self.run_deck("chebyshev", deck_text)
        self.assertEqual(result.returncode, 0, result.stderr)
        record = json.loads((run_dir / "preflight.json").read_text())
        self.assertEqual(record["status"], "RUN_COMPLETED")
        self.assertEqual(record["preflight_result"]["recur"], "chebyshev")
        self.assertEqual(self.marker.read_text(), "stub launch\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
