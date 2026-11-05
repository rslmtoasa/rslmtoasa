import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[3] / "docs/validation/fixed_ldau/rs/C5_4_PREFLIGHT.py"
spec = importlib.util.spec_from_file_location("c5_4_preflight", SCRIPT)
preflight = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preflight)


class HostFermiPreflightTest(unittest.TestCase):
    def test_matches_host_output_fermi_within_tolerance(self):
        with tempfile.TemporaryDirectory() as directory:
            trace = Path(directory) / "host_trace.csv"
            trace.write_text("kind,scalar\nfermi_level,-0.082949627451487762\n")
            result = preflight.check_host_fermi(
                {"energy": {"fermi": -0.0829497}}, trace
            )
            self.assertEqual(result["status"], "PASS")
            self.assertLessEqual(abs(result["delta_ry"]), 1e-6)
            self.assertEqual(result["host_ef_trace"], str(trace))

    def test_rejects_input_fermi_outside_tolerance(self):
        with tempfile.TemporaryDirectory() as directory:
            trace = Path(directory) / "host_trace.csv"
            trace.write_text("kind,scalar\nfermi_level,-0.082949627451487762\n")
            with self.assertRaises(SystemExit):
                preflight.check_host_fermi({"energy": {"fermi": -0.085152}}, trace)


if __name__ == "__main__":
    unittest.main()
