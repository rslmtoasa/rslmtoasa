import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("c2_hubbard_operator_audit.py")
SPEC = importlib.util.spec_from_file_location("c2_hubbard_operator_audit", MODULE_PATH)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def emit_record(assembler, site, correction_shift=0.0):
    lines = [f"META {assembler} 1 {site} 2 Fe 1.39650987654321000E-01 2.94107342100000000E-02 1"]
    for spin in (1, 2):
        for i in range(1, 6):
            for j in range(1, 6):
                value = (spin * 100 + i * 10 + j) / 1000
                lines.append(f"OCC {assembler} 1 {site} 2 {spin} {i} {j} {value:.17E} 0.0")
    for spin in (1, 2):
        for i in range(1, 6):
            for j in range(1, 6):
                correction = (i + j + spin) / 100 + correction_shift
                before = (i - j) / 100
                after = before + correction
                lines.append(
                    f"BLOCK {assembler} 1 {site} 2 {spin} {i} {j} 1 "
                    f"{before:.17E} 0.0 {correction:.17E} 0.0 {after:.17E} 0.0"
                )
    return lines


class C2AuditParserTest(unittest.TestCase):
    def load_pair(self, local_shift=0.0):
        contents = emit_record("BULK_EE", 4) + emit_record("LOCAL_HALL", 9, local_shift)
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "audit.dat"
            path.write_text("\n".join(contents) + "\n")
            records = AUDIT.read_records(path)
        bulk = AUDIT.choose_record(records, "BULK_EE", 4, 1)
        local = AUDIT.choose_record(records, "LOCAL_HALL", 9, 1)
        return bulk, local

    def test_equal_bulk_and_local_operator_records_pass(self):
        bulk, local = self.load_pair()
        result = AUDIT.compare_records(bulk, local, 1.0e-12)
        self.assertEqual(result["insertion_counts_match"], 1)
        self.assertEqual(result["correction_max_abs_diff"], 0.0)

    def test_operator_difference_is_rejected(self):
        bulk, local = self.load_pair(local_shift=1.0e-5)
        with self.assertRaises(AssertionError):
            AUDIT.compare_records(bulk, local, 1.0e-12)


if __name__ == "__main__":
    unittest.main()
