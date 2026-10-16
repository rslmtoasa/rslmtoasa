#!/usr/bin/env python3
"""The fourteen preregistered W1 cases must use native KS and proper spin selectors."""
from pathlib import Path
import re
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a5_spin_w1_probe import case_specs, configure_case_deck, reconstruct_local_free_energy


class SpinW1DeckConstructorTests(unittest.TestCase):
    def test_external_field_free_energy_uses_same_state_trace_and_entropy(self):
        self.assertAlmostEqual(reconstruct_local_free_energy(-10.0, -0.005, 4.0, 0.2), -10.22)

    def test_all_fourteen_cases_are_spin_selective_native_ks(self):
        base = (Path(__file__).resolve().parent / "data/a5_common_seed/input.nml").read_text()
        specs = case_specs()
        self.assertEqual(len(specs), 14)
        self.assertEqual(len({case["name"] for case in specs}), 14)
        for spec in specs:
            with self.subTest(case=spec["name"]):
                deck = configure_case_deck(base, spec)
                self_block = re.search(r"(?ims)^\s*&self\b(.*?)/", deck).group(1)
                ham_block = re.search(r"(?ims)^\s*&hamiltonian\b(.*?)/", deck).group(1)
                self.assertRegex(self_block, r"(?im)^\s*use_kspace\s*=\s*\.true\.\s*(?:!.*)?$")
                self.assertRegex(ham_block, rf"(?im)^\s*hubbard_test_spin\s*=\s*{spec['spin']}\s*(?:!.*)?$")
                self.assertRegex(ham_block, r"(?im)^\s*hubbard_u_general\(1,3\)\s*=\s*0\.0\s*(?:!.*)?$")
                self.assertRegex(ham_block, r"(?im)^\s*hubbard_j_general\(1,3\)\s*=\s*0\.0\s*(?:!.*)?$")


if __name__ == "__main__":
    unittest.main()
