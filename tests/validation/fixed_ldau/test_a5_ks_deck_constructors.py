#!/usr/bin/env python3
"""Every main-campaign case must explicitly construct a native-KS deck."""
from pathlib import Path
import re
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a5_ks_m22_m23_runner import cases, configure_case_deck


class MainDeckConstructorTests(unittest.TestCase):
    def test_every_case_selects_native_kspace(self):
        base = (Path(__file__).resolve().parent / "data" / "a5_common_seed" / "input.nml").read_text()
        for cfg in cases():
            with self.subTest(case=cfg["name"]):
                deck = configure_case_deck(base, cfg)
                self_block = re.search(r"(?ims)^\s*&self\b(.*?)/", deck)
                reciprocal_block = re.search(r"(?ims)^\s*&reciprocal\b(.*?)/", deck)
                self.assertIsNotNone(self_block)
                self.assertRegex(self_block.group(1), r"(?im)^\s*use_kspace\s*=\s*\.true\.\s*(?:!.*)?$")
                self.assertIsNotNone(reciprocal_block)
                self.assertRegex(reciprocal_block.group(1), r"(?im)^\s*nk1\s*=\s*12\s*(?:!.*)?$")


if __name__ == "__main__":
    unittest.main()
