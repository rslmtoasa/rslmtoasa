#!/usr/bin/env python3
"""Regression checks for native `Converged!` marker parsing."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a5_ks_m22_m23_runner import (
    fortran_true,
    has_converged_marker as main_marker,
)
from a5_local_w1_controls import has_converged_marker as w1_marker


class ConvergenceMarkerTests(unittest.TestCase):
    def test_r05_fortran_true_spellings(self):
        self.assertTrue(fortran_true(".true."))
        self.assertTrue(fortran_true("true"))
        self.assertFalse(fortran_true(".false."))
        self.assertFalse(fortran_true("false"))

    def test_not_converged_is_not_a_positive_marker(self):
        log = "Not converged! Diff=1.0e-8\n"
        for parser in (main_marker, w1_marker):
            with self.subTest(parser=parser.__module__):
                self.assertFalse(parser(log))

    def test_actual_converged_marker_is_positive(self):
        log = "Converged! Diff=1.0e-10\n"
        for parser in (main_marker, w1_marker):
            with self.subTest(parser=parser.__module__):
                self.assertTrue(parser(log))

    def test_mixed_log_requires_and_recognizes_positive_marker(self):
        log = "Not converged! Diff=1.0e-8\nConverged! Diff=1.0e-10\n"
        for parser in (main_marker, w1_marker):
            with self.subTest(parser=parser.__module__):
                self.assertTrue(parser(log))


if __name__ == "__main__":
    unittest.main()
