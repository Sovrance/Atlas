import math
import os
import sys
import unittest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))

import cells
import core
import scalar


class TestPrimePowerCells(unittest.TestCase):
    def test_split_log3_log4(self):
        cs = cells.split_cells(math.log(3), math.log(4))
        self.assertEqual(len(cs), 1)
        c = cs[0]
        self.assertAlmostEqual(c.L_left, math.log(3))
        self.assertAlmostEqual(c.L_right, math.log(4))

    def test_breaks_include_3_and_4(self):
        br = cells.prime_powers_in_log_interval(math.log(3), math.log(4))
        qs = [q for q, _, _ in br]
        self.assertIn(3, qs)
        self.assertIn(4, qs)


class TestScalarVerifier(unittest.TestCase):
    def test_current_cell_e0(self):
        report = scalar.verify_scalar_cell()
        self.assertTrue(report.w00_second_positive)
        self.assertEqual(report.evidence_class, "E0")
        self.assertFalse(report.rh_proof_claim)
        self.assertTrue(report.at_most_one_interior_minimizer)
        self.assertIsNotNone(report.left_jump)
        self.assertIsNotNone(report.right_jump)
        self.assertLess(report.left_jump, 0)
        self.assertLess(report.right_jump, 0)

    def test_jump_formula_matches_core(self):
        self.assertAlmostEqual(core.von_mangoldt_jump(3, 3), 2 * math.log(3) / math.sqrt(3))
        self.assertAlmostEqual(core.von_mangoldt_jump(4, 2), 2 * math.log(2) / 2)

    def test_grid_positive_diagnostic(self):
        vals = scalar.sample_curvature_grid(n=11)
        self.assertTrue(all(v > 0 for v in vals))


if __name__ == "__main__":
    unittest.main()
