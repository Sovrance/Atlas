import math
import os
import sys
import unittest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))
import core


class RHWeilExactTests(unittest.TestCase):
    def test_normalization_audit_constant(self):
        self.assertEqual(core.NORMALIZATION, "G = G0 - Gp + Ginf")
        self.assertIn("no_rh", core.CLAIM_BOUNDARY)

    def test_overlap_low_degree(self):
        for L, a in [(1.3, 0.7), (2.0, 0.4), (1.1, 0.2)]:
            self.assertAlmostEqual(core.overlap_c(0, 0, a, L) * 2, core.kernel_00(a, L))
            k01 = core.overlap_c(0, 1, a, L) + core.overlap_c(1, 0, a, L)
            self.assertAlmostEqual(k01, core.kernel_01(a, L))
            self.assertAlmostEqual(2 * core.overlap_c(1, 1, a, L), core.kernel_11(a, L))
            self.assertAlmostEqual(core.kernel_ij(0, 1, a, L), core.kernel_01(a, L))

    def test_q1_formula_by_basis_transform(self):
        L, a = 1.31, 0.69
        self.assertAlmostEqual(core.midpoint_reflection_q1(L, a), core.kernel_q1q1(a, L))

    def test_q1_sign_threshold(self):
        thr = core.q1_sign_threshold()
        L = 2.0
        a_neg = thr * L + 1e-9
        a_pos = thr * L - 1e-9
        self.assertLess(core.kernel_q1q1(a_neg, L), 0)
        self.assertGreater(core.kernel_q1q1(a_pos, L), 0)

    def test_bubble_det(self):
        L, a = 1.31, 0.69
        lhs = core.kernel_00(a, L) * core.kernel_bb(a, L) - core.kernel_0b(a, L) ** 2
        self.assertAlmostEqual(lhs, core.kernel_bubble_det(a, L))

    def test_bubble_det_threshold(self):
        thr = core.bubble_det_threshold()
        L = 2.0
        self.assertLess(core.kernel_bubble_det(thr * L + 1e-9, L), 0)
        self.assertGreater(core.kernel_bubble_det(thr * L - 1e-9, L), 0)

    def test_current_cell_q1_prime_sign(self):
        for a in [math.log(2), math.log(3)]:
            for L in [math.log(3), 1.2, math.log(4)]:
                if a < L:
                    self.assertLess(core.kernel_q1q1(a, L), 0)

    def test_scalar_curvature_positive_current_cell(self):
        for L in [math.log(3), 1.2, math.log(4)]:
            self.assertGreater(core.scalar_curvature(L), 0)

    def test_degree2_parity_factorization(self):
        # Synthetic positive pivots exercising the algebraic identities.
        g00, o1, e2, L = 1.5, 0.7, 2.25, 1.3
        d2 = core.degree2_raw_det(g00, o1, e2, L)
        self.assertAlmostEqual(d2, e2 + L * L * g00 * o1)
        self.assertAlmostEqual(core.degree2_full_det(o1, e2), o1 * e2)


if __name__ == "__main__":
    unittest.main()
