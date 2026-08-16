import math, os, sys, unittest
HERE=os.path.dirname(__file__)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'src'))
import core

class RHWeilExactTests(unittest.TestCase):
    def test_overlap_low_degree(self):
        for L,a in [(1.3,.7),(2.0,.4),(1.1,.2)]:
            self.assertAlmostEqual(core.overlap_c(0,0,a,L)*2, core.kernel_00(a,L))
            k01=core.overlap_c(0,1,a,L)+core.overlap_c(1,0,a,L)
            self.assertAlmostEqual(k01, core.kernel_01(a,L))
            self.assertAlmostEqual(2*core.overlap_c(1,1,a,L), core.kernel_11(a,L))
    def test_q1_formula_by_basis_transform(self):
        L,a=1.31,.69
        k00=core.kernel_00(a,L); k01=core.kernel_01(a,L); k11=core.kernel_11(a,L)
        transformed=k11-L*k01+(L*L/4)*k00
        self.assertAlmostEqual(transformed, core.kernel_q1q1(a,L))
    def test_bubble_det(self):
        L,a=1.31,.69
        lhs=core.kernel_00(a,L)*core.kernel_bb(a,L)-core.kernel_0b(a,L)**2
        self.assertAlmostEqual(lhs, core.kernel_bubble_det(a,L))
    def test_current_cell_q1_prime_sign(self):
        for a in [math.log(2), math.log(3)]:
            for L in [math.log(3),1.2,math.log(4)]:
                if a < L:
                    self.assertLess(core.kernel_q1q1(a,L),0)
    def test_scalar_curvature_positive_current_cell(self):
        for L in [math.log(3),1.2,math.log(4)]:
            self.assertGreater(core.scalar_curvature(L),0)

if __name__=='__main__': unittest.main()
