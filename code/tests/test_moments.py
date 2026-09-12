"""Independent integration checks for the finite moment-repair substep."""
import unittest
import numpy as np
from scipy.integrate import quad
from nscomp.moments import SmoothBump, moment_matrices, solve_moments, original_moments, repaired_moments


class MomentChecks(unittest.TestCase):
    def setUp(self):
        self.u=[SmoothBump(c,.12) for c in (2,3)]
        self.e=[SmoothBump(c,.12) for c in (4,5,6)]
        self.du=np.array([1,-.4])
        self.de=np.array([.3,-.2,.1])

    def integrate(self,bump,weight):
        return quad(lambda r:float(bump.value(r))*weight(r),bump.center-bump.halfwidth,
                    bump.center+bump.halfwidth,epsabs=1e-12,epsrel=1e-12)[0]

    def test_bump_normalization_and_support(self):
        b=self.u[0]
        self.assertAlmostEqual(self.integrate(b,lambda r:1),1,delta=2e-13)
        np.testing.assert_array_equal(b.value(np.array([0,b.center-b.halfwidth,b.center+b.halfwidth,9])),0)

    def test_powers_against_adaptive_quadrature(self):
        for b in (self.u[0],self.e[-1]):
            for power in (1,.8,2,-2.2,-.2):
                exact=self.integrate(b,lambda r:r**power)
                self.assertAlmostEqual(b.moment(power)/exact,1,delta=2e-13)

    def test_all_five_physical_moment_changes(self):
        lam,ef=.1,1.3
        alpha,beta=solve_moments(lam,self.u,self.e,self.du,self.de)
        # Integrate the reserved-patch formulas directly, independently of BU/BE.
        E0=lambda r:ef*r**(-1-2*lam)
        d1=sum(a*self.integrate(b,lambda r:r) for a,b in zip(alpha,self.u))
        d2=sum(a*self.integrate(b,lambda r:r*r) for a,b in zip(beta,self.e))
        d3=sum(a*self.integrate(b,lambda r:2*E0(r)/r) for a,b in zip(beta,self.e))
        d4=sum(a*self.integrate(b,lambda r:r*r*E0(r)) for a,b in zip(alpha,self.u))
        d5=sum(a*self.integrate(b,lambda r:-r*E0(r)) for a,b in zip(beta,self.e))
        np.testing.assert_allclose(original_moments(self.du,self.de,ef)+[d1,d2,d3,d4,d5],0,atol=3e-12)

    def test_linear_response_to_discrepancies(self):
        a,b=solve_moments(.1,self.u,self.e,self.du,self.de)
        aa,bb=solve_moments(.1,self.u,self.e,-2*self.du,-2*self.de)
        np.testing.assert_allclose(aa,-2*a,rtol=1e-13)
        np.testing.assert_allclose(bb,-2*b,rtol=1e-13)

    def test_reference_moments_cancel(self):
        alpha,beta=solve_moments(.1,self.u,self.e,self.du,self.de,96)
        BU,BE=moment_matrices(.1,self.u,self.e,384)
        after=repaired_moments(BU,BE,alpha,beta,original_moments(self.du,self.de))
        self.assertLess(float(np.max(abs(after))),1e-12)

    def test_distinct_powers_can_be_ill_conditioned(self):
        b1,_=moment_matrices(.1,self.u,self.e)
        b2,_=moment_matrices(1e-6,self.u,self.e)
        self.assertGreater(np.linalg.cond(b2)/np.linalg.cond(b1),1e4)

    def test_invalid_degenerate_setup(self):
        with self.assertRaises(ValueError): moment_matrices(0,self.u,self.e)
        with self.assertRaises(ValueError): moment_matrices(.1,[self.u[0],self.u[0]],self.e)


if __name__=='__main__':
    unittest.main()
