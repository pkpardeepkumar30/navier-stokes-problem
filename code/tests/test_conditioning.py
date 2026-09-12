"""Accuracy against high-precision integrals and a controlled input perturbation."""
import unittest
import numpy as np
import mpmath as mp
from nscomp.moments import SmoothBump
from nscomp.conditioning import (axial_system,stable_axial_solve,transformed_axial_system,
    high_precision_solve,coefficient_relative_error,reference_moment_remainder)


class ConditioningChecks(unittest.TestCase):
    def setUp(self):
        self.bumps=[SmoothBump(c,.12) for c in (2.,3.)]
        self.d=np.array([1.,-.4])

    def test_transformation_agrees_when_original_is_resolved(self):
        original=np.linalg.solve(axial_system(.1,self.bumps),-self.d)
        stable=stable_axial_solve(.1,self.bumps,self.d)
        np.testing.assert_allclose(original,stable,rtol=3e-13)

    def test_small_lambda_accuracy_against_independent_reference(self):
        _,reference=high_precision_solve(1e-12,(2.,3.),.12,self.d,60)
        stable=stable_axial_solve(1e-12,self.bumps,self.d)
        original=np.linalg.solve(axial_system(1e-12,self.bumps),-self.d)
        good=coefficient_relative_error(stable,reference,60)
        bad=coefficient_relative_error(original,reference,60)
        self.assertLess(good,2e-13)
        self.assertGreater(bad,good*1000)

    def test_more_precision_resolves_original_moment_rows(self):
        _,reference=high_precision_solve(1e-12,(2.,3.),.12,self.d,60)
        _,coarse=high_precision_solve(1e-12,(2.,3.),.12,self.d,20)
        _,fine=high_precision_solve(1e-12,(2.,3.),.12,self.d,40)
        self.assertLess(coefficient_relative_error(fine,reference,60),1e-25)
        self.assertGreater(coefficient_relative_error(coarse,reference,60),1e-13)

    def test_transformed_rows_remain_well_conditioned(self):
        for lam in (1e-8,1e-12,1e-18):
            B,_=transformed_axial_system(lam,self.bumps,self.d)
            self.assertLess(np.linalg.cond(B),20)

    def test_reference_precision_and_quadrature_refinement(self):
        _,coarse=high_precision_solve(1e-18,(2.,3.),.12,self.d,70)
        _,fine=high_precision_solve(1e-18,(2.,3.),.12,self.d,90)
        self.assertLess(coefficient_relative_error(coarse,fine,90),1e-45)
        a=stable_axial_solve(1e-12,self.bumps,self.d,128)
        b=stable_axial_solve(1e-12,self.bumps,self.d,256)
        np.testing.assert_allclose(a,b,rtol=2e-13)

    def test_accurate_algorithm_does_not_remove_input_sensitivity(self):
        d=np.array([1.,1.])
        baseline=stable_axial_solve(1e-10,self.bumps,d)
        perturbed=stable_axial_solve(1e-10,self.bumps,d+[0,1e-10])
        self.assertGreater(np.linalg.norm(perturbed-baseline)/np.linalg.norm(baseline),.1)

    def test_rounding_coefficients_changes_physical_moment_error(self):
        matrix,reference=high_precision_solve(1e-12,(2.,3.),.12,self.d,60)
        precise=reference_moment_remainder(matrix,reference,self.d,60)
        rounded=reference_moment_remainder(matrix,np.array([float(x) for x in reference]),self.d,60)
        self.assertLess(precise,1e-40)
        self.assertGreater(rounded,1e-8)


if __name__=='__main__':
    unittest.main()
