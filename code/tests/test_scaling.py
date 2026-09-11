"""Checks for the scientific relations, units and screening diagnostics."""

import unittest
import numpy as np
import sympy as sp
from nscomp.scaling import (
    CoreScaling, dimensionless_viscosity, exponent_envelope,
    knudsen_number, mean_spacing_m,
)


class ScalingChecks(unittest.TestCase):
    def test_time_elimination_against_independent_parameterization(self):
        tau = np.geomspace(1e-24, 1.0, 81)
        for h in [0, 0.001, 0.009, 0.01]:
            model = CoreScaling(2e-3, 7.0, h)
            speeds = 7 * tau ** (-0.5-h)
            expected_radii = 2e-3 * np.sqrt(tau)
            np.testing.assert_allclose(model.radius(speeds), expected_radii, rtol=2e-14)
            np.testing.assert_allclose(model.remaining_time_ratio(speeds), tau, rtol=2e-14)

    def test_symbolic_logarithmic_slope(self):
        s, h = sp.symbols("s h", real=True)
        slope = sp.simplify(sp.diff(s/2, s) / sp.diff(-(sp.Rational(1, 2)+h)*s, s))
        self.assertEqual(sp.simplify(slope + 1/(1+2*h)), 0)

    def test_inverse_over_many_decades(self):
        speed = np.geomspace(0.1, 3e8, 111)
        for h in [0, 0.005, 0.01]:
            model = CoreScaling(1e-3, 1, h)
            np.testing.assert_allclose(model.speed(model.radius(speed)), speed, rtol=1e-14)

    def test_anchor_and_monotonicity(self):
        model = CoreScaling(1e-4, 12, 0.009)
        self.assertEqual(float(model.radius(12)), 1e-4)
        self.assertTrue(np.all(np.diff(model.radius(np.geomspace(0.1, 1e8, 50))) < 0))

    def test_inverse_limit_known_hand_calculation(self):
        self.assertAlmostEqual(float(CoreScaling().radius(1e6)), 1e-9, delta=1e-23)

    def test_envelope_handles_both_sides_of_reference(self):
        speeds = np.array([0.1, 1, 10, 1e8])
        lo, hi = exponent_envelope(speeds)
        mid = CoreScaling(h=0.005).radius(speeds)
        self.assertTrue(np.all((lo <= mid) & (mid <= hi)))

    def test_units_rescale_consistently(self):
        # A change of numerical length and speed units must preserve the curve.
        model = CoreScaling(0.002, 3, 0.005)
        rescaled = CoreScaling(0.2, 300, 0.005)
        np.testing.assert_allclose(rescaled.radius(50000), 100*model.radius(500), rtol=1e-14)
        self.assertAlmostEqual(float(dimensionless_viscosity(1e-6, 1e-3, 1)), 1)

    def test_physical_length_diagnostics(self):
        self.assertAlmostEqual(float(knudsen_number(67.3e-9, 6.73e-6)), 0.01)
        spacing = float(mean_spacing_m(997.047013, 0.0180153, 6.02214076e23))
        self.assertTrue(0.30e-9 < spacing < 0.32e-9)

    def test_invalid_inputs_fail_explicitly(self):
        for value in [0, -1, np.nan, np.inf]:
            with self.assertRaises(ValueError):
                CoreScaling().radius(value)
        for h in [-0.1, 0.011, np.nan]:
            with self.assertRaises(ValueError):
                CoreScaling(h=h)


if __name__ == "__main__":
    unittest.main()
