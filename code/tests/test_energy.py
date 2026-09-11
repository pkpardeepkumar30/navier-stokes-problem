"""Independent energy integrals, volume preservation, and scaling checks."""
import math
import unittest
import numpy as np
import sympy as sp
from nscomp.energy import (
    CoreEnergyScaling, cylindrical_energy_quadrature, gaussian_energy_fraction,
    gaussian_swirl, material_map, narrow_gaussian, narrow_gaussian_squared_integral,
    swirl_cartesian,
)


class EnergyChecks(unittest.TestCase):
    def test_cylindrical_integral_against_analytic_value(self):
        for h in (0, 1e-6, 0.1):
            model = CoreEnergyScaling(h=h, a0=1.7, b0=3.1, u0=2.3, density=0.8)
            for s in (1, 0.03, 1e-8):
                exact = model.toy_energy(s)*gaussian_energy_fraction(8)
                self.assertAlmostEqual(cylindrical_energy_quadrature(model, s)/exact,
                                       1, delta=2e-12)

    def test_cartesian_and_cylindrical_integrals_agree(self):
        # A separately assembled Cartesian tensor-product integral catches a
        # missing cylindrical r factor or a mistaken circumference factor.
        model = CoreEnergyScaling(a0=0.7, b0=1.9, u0=2.1, density=1.3)
        q, w = np.polynomial.legendre.leggauss(72)
        x, y, z = q*8*model.a0, q*8*model.a0, q*8*model.b0
        u = swirl_cartesian(x[:, None, None], y[None, :, None], z[None, None, :],
                            model.a0, model.b0, model.u0)
        weights = w[:, None, None]*w[None, :, None]*w[None, None, :]
        integral = model.density/2*np.sum(np.sum(u*u, axis=-1)*weights)
        integral *= (8*model.a0)**2*(8*model.b0)
        self.assertAlmostEqual(integral/model.toy_energy(1), 1, delta=3e-12)

    def test_symbolic_jacobian_and_energy_power(self):
        s, h = sp.symbols("s h", positive=True)
        volume = s**sp.Rational(1,2)*s**sp.Rational(1,2)*s**(sp.Rational(1,2)-h)
        energy = volume*(s**(-sp.Rational(1,2)-h))**2
        self.assertEqual(sp.simplify(energy/s**(sp.Rational(1,2)-3*h)), 1)

    def test_energy_regimes_are_not_confused_with_peak(self):
        s = np.geomspace(1, 1e-12, 100)
        for h, direction in ((0, -1), (0.1, -1), (1/6, 0), (0.2, 1)):
            ratios = CoreEnergyScaling(h=h).ratios(s)
            self.assertTrue(np.all(np.diff(ratios["speed"]) > 0))
            delta = np.diff(ratios["energy"])
            if direction == 0:
                np.testing.assert_allclose(ratios["energy"], 1, rtol=1e-14)
            else:
                self.assertTrue(np.all(direction*delta > 0))

    def test_one_dimensional_integral(self):
        q, w = np.polynomial.legendre.leggauss(80)
        for epsilon in (1, 0.01, 1e-5):
            for exponent in (0.25, 0.5, 0.75):
                x = epsilon*8*q
                numeric = np.sum(w*narrow_gaussian(x, epsilon, exponent)**2)*epsilon*8
                exact = narrow_gaussian_squared_integral(epsilon, exponent)
                self.assertAlmostEqual(numeric/exact, 1, delta=2e-12)

    def test_swirl_is_smooth_on_axis_and_divergence_free(self):
        x, y, z, a, b, U = sp.symbols("x y z a b U", real=True)
        factor = U/a*sp.exp((1-(x*x+y*y)/a**2-z*z/b**2)/2)
        divergence = sp.diff(-y*factor, x)+sp.diff(x*factor, y)
        self.assertEqual(sp.simplify(divergence), 0)
        np.testing.assert_array_equal(swirl_cartesian(0, 0, 0.7, 1, 2, 3), [0, 0, 0])
        self.assertAlmostEqual(float(gaussian_swirl(1.7, 0, 1.7, 2, 3.2)), 3.2)

    def test_material_volume_and_pathline_equation(self):
        point = np.array([0.8, -0.4, 0.3])
        for s in (1, 0.3, 0.01):
            transform = material_map(np.eye(3), s)
            self.assertAlmostEqual(np.linalg.det(transform), 1, delta=5e-15)
            dt = 1e-5*s
            # Increasing physical time decreases s=1-t.
            derivative = (material_map(point, s-dt)-material_map(point, s+dt))/(2*dt)
            position = material_map(point, s)
            expected = position*np.array([-0.5, -0.5, 1])/s
            np.testing.assert_allclose(derivative, expected, rtol=3e-10)

    def test_invalid_scales_and_time(self):
        for bad in (0, -1, np.nan, np.inf):
            with self.assertRaises(ValueError):
                CoreEnergyScaling(a0=bad)
            with self.assertRaises(ValueError):
                CoreEnergyScaling().ratios(bad)
        with self.assertRaises(ValueError):
            CoreEnergyScaling(h=0.5)
        with self.assertRaises(ValueError):
            gaussian_swirl(-1, 0, 1, 1, 1)


if __name__ == "__main__":
    unittest.main()
