"""Independent checks of scale elimination, threshold ordering, and physical ratios."""
import unittest
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from nscomp.scaling import CoreScaling
from nscomp.regimes import (GasReference, gas_diagnostics, kn_crossing_speed,
    rate_crossing_speed, carrier_power_length, averaging_cube_count)


class RegimeChecks(unittest.TestCase):
    def setUp(self):
        self.gas = GasReference(296.15, 101300., 287.05, 1.4, 1.83245e-5, 6.73e-8)
        self.model = CoreScaling()

    def test_thermal_speed_against_maxwell_integral(self):
        rt = self.gas.gas_constant_J_kg_K * self.gas.temperature_K
        # Integrate the normalized dimensionless Maxwell speed density.
        norm = quad(lambda x: 4/np.sqrt(np.pi)*x*x*np.exp(-x*x), 0, np.inf)[0]
        mean = quad(lambda x: 4/np.sqrt(np.pi)*x**3*np.exp(-x*x), 0, np.inf)[0]
        self.assertAlmostEqual(norm, 1, places=12)
        self.assertAlmostEqual(mean*np.sqrt(2*rt)/self.gas.mean_thermal_speed, 1, places=12)

    def test_carrier_against_independent_scale_parameter(self):
        q = np.geomspace(1e-18, 1, 27)
        for h in (0, 1e-6, np.exp(-10)):
            model = CoreScaling(.002, 3., h)
            speed = 3*q**(-.5-h)
            expected = .002*.03*q**(.5+h/2)
            np.testing.assert_allclose(carrier_power_length(model, speed, .03), expected, rtol=6e-15)

    def test_screen_crossing_against_root_search(self):
        model = CoreScaling(.003, 2., 1e-6)
        for beta in (1., .1, .001):
            root = brentq(lambda logu: self.gas.mean_free_path_m/(beta*model.radius(np.exp(logu)))-.01, -20, 20)
            self.assertAlmostEqual(kn_crossing_speed(model,beta,self.gas.mean_free_path_m,.01)/np.exp(root),1,places=11)
            rate = rate_crossing_speed(model,beta,self.gas.collision_time_proxy,.1)
            d = gas_diagnostics(rate,beta*model.radius(rate),self.gas)
            self.assertAlmostEqual(float(d["collision_rate_proxy"]),.1,places=13)

    def test_independent_timescale_ratios(self):
        u = np.geomspace(.1, 1e4, 19)
        length = .1*self.model.radius(u)
        d = gas_diagnostics(u,length,self.gas)
        np.testing.assert_allclose(d["sound_crossing_time_s"]/d["advection_time_proxy_s"],d["mach"])
        np.testing.assert_allclose(d["diffusion_time_proxy_s"]/d["advection_time_proxy_s"],d["reynolds_proxy"])
        np.testing.assert_allclose(d["collision_rate_proxy"],d["knudsen_proxy"]*u/self.gas.mean_thermal_speed)

    def test_gradient_assumption_reverses_screen_order(self):
        mach_speed = .3*self.gas.sound_speed
        self.assertGreater(kn_crossing_speed(self.model,1,self.gas.mean_free_path_m,.01),mach_speed)
        self.assertLess(kn_crossing_speed(self.model,.1,self.gas.mean_free_path_m,.01),mach_speed)
        self.assertLess(kn_crossing_speed(self.model,.001,self.gas.mean_free_path_m,.01),self.model.u0_m_s)

    def test_volume_count_and_unit_invariance(self):
        d = 3.1e-10
        self.assertAlmostEqual(float(averaging_cube_count(10*d,d,.1)),1)
        self.assertAlmostEqual(float(averaging_cube_count(100*d,d,.1)),1000)
        np.testing.assert_allclose(averaging_cube_count(100*d,d,.03),
                                   averaging_cube_count(100*d*1e9,d*1e9,.03))

    def test_inverse_limit_constant_reynolds_and_carrier_ratio(self):
        u = np.geomspace(.01, 3e7, 47)
        d = gas_diagnostics(u,self.model.radius(u),self.gas)
        np.testing.assert_allclose(d["reynolds_proxy"],self.model.r0_m*self.model.u0_m_s/self.gas.kinematic_viscosity)
        np.testing.assert_allclose(carrier_power_length(self.model,u,.001)/self.model.radius(u),.001)

    def test_invalid_physical_inputs(self):
        for value in (0., -1., np.nan, np.inf):
            with self.assertRaises(ValueError):
                gas_diagnostics(1.,value,self.gas)
            with self.assertRaises(ValueError):
                carrier_power_length(self.model,1.,value)
        with self.assertRaises(ValueError):
            averaging_cube_count(1.,.01,1.)


if __name__ == "__main__":
    unittest.main()
