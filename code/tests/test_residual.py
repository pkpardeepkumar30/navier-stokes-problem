"""Residual checks using independent integration and physical finite differences."""
import unittest
import numpy as np
import mpmath as mp
from nscomp.energy import CoreEnergyScaling, gaussian_swirl
from nscomp.residual import HeatExterior, gaussian_pressure, gaussian_terms


class ResidualChecks(unittest.TestCase):
    def test_profile_against_high_precision_integral(self):
        mp.mp.dps = 45
        for h in (1e-6, 0.005):
            model = HeatExterior(h=h, order=128)
            mh = mp.mpf(str(h))
            for zeta in (0, 0.1, 1.0):
                mz = mp.mpf(str(zeta))
                for n in (0, 1, 2):
                    integral = mp.quad(lambda v: mp.exp(-v)*v**(mh+n)*(1+mz*v)**(-mh-n),
                                       [0, 1, mp.inf])
                    factor = (-1)**n*(mp.rf(mh, n) if n else 1)/mp.gamma(1+mh)
                    exact = float(factor*integral)
                    self.assertAlmostEqual(float(model.profile(zeta, n))/exact, 1, delta=2e-10)

    def test_radial_and_time_derivatives_independently(self):
        model = HeatExterior(order=128)
        r, tau = 2.1, 0.12
        dr, dt = 2e-4*r, 1e-3*tau
        exact = model.derivatives(r, tau)
        # A fourth-order first derivative, with physical time t increasing
        # opposite to tau, independent of the analytic chain-rule operators.
        ur = (model.velocity(r-2*dr,tau)-8*model.velocity(r-dr,tau)
              +8*model.velocity(r+dr,tau)-model.velocity(r+2*dr,tau))/(12*dr)
        ut = -(model.velocity(r,tau-2*dt)-8*model.velocity(r,tau-dt)
               +8*model.velocity(r,tau+dt)-model.velocity(r,tau+2*dt))/(12*dt)
        self.assertAlmostEqual(float(ur/exact["radial"]), 1, delta=1e-10)
        self.assertAlmostEqual(float(ut/exact["time"]), 1, delta=2e-5)

    def test_pressure_integral_derivative(self):
        model = HeatExterior(h=1e-6, order=128)
        r, tau = 1.7, 0.2
        step = 2e-4*r
        gradient = (model.pressure(r-2*step,tau)-8*model.pressure(r-step,tau)
                    +8*model.pressure(r+step,tau)-model.pressure(r+2*step,tau))/(12*step)
        expected = model.velocity(r,tau)**2/r
        self.assertAlmostEqual(float(gradient/expected), 1, delta=2e-10)

    def test_heat_residual_is_small_component_by_component(self):
        model = HeatExterior(order=128)
        for tau in (1, 0.01, 1e-6):
            radii = np.sqrt(2*tau*np.array([4,16,64,256]))
            terms = model.terms(radii, tau)
            self.assertLess(float(np.max(terms["component_normalized"])), 2e-10)
            self.assertGreater(float(np.min(abs(terms["time"][:,1]))), 0)

    def test_h_zero_limit_is_steady_inverse_radius(self):
        model = HeatExterior(h=0, amplitude=2)
        r = np.array([0.5, 1, 3])
        np.testing.assert_allclose(model.velocity(r,0.1), 2*np.sqrt(2)/r, rtol=1e-14)
        np.testing.assert_allclose(model.pressure(r,0.1), -4/(r*r), rtol=2e-14)
        np.testing.assert_allclose(model.terms(r,0.1)["residual"], 0, atol=1e-14)

    def test_gaussian_time_laplacian_and_pressure(self):
        model = CoreEnergyScaling(h=0.003)
        r, z, tau = 0.6, 0.4, 0.3
        def u(r,z,tau):
            a,b,U=model.widths_and_speed(tau)
            return gaussian_swirl(r,z,a,b,U)
        dr,dz,dt=1e-4,1e-4,1e-5
        value=u(r,z,tau)
        ur=(u(r+dr,z,tau)-u(r-dr,z,tau))/(2*dr)
        urr=(u(r+dr,z,tau)-2*value+u(r-dr,z,tau))/dr**2
        uzz=(u(r,z+dz,tau)-2*value+u(r,z-dz,tau))/dz**2
        lap=urr+ur/r-value/r**2+uzz
        ut=-(u(r,z,tau+dt)-u(r,z,tau-dt))/(2*dt)
        dpz=(gaussian_pressure(r,z+dz,tau,model)-gaussian_pressure(r,z-dz,tau,model))/(2*dz)
        terms=gaussian_terms(r,z,tau,model)
        self.assertAlmostEqual(float(ut/terms["time"][1]),1,delta=1e-8)
        self.assertAlmostEqual(float(-lap/terms["minus_viscosity"][1]),1,delta=2e-7)
        self.assertAlmostEqual(float(dpz/terms["pressure"][2]),1,delta=1e-8)

    def test_gaussian_force_growth_at_fixed_similarity_marker(self):
        model=CoreEnergyScaling()
        magnitudes=[]
        times=np.array([1,1e-2,1e-4])
        for tau in times:
            a,b,U=model.widths_and_speed(tau)
            magnitudes.append(np.linalg.norm(gaussian_terms(a,0.5*b,tau,model)["residual"]))
        np.testing.assert_allclose(np.array(magnitudes)/magnitudes[0],times**(-1.5),rtol=2e-14)

    def test_invalid_exterior_domain(self):
        for r in (0,-1,np.nan,np.inf):
            with self.assertRaises(ValueError):
                HeatExterior().velocity(r,1)
        with self.assertRaises(ValueError):
            HeatExterior(h=-1)


if __name__ == "__main__":
    unittest.main()
