"""Validate EOS, conservative equations, nonlinear radial solve, and local scale audit."""
import unittest
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from nscomp.compressible import (IsentropicVortex,integrate_density,density_deficit_mach,local_gas_audit)
from nscomp.regimes import GasReference


class CompressibleChecks(unittest.TestCase):
    def test_independent_enthalpy_integral(self):
        v = IsentropicVortex(.7)
        for r in (0,.2,1,2,4):
            potential = quad(lambda x: v.peak_mach**2*x*np.exp(1-x*x),r,np.inf,
                             epsabs=1e-13,epsrel=1e-13)[0]
            theta = 1-(v.gamma-1)*potential
            self.assertAlmostEqual(theta,float(v.state(r)["temperature_ratio"]),places=13)

    def test_eos_entropy_and_peak_normalization(self):
        v = IsentropicVortex(.8)
        s = v.state(np.linspace(0,6,1001))
        np.testing.assert_allclose(s["pressure_ratio"],s["density_ratio"]*s["temperature_ratio"],rtol=6e-15)
        np.testing.assert_allclose(s["pressure_ratio"]/s["density_ratio"]**v.gamma,1,rtol=6e-15)
        self.assertEqual(float(v.state(0)["swirl_over_exterior_sound_speed"]),0)
        self.assertAlmostEqual(float(v.state(1)["swirl_over_exterior_sound_speed"]),.8)

    def test_all_four_conservative_euler_equations(self):
        v = IsentropicVortex(.7)
        h = 2e-4
        for x,y in ((0,0),(.3,.4),(1,.2),(-.7,1.2),(2,2)):
            def derivative(axis):
                values = []
                for shift in (-2,-1,1,2):
                    xx,yy = (x+shift*h,y) if axis==0 else (x,y+shift*h)
                    values.append(v.conservative_flux(xx,yy)[axis])
                return (values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
            np.testing.assert_allclose(derivative(0)+derivative(1),0,atol=3e-12)

    def test_nonlinear_radial_solver_converges_fourth_order(self):
        v = IsentropicVortex(.6)
        errors = []
        for n in (32,64,128,256):
            r,rho = integrate_density(v,4,n)
            errors.append(np.max(np.abs(rho-v.state(r)["density_ratio"])))
        self.assertLess(errors[-1],2e-10)
        self.assertTrue(all(13 < x/y < 19 for x,y in zip(errors[:-1],errors[1:])),errors)

    def test_outer_boundary_error_is_separate(self):
        v = IsentropicVortex(.6)
        for outer in (2,3,4):
            r,rho = integrate_density(v,outer,1024,boundary="far_field")
            expected = (1-v.depression*(1-np.exp(-outer*outer)))**(1/(v.gamma-1))
            self.assertAlmostEqual(rho[0],expected,delta=2e-12)
            self.assertGreater(rho[0],float(v.state(0)["density_ratio"]))

    def test_small_mach_density_and_pressure_orders(self):
        # Stable deficit formulas avoid subtracting density close to one.
        for m in (.01,.005,.0025):
            v = IsentropicVortex(m)
            s = v.state(0)
            expected = .5*np.e*m*m
            self.assertLess(abs(float(s["density_deficit"])/expected-1),5e-5)
        errors = []
        for m in (.04,.02,.01):
            v = IsentropicVortex(m)
            s = v.state(0)
            errors.append(float(s["pressure_ratio"]-s["incompressible_pressure_ratio"]))
        self.assertTrue(all(15.9 < x/y < 16.1 for x,y in zip(errors[:-1],errors[1:])))

    def test_density_marker_and_local_gradient(self):
        m = float(density_deficit_mach(.05))
        root = brentq(lambda x: float(IsentropicVortex(x).state(0)["density_deficit"])-.05,.01,.5)
        self.assertAlmostEqual(m,root,places=11)
        gas = GasReference(296.15,101300,287.05,1.4,1.83245e-5,6.73e-8)
        v = IsentropicVortex(.6)
        r,h = .7,1e-5
        a = .001
        audit = local_gas_audit(v,r,a,gas)
        derivative = (np.log(v.state(r+h)["density_ratio"])-np.log(v.state(r-h)["density_ratio"]))/(2*h*a)
        self.assertAlmostEqual(float(audit["kn_density"])/(float(audit["mean_free_path_m"])*derivative),1,places=9)
        swirl_plus=v.state(r+h)["swirl_over_exterior_sound_speed"]
        swirl_minus=v.state(r-h)["swirl_over_exterior_sound_speed"]
        numerical_shear=gas.sound_speed/a*abs((swirl_plus-swirl_minus)/(2*h)-v.state(r)["swirl_over_exterior_sound_speed"]/r)
        self.assertAlmostEqual(float(audit["shear_rate_s_inv"])/numerical_shear,1,places=9)
        outer = local_gas_audit(v,8,a,gas)
        self.assertAlmostEqual(float(outer["mean_free_path_m"])/gas.mean_free_path_m,1,places=14)

    def test_rest_limit_and_positivity_rejection(self):
        v = IsentropicVortex(0)
        np.testing.assert_array_equal(v.state([0,1,3])["density_ratio"],1)
        r,rho = integrate_density(v,4,16)
        np.testing.assert_array_equal(rho,1)
        for m in (-.1,np.nan,1.5):
            with self.assertRaises(ValueError):
                IsentropicVortex(m)
        with self.assertRaises(ValueError):
            IsentropicVortex(.3,gamma=1.)


if __name__=="__main__":
    unittest.main()
