"""Independent identities and resolved-grid checks for D4's momentum flux."""
import unittest
import numpy as np
import sympy as sp
from nscomp.stress import (periodic_grid, derivative, two_mode_velocity, two_mode_target,
    covariance, divergence, envelope, enveloped_velocity, envelope_covariance,
    envelope_mean_force, averaged_transport, force_from_mean_covariance, project_mean_force)


class StressChecks(unittest.TestCase):
    def test_sine_averages(self):
        x, _, _ = periodic_grid(128)
        self.assertAlmostEqual(float(np.sin(8*x).mean()), 0, delta=1e-14)
        self.assertAlmostEqual(float((np.sin(8*x)**2).mean()), 0.5, delta=1e-14)

    def test_prescribed_covariance_and_energy(self):
        _, x, y = periodic_grid(128)
        for lp, lm in ((3,1), (0,2), (4,0)):
            w = two_mode_velocity(x,y,lp,lm)
            q = covariance(w, (0,1))
            np.testing.assert_allclose(w.mean(axis=(0,1)), 0, atol=2e-14)
            np.testing.assert_allclose(q, two_mode_target(lp,lm), atol=2e-13)
            np.testing.assert_allclose(np.linalg.eigvalsh(q), sorted((lp,lm)), atol=3e-13)
            self.assertAlmostEqual(float(np.trace(q)/2), float(np.mean(np.sum(w*w,axis=-1))/2), delta=2e-13)

    def test_exact_symbolic_streamfunction(self):
        x,y,k,m = sp.symbols("x y k m", real=True, nonzero=True)
        a = sp.Function("a")(y)
        psi = a*sp.cos(k*(x+m*y))/k
        wx, wy = sp.diff(psi,y), -sp.diff(psi,x)
        self.assertEqual(sp.simplify(sp.diff(wx,x)+sp.diff(wy,y)),0)
        # Dropping the envelope correction creates a nonzero divergence.
        incomplete = -m*a*sp.sin(k*(x+m*y))
        defect = sp.simplify(sp.diff(incomplete,x)+sp.diff(wy,y))
        self.assertEqual(sp.simplify(defect-sp.diff(a,y)*sp.sin(k*(x+m*y))),0)

    def test_resolved_fourier_divergence(self):
        _,x,y = periodic_grid(128)
        for w in (two_mode_velocity(x,y), enveloped_velocity(x,y)):
            # Normalize by carrier wavenumber times the sampled speed scale.
            relative = float(np.max(abs(divergence(w))))/(9*np.max(np.linalg.norm(w,axis=-1)))
            self.assertLess(relative, 1e-13)

    def test_spatial_covariance_and_positive_eigenvalues(self):
        points,x,y = periodic_grid(128)
        w = enveloped_velocity(x,y)
        exact,_ = envelope_covariance(points)
        q = covariance(w,0)
        np.testing.assert_allclose(w.mean(axis=0),0,atol=2e-14)
        np.testing.assert_allclose(q,exact,atol=2e-14)
        self.assertGreaterEqual(float(np.linalg.eigvalsh(q).min()), -1e-14)

    def test_mean_force_from_independent_transport(self):
        points,x,y = periodic_grid(128)
        w = enveloped_velocity(x,y)
        exact = envelope_mean_force(points)
        np.testing.assert_allclose(-averaged_transport(w),exact,atol=8e-13)
        np.testing.assert_allclose(force_from_mean_covariance(covariance(w,0)),exact,atol=8e-13)
        np.testing.assert_allclose(exact.mean(axis=0),0,atol=1e-15)

    def test_pressure_removes_only_gradient(self):
        points,_,_ = periodic_grid(128)
        a,ap,app = envelope(points)
        force = envelope_mean_force(points)
        pressure = -a*a/2
        pressure_force = np.stack((np.zeros_like(a),-derivative(pressure,0)),axis=-1)
        projected = project_mean_force(force)
        np.testing.assert_allclose(force+pressure_force,projected,atol=3e-14)
        np.testing.assert_allclose(-derivative(projected[:,0],0),-(ap*ap+a*app),atol=2e-14)
        self.assertGreater(float(np.max(abs(projected[:,0]))), 0.4)

    def test_finite_frequency_correction(self):
        points,_,_ = periodic_grid(256)
        corrections=[]
        for k in (2,4,8,16,32):
            q,leading = envelope_covariance(points,k=k)
            corrections.append(float(np.max(abs(q-leading))))
        np.testing.assert_allclose(np.array(corrections[:-1])/corrections[1:],4,rtol=2e-12)

    def test_nyquist_sampling_can_lose_stress(self):
        _,x,y = periodic_grid(16)
        q = covariance(two_mode_velocity(x,y),(0,1))
        self.assertGreater(float(np.linalg.norm(q-two_mode_target())),2.9)
        # N=20 captures the mean accurately but cannot resolve all quadratic modes.
        _,x,y = periodic_grid(20)
        np.testing.assert_allclose(covariance(two_mode_velocity(x,y),(0,1)),two_mode_target(),atol=1e-14)

    def test_invalid_periodic_carriers(self):
        for kwargs in ({"n":1.5},{"lambda_plus":-1}):
            with self.assertRaises(ValueError): two_mode_velocity(0.,0.,**kwargs)
        with self.assertRaises(ValueError): enveloped_velocity(0.,0.,slope=0.5)


if __name__ == "__main__":
    unittest.main()
