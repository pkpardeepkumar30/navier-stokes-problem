"""Validate D5 against energy identities, primitive momentum, and an ODE solver."""
import unittest
import numpy as np
from scipy.integrate import quad
from nscomp.shear import KelvinWave, integrate_amplitude


class ShearChecks(unittest.TestCase):
    def test_damping_integral_against_quadrature(self):
        wave=KelvinWave()
        for t in (0.1,4,12,20):
            integral=quad(lambda s:wave.kx**2+(wave.ky0-wave.shear*wave.kx*s)**2,0,t)[0]
            self.assertAlmostEqual(float(wave.state(t)['damping_integral'])/integral,1,delta=5e-15)

    def test_heat_limit(self):
        wave=KelvinWave(shear=0)
        state=wave.state(np.linspace(0,20,101))
        np.testing.assert_allclose(state['speed'],state['heat_reference_speed'],rtol=5e-15)
        np.testing.assert_allclose(state['production'],0,atol=0)

    def test_inviscid_vorticity_and_peak(self):
        wave=KelvinWave(viscosity=0)
        state=wave.state(np.array([0,4,8]))
        np.testing.assert_allclose(state['vorticity_amplitude'],np.sqrt(17),rtol=1e-15)
        np.testing.assert_allclose(state['speed'],[1,np.sqrt(17),1],rtol=1e-15)

    def test_incompressibility_and_phase_average(self):
        wave=KelvinWave()
        x=2*np.pi*np.arange(128)/128
        for t in (0.2,3.9,12.1):
            s=wave.state(t)
            b=s['velocity_amplitude']
            self.assertAlmostEqual(wave.kx*b[0]+s['ky']*b[1],0,delta=1e-15)
            values=wave.field(x,0.371,t)
            np.testing.assert_allclose(values.mean(axis=0),0,atol=3e-16)
            self.assertAlmostEqual(float((values[:,0]*values[:,1]).mean()),float(s['qxy']),delta=2e-15)
            self.assertAlmostEqual(float(np.mean(np.sum(values*values,axis=1))/2),float(s['energy']),delta=2e-15)

    def test_energy_budget_by_time_differences(self):
        wave=KelvinWave()
        dt=2e-4
        for t in (0.4,3.5,4.4,12):
            E=lambda t:float(wave.state(t)['energy'])
            derivative=(E(t-2*dt)-8*E(t-dt)+8*E(t+dt)-E(t+2*dt))/(12*dt)
            self.assertAlmostEqual(derivative,float(wave.state(t)['energy_rate']),delta=2e-11)

    def test_primitive_momentum_with_physical_finite_differences(self):
        wave=KelvinWave()
        x,y=0.31,0.47
        step=2e-4
        def first(f,s): return (f(s-2*step)-8*f(s-step)+8*f(s+step)-f(s+2*step))/(12*step)
        def second(f,s): return (-f(s+2*step)+16*f(s+step)-30*f(s)+16*f(s-step)-f(s-2*step))/(12*step**2)
        for t in (0.7,3.9,12.3):
            u=wave.field(x,y,t)
            dt=first(lambda s:wave.field(x,y,s),t)
            dx=first(lambda s:wave.field(s,y,t),x)
            dy=first(lambda s:wave.field(x,s,t),y)
            lap=second(lambda s:wave.field(s,y,t),x)+second(lambda s:wave.field(x,s,t),y)
            gradp=np.array([first(lambda s:wave.pressure(s,y,t),x),first(lambda s:wave.pressure(x,s,t),y)])
            transport=(wave.shear*y+u[0])*dx+u[1]*dy+np.array([wave.shear*u[1],0.])
            residual=dt+transport-wave.viscosity*lap+gradp
            scale=np.linalg.norm(dt)+np.linalg.norm(transport)+np.linalg.norm(wave.viscosity*lap)+np.linalg.norm(gradp)
            self.assertLess(float(np.linalg.norm(residual)/scale),2e-8)

    def test_rk4_refinement_and_divergence_constraint(self):
        wave=KelvinWave()
        errors=[]
        for step in (.2,.1,.05):
            times,values=integrate_amplitude(wave,20,step)
            exact=wave.state(times)
            errors.append(float(np.max(np.linalg.norm(values-exact['velocity_amplitude'],axis=1))))
        self.assertGreater(errors[0]/errors[1],12)
        self.assertGreater(errors[1]/errors[2],12)
        self.assertLess(errors[-1],1e-5)
        self.assertLess(float(np.max(abs(wave.kx*values[:,0]+exact['ky']*values[:,1]))),1e-5)

    def test_wavelength_and_late_decay(self):
        wave=KelvinWave()
        s=wave.state(np.array([0,4,20]))
        self.assertGreater(s['wavelength'][1],s['wavelength'][0])
        self.assertLess(s['wavelength'][2],s['wavelength'][0])
        self.assertLess(s['speed'][2],1e-3)
        self.assertLess(s['qxy'][0],0)
        self.assertGreater(s['qxy'][2],0)


if __name__=='__main__':
    unittest.main()
