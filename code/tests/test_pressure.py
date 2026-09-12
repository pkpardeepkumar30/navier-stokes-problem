"""Independent checks for the source's ideal-schedule axis-pressure integral."""
import math
import unittest
import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp
from nscomp.pressure import Schedule, pressure_trace, smooth_step, step_derivative
from nscomp.inner import AxisData, inner_series


class PressureChecks(unittest.TestCase):
    def test_smooth_step_symmetry_and_derivative(self):
        x=np.linspace(.05,.95,41)
        np.testing.assert_allclose(smooth_step(x)+smooth_step(1-x),1,atol=2e-15)
        np.testing.assert_allclose((smooth_step(x+1e-6)-smooth_step(x-1e-6))/(2e-6),
                                   step_derivative(x),atol=2e-8,rtol=2e-7)
        self.assertAlmostEqual(quad(lambda x:float(smooth_step(x)),0,1,
                                    epsabs=1e-13,epsrel=1e-13)[0],.5,places=13)

    def test_independent_adaptive_pressure_integral(self):
        trace=pressure_trace()
        total=5.
        logc=0.
        for segment in trace.segments[:5]:
            slope=segment.ell_start-.5
            if segment.kind=="slope":
                delta=segment.ell_end-segment.ell_start
                def relative(y):
                    integral=quad(lambda v:float(smooth_step(v)),0,y,
                                  epsabs=2e-13,epsrel=2e-13)[0]
                    return slope*y+delta*integral
                local=quad(lambda y:math.exp(2*relative(y)),0,1,
                           epsabs=2e-12,epsrel=2e-12)[0]
                total+=math.exp(2*logc)*local
                logc+=relative(1)
            else:
                # Integrate a transformed semi-infinite exponential; account
                # for its exact removed tail at the declared finite endpoint.
                decay=-2*slope
                integral=quad(lambda s:math.exp(-s),0,np.inf,
                              epsabs=2e-13,epsrel=2e-13)[0]/decay
                total+=math.exp(2*logc)*integral*(-math.expm1(-decay*segment.length))
                logc+=slope*segment.length
        omitted=math.exp(2*logc)/(1+2*trace.schedule.intermediate_lambda)
        self.assertLess(omitted,1e-113)
        self.assertLess(abs(float(trace.normalized(0))+total/2),2e-13)

    def test_quadrature_refinement_of_pressure_and_terminal_datum(self):
        reference=pressure_trace(order=256)
        pressure_errors=[]
        q_errors=[]
        for order in (8,16,32,64,128):
            t=pressure_trace(order=order)
            pressure_errors.append(abs(float(t.normalized(0)-reference.normalized(0))))
            q_errors.append(abs(t.Q_target/reference.Q_target-1))
        self.assertTrue(all(b<a/10 for a,b in zip(pressure_errors[:3],pressure_errors[1:4])))
        self.assertLess(pressure_errors[-1],1e-13)
        self.assertLess(q_errors[-1],1e-13)
        with mp.workdps(45):
            h=mp.mpf(str(reference.schedule.h))
            rho=h*mp.mpf(str(reference.schedule.terminal_c))
            def integrand(u):
                if u==0 or u==1:
                    return mp.mpf(0)
                s=1/(1+mp.exp(1/u**2-1/(1-u)**2))
                ds=s*(1-s)*(2/u**3+2/(1-u)**3)
                return mp.exp((1-h)*(1+2*u))*ds
            exact=rho/(1-rho)*mp.quad(integrand,[0,.25,.5,.75,1])
            self.assertLess(abs(reference.Q_target/float(exact)-1),2e-13)

    def test_terminal_schedule_against_independent_ode(self):
        trace=pressure_trace()
        lam,h=trace.schedule.intermediate_lambda,trace.schedule.h
        Q=(lam-h)/(1-lam)
        for ell0,ell1 in ((-lam,-1),(-1,-h)):
            def rhs(y,q):
                ell=ell0+(ell1-ell0)*float(smooth_step(y))
                return [-(1+ell)*q[0]-ell-h]
            result=solve_ivp(rhs,(0,1),[Q],method="DOP853",rtol=2e-13,atol=2e-14)
            self.assertTrue(result.success)
            Q=result.y[0,-1]
            if ell1==-1:
                Q+=(1-h)*4*math.log(1/h)
        self.assertLess(abs(Q/trace.Q_before_wait-1),2e-12)
        self.assertLess(abs(Q*math.exp(-(1-h)*trace.waiting_length)/trace.Q_target-1),2e-12)
        rho=trace.schedule.terminal_c*h
        def terminal_rhs(y,q):
            f=1-rho*(1-float(smooth_step((y-1)/2)))
            ratio=rho*float(step_derivative((y-1)/2))/(2*f)
            return [-(1-h+ratio)*q[0]-ratio]
        result=solve_ivp(terminal_rhs,(0,3),[trace.Q_target],method="DOP853",
                         rtol=2e-12,atol=1e-23,max_step=.08)
        self.assertLess(abs(result.y[0,-1]),1e-22)

    def test_pressure_sign_bound_and_eta_derivatives(self):
        trace=pressure_trace()
        eta=np.linspace(0,1,31)
        np.testing.assert_allclose(trace.normalized(eta),trace.normalized(-eta),atol=0)
        self.assertTrue(np.all(trace.normalized(eta)<-2.5/(1+eta**2)**2))
        self.assertTrue(np.all(trace.normalized(eta[1:],1)>0))
        e,step=.3,2e-4
        f=lambda x:float(trace.normalized(x))
        first=(f(e-2*step)-8*f(e-step)+8*f(e+step)-f(e+2*step))/(12*step)
        second=(-f(e+2*step)+16*f(e+step)-30*f(e)+16*f(e-step)-f(e-2*step))/(12*step**2)
        self.assertLess(abs(first-trace.normalized(e,1)),2e-11)
        self.assertLess(abs(second-trace.normalized(e,2)),1e-7)

    def test_eta_coefficients_against_direct_high_precision_differentiation(self):
        trace=pressure_trace()
        with mp.workdps(55):
            e=mp.mpf(".2")
            direct=lambda x:-mp.fsum(mp.exp(mp.mpf(str(r["log_weight"])))
                    *(1+x*x)**(-2*mp.mpf(str(r["theta"]))) for r in trace.components)/2
            coefficients=trace.eta_coefficients(.2,4,55,normalized=True)
            for n in range(5):
                expected=mp.diff(direct,e,n)/mp.factorial(n)
                self.assertLess(abs(coefficients[n]-expected),mp.mpf("1e-48"))

    def test_source_trace_enters_inner_equations(self):
        trace=pressure_trace()
        data=AxisData(h=trace.schedule.h,regularization_sigma=1e-4,
                      radial_scale_Lambda=100*math.exp(2*trace.schedule.log_P))
        series=inner_series(data,.2,12,70,axis_pressure=trace)
        axis=series.evaluate(0)
        with mp.workdps(70):
            expected=trace.eta_coefficients(.2,0,70)[0]
            self.assertLess(abs((axis["Pi"]-expected)/expected),mp.mpf("1e-60"))
        results=[series.evaluate(1,order=n) for n in (4,8,12)]
        errors=[max(float(v[k]) for k in ("angular_relative","axial_relative","pressure_relative"))
                for v in results]
        self.assertTrue(all(b<a/100 for a,b in zip(errors[:-1],errors[1:])),errors)
        self.assertLess(errors[-1],1e-10)

    def test_axis_condition_bound_and_scale_independence(self):
        trace=pressure_trace()
        bound=trace.axis_condition_bound()
        self.assertTrue(bound["condition_sufficient"])
        self.assertLess(float(bound["eta_band_bound"]),4e-15)
        self.assertGreater(float(bound["chi_lower"]),.9999)
        # Independent physical-radius integration on the finite initial branch.
        # E/P=(X/X_R)^.1*f, including its integrable axis endpoint.
        eta=.4
        for radius in (.01,1,100):
            value=quad(lambda x:(x/radius)**.2/(2*x*(1+eta**2)**2),
                       0,radius,epsabs=1e-11,epsrel=1e-11)[0]
            self.assertAlmostEqual(value,2.5/(1+eta**2)**2,places=11)

    def test_invalid_schedule_and_inputs(self):
        for value in (0,-1,float("nan")):
            with self.assertRaises(ValueError):
                Schedule(h=value)
        with self.assertRaises(ValueError):
            Schedule(h=1e-3)
        with self.assertRaises(ValueError):
            pressure_trace(order=2)
        with self.assertRaises(ValueError):
            pressure_trace().normalized(2)


if __name__=="__main__":
    unittest.main()
