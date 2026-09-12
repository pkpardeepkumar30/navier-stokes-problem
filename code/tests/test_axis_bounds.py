"""Amplitude-envelope proof arithmetic and independent inner-exit diagnostics."""
from dataclasses import replace
from fractions import Fraction
import math
import unittest
import mpmath as mp
from nscomp.axis_bounds import (complex_amplitude_bound,axis_primitive,
    logarithmic_axis_value,feature_scales,exit_diagnostics)
from nscomp.inner import AxisData,inner_series,_Jet
from nscomp.pressure import pressure_trace


def data_for(ratio=1e24):
    trace=pressure_trace()
    data=AxisData(h=trace.schedule.h,regularization_sigma=1e-4,
                  radial_scale_Lambda=ratio*math.exp(2*trace.schedule.log_P))
    return replace(data,log_amplitude_C=complex_amplitude_bound(data)["log_C_upper"])


class AxisBoundsChecks(unittest.TestCase):
    def test_exact_rectangle_bound_and_upward_rounding(self):
        data=data_for()
        b=complex_amplitude_bound(data)
        self.assertGreaterEqual(Fraction(b["log_C_upper"]),b["required_log_C"])
        self.assertLess(b["H_derivative_bound"],17)
        self.assertLessEqual(17*b["half_width"],Fraction(str(data.regularization_sigma))/4)
        with mp.workdps(70):
            width=mp.mpf(b["half_width"].numerator)/b["half_width"].denominator
            bound=mp.mpf(b["zeta_modulus_bound"].numerator)/b["zeta_modulus_bound"].denominator
            for x in (-1-width,-.00445,0,1+width):
                for y in (-width,0,width):
                    z=mp.mpc(x,y)
                    H=(mp.mpf(".5")-mp.mpf(str(data.h)))*z+(1-z*z)*(4*z+mp.mpf(".02"))
                    zeta=-(1-2*mp.mpf(str(data.h))*z*z)*H/(H*H+mp.mpf("1e-8"))
                    self.assertLess(abs(zeta),bound)

    def test_primitive_against_independent_split_quadrature(self):
        data=data_for()
        with mp.workdps(65):
            feature=feature_scales(data,pressure_trace(),65)
            center,width=feature["center"],feature["chi_width"]
            h=mp.mpf(str(data.h))
            def integrand(e):
                H=(mp.mpf(".5")-h)*e+(1-e*e)*(4*e+mp.mpf(".02"))
                return -(1-2*h*e*e)*H/(H*H+mp.mpf("1e-8"))
            for eta in (mp.mpf("-.6"),mp.mpf(".2")):
                points=[eta,mp.mpf(0)]
                points.extend(center+k*width for k in (-32,-8,-2,0,2,8,32)
                              if min(eta,0)<center+k*width<max(eta,0))
                integral=mp.quad(integrand,sorted(points))*(1 if eta>0 else -1)
                self.assertLess(abs(axis_primitive(data,eta,65)-integral),mp.mpf("1e-50"))

    def test_logarithmic_amplitude_retains_precision(self):
        data=data_for()
        with mp.workdps(110):
            a,ga=logarithmic_axis_value(data,-.6,60)
            b,gb=logarithmic_axis_value(data,-.6,90)
            self.assertLess(abs((ga-gb)/gb),mp.mpf("1e-60"))
            self.assertLess(a,-1)
            self.assertTrue(mp.isfinite(ga) and ga>0)

    def test_log_input_matches_moderate_original_amplitude(self):
        old=AxisData()
        with mp.workdps(75):
            new=replace(old,log_amplitude_C=mp.nstr(mp.log(10000),70))
            a=inner_series(old,.05,8,65).evaluate(.5)
            b=inner_series(new,.05,8,65).evaluate(.5)
            for key in ("Phi","U","F"):
                self.assertLess(abs((a[key]-b[key])/a[key]),mp.mpf("1e-58"))

    def test_centered_jet_values_against_polynomial_derivatives(self):
        with mp.workdps(50):
            coefficients=[mp.mpf(n*n+1)/7 for n in range(9)]
            jet=_Jet(coefficients)
            polynomial=lambda x:mp.fsum(c*x**n for n,c in enumerate(coefficients))
            for n in range(10):
                self.assertLess(abs(jet.value(0,n)-mp.diff(polynomial,0,n)),mp.mpf("1e-42"))

    def test_shear_diagnostics_against_independent_radial_difference(self):
        data=data_for()
        s=inner_series(data,.2,18,80,axis_pressure=pressure_trace())
        with mp.workdps(80):
            Y,step=mp.mpf(4),mp.mpf("1e-4")
            def difference(key):
                v=[s.evaluate(Y+k*step)[key] for k in (-2,-1,1,2)]
                return (v[0]-8*v[1]+8*v[2]-v[3])/(12*step)
            v=s.evaluate(Y)
            self.assertLess(abs(-2*Y*difference("Phi")/v["Phi"]-v["p1"]),mp.mpf("1e-17"))
            expected_ns=-2*mp.mpf(str(data.radial_scale_Lambda))*difference("U")
            self.assertLess(abs((expected_ns-v["n_s"])/v["n_s"]),mp.mpf("1e-17"))

    def test_insufficient_scale_failure_survives_refinement(self):
        data=data_for(1e16)
        with mp.workdps(95):
            f=feature_scales(data,pressure_trace(),95)
            e=f["center"]+f["pressure_width"]/2
            a=exit_diagnostics(inner_series(data,e,12,65,axis_pressure=pressure_trace()))
            b=exit_diagnostics(inner_series(data,e,18,90,axis_pressure=pressure_trace()))
            self.assertLess(a["scaled_p1"],-1e5)
            self.assertIsNone(a["exit_lower"])
            self.assertLess(abs((a["scaled_p1"]-b["scaled_p1"])/b["scaled_p1"]),mp.mpf("1e-28"))

    def test_repaired_critical_band_and_pressure_sensitivity(self):
        data=data_for()
        with mp.workdps(90):
            f=feature_scales(data,pressure_trace(),90)
            for c in (0,mp.mpf(2)/3,20):
                e=f["center"]+c*f["pressure_width"]
                s=inner_series(data,e,12,80,axis_pressure=pressure_trace())
                d=exit_diagnostics(s)
                self.assertGreater(d["scaled_p1"],5.99)
                self.assertGreater(d["Sq_margin"],.35)
                self.assertGreater(d["exit_lower"],mp.mpf("2.2"))
            other=exit_diagnostics(inner_series(data,e,12,80,axis_pressure=pressure_trace(order=256)))
            self.assertLess(abs(other["Sq_margin"]-d["Sq_margin"]),mp.mpf("1e-12"))

    def test_global_and_transition_samples(self):
        data=data_for()
        with mp.workdps(85):
            f=feature_scales(data,pressure_trace(),85)
            centers=[mp.mpf(-1),mp.mpf(0),mp.mpf(1)]
            centers.extend(f["center"]+x*f["chi_width"] for x in (-12,-2,2,12))
            for e in centers:
                s=inner_series(data,e,18,75,axis_pressure=pressure_trace())
                d=exit_diagnostics(s)
                self.assertGreater(d["p1"],0)
                self.assertGreater(d["exit_lower"],mp.mpf("2.2"))
                self.assertGreater(s.evaluate("4.1")["Phi"],0)

    def test_invalid_logarithm_and_bound(self):
        for value in ("nan","-1","inf"):
            with self.assertRaises(ValueError):
                AxisData(log_amplitude_C=value)
        with self.assertRaises(ValueError):
            complex_amplitude_bound(AxisData(),margin=0)
        with self.assertRaises(ValueError):
            axis_primitive(AxisData(),2)


if __name__=="__main__":
    unittest.main()
