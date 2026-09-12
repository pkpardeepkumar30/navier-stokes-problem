"""Independent checks for the source-defined leading inner recurrence."""
import unittest
import mpmath as mp
from nscomp.inner import AxisData,inner_series


class InnerChecks(unittest.TestCase):
    def test_axis_values_and_first_coefficients(self):
        data=AxisData(radial_scale_Lambda=512)
        s=inner_series(data,.05,10,65)
        with mp.workdps(65):
            eta=mp.mpf(".05")
            h=mp.mpf(str(data.h)); A=mp.mpf(".5")+h; D=mp.mpf(".5")-h
            d,L=1-eta*eta,1-2*h*eta*eta
            U=4*eta+mp.mpf(str(data.offset_j))
            H=D*eta+d*U
            W=1-4*d-2*D*eta*U
            pi=-100/(1+eta*eta)**2
            pi_eta=400*eta/(1+eta*eta)**3
            xi=-data.radial_scale_Lambda*L*H/(H*H+mp.mpf(".2")**2)
            expected_phi=(W+h*(1-2*eta*U)+H*xi)/(4*data.radial_scale_Lambda*L)
            expected_u=(A*(1-2*eta*U)*U+4*H+d*pi_eta-4*A*eta*pi)/(2*data.radial_scale_Lambda*L)
            self.assertLess(abs(s.phi[1].c[0]-expected_phi),mp.mpf("1e-60"))
            self.assertLess(abs(s.U[1].c[0]-expected_u),mp.mpf("1e-60"))
            axis=s.evaluate(0)
            self.assertEqual(axis["Phi"],1)
            self.assertLess(abs(axis["U"]-U),mp.mpf("1e-60"))
            self.assertGreater(axis["F"],0)
            self.assertLess(abs(s.Pi[1].c[0]-s.g.c[0]**2/data.radial_scale_Lambda),mp.mpf("1e-60"))

    def test_degree_refinement_of_direct_balances(self):
        s=inner_series(AxisData(radial_scale_Lambda=512),.2,18,70)
        errors=[]
        for n in (4,8,12,18):
            d=s.evaluate(2,order=n)
            errors.append(max(float(d[k]) for k in ("angular_relative","axial_relative","pressure_relative")))
        self.assertTrue(all(b<a/100 for a,b in zip(errors[:-1],errors[1:])),errors)
        self.assertLess(errors[-1],1e-17)

    def test_precision_refinement(self):
        data=AxisData(radial_scale_Lambda=512)
        low=inner_series(data,.2,12,45).evaluate(1)
        high=inner_series(data,.2,12,80).evaluate(1)
        with mp.workdps(80):
            for key in ("Phi","U","Pi","F"):
                self.assertLess(abs((low[key]-high[key])/high[key]),mp.mpf("1e-40"))

    def test_independently_centered_eta_expansions(self):
        data=AxisData(radial_scale_Lambda=512)
        first=inner_series(data,.05,12,70).evaluate(.5,eta_offset=.0001)
        second=inner_series(data,.0501,12,70).evaluate(.5)
        with mp.workdps(70):
            for key in ("Phi","U","Pi","F"):
                self.assertLess(abs((first[key]-second[key])/second[key]),mp.mpf("1e-18"))

    def test_large_lambda_comparison_has_expected_orders(self):
        with mp.workdps(60):
            eta,Y=mp.mpf(".2"),mp.mpf(1)
            h=mp.mpf("1e-6"); A=mp.mpf(".5")+h; D=mp.mpf(".5")-h
            d,L=1-eta*eta,1-2*h*eta*eta
            Ustar=4*eta+mp.mpf(".02"); H=D*eta+d*Ustar
            Z=-A*(1-2*eta*Ustar)*Ustar-4*H-d*400*eta/(1+eta*eta)**3+4*A*eta*(-100/(1+eta*eta)**2)
            chi=H*H/(H*H+mp.mpf(".2")**2)
            phi_limit=mp.hyper([],[2],-Y*chi/2)
            ep,eu=[],[]
            for lam in (256,512,1024):
                v=inner_series(AxisData(radial_scale_Lambda=lam),.2,12,60).evaluate(Y)
                ep.append(abs(v["Phi"]-phi_limit))
                eu.append(abs(v["U"]-(Ustar-Y*Z/(2*lam*L))))
            self.assertTrue(all(mp.mpf("1.8")<a/b<mp.mpf("2.2") for a,b in zip(ep[:-1],ep[1:])))
            self.assertTrue(all(mp.mpf("3.5")<a/b<mp.mpf("4.5") for a,b in zip(eu[:-1],eu[1:])))

    def test_physical_coordinate_divergence_and_leading_momentum(self):
        with mp.workdps(70):
            s=inner_series(AxisData(),.05,12,70)
            Y,q,eta=mp.mpf(".4"),mp.mpf(".8"),mp.mpf(".05")
            D=mp.mpf(".5")-mp.mpf(str(s.data.h))
            r=mp.sqrt(2*q*Y/s.data.radial_scale_Lambda)
            z=eta*q**D; t=1-q*(1-eta*eta)
            values=s.physical_field(r,z,t)
            step=mp.mpf("2e-5")
            def derivative(axis,part):
                samples=[]
                for shift in (-2,-1,1,2):
                    coords=[r,z,t]; coords[axis]+=shift*step
                    samples.append(s.physical_field(*coords)[part])
                return (samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*step)
            def radial_second(part):
                samples=[s.physical_field(r+shift*step,z,t)[part] for shift in (-2,-1,1,2)]
                return (-samples[3]+16*samples[2]-30*values[part]+16*samples[1]-samples[0])/(12*step**2)
            ur,w,u,p=values
            angular=[derivative(2,1),ur*derivative(0,1),u*derivative(1,1),ur*w/r,
                     -radial_second(1),-derivative(0,1)/r,w/r**2]
            axial=[derivative(2,2),ur*derivative(0,2),u*derivative(1,2),
                   derivative(1,3),-radial_second(2),-derivative(0,2)/r]
            self.assertLess(abs(mp.fsum(angular))/mp.fsum(abs(x) for x in angular),mp.mpf("2e-11"))
            self.assertLess(abs(mp.fsum(axial))/mp.fsum(abs(x) for x in axial),mp.mpf("1e-14"))
            self.assertLess(abs(derivative(0,0)+ur/r+derivative(1,2)),mp.mpf("1e-14"))
            self.assertLess(abs(derivative(0,3)-w*w/r),mp.mpf("1e-26"))

    def test_regular_axis_and_sampled_positivity(self):
        data=AxisData(radial_scale_Lambda=512)
        for eta in (-.6,0,.6):
            s=inner_series(data,eta,18,70)
            for Y in (0,.1,1,2,4.1):
                v=s.evaluate(Y)
                self.assertGreater(v["Phi"],0)
                self.assertTrue(mp.isfinite(v["V0_over_X"]))

    def test_invalid_inputs(self):
        for bad in (0,-1,float("nan")):
            with self.assertRaises(ValueError):
                AxisData(regularization_sigma=bad)
        with self.assertRaises(ValueError):
            inner_series(order=0)
        with self.assertRaises(ValueError):
            inner_series().evaluate(-1)


if __name__=="__main__":
    unittest.main()
