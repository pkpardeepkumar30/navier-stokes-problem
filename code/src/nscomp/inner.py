"""Regular radial coefficients for the source's leading inner profile equations.

The solver accepts explicit analytic axis data. The built-in pressure trace is a
validation datum, NOT the pressure trace of an assembled Lemma 4.8 outer profile.
Radial order is not a correction order in the paper's q^(2nh) expansion.
"""
from dataclasses import dataclass
from functools import lru_cache
import math
import mpmath as mp
from .axis_bounds import logarithmic_axis_value


class _Jet:
    """Truncated ordinary Taylor coefficients in eta around one evaluation point."""
    def __init__(self, values):
        self.c = list(values)

    @classmethod
    def constant(cls, value, order):
        return cls([mp.mpf(value)]+[mp.mpf(0)]*order)

    def _other(self, value):
        return value if isinstance(value,_Jet) else self.constant(value,len(self.c)-1)

    def __add__(self, value):
        value=self._other(value)
        return _Jet([a+b for a,b in zip(self.c,value.c)])

    __radd__=__add__

    def __neg__(self):
        return _Jet([-x for x in self.c])

    def __sub__(self,value):
        return self+-self._other(value)

    def __rsub__(self,value):
        return self._other(value)+-self

    def __mul__(self,value):
        value=self._other(value)
        n=len(self.c)
        return _Jet([mp.fsum(self.c[j]*value.c[i-j] for j in range(i+1)) for i in range(n)])

    __rmul__=__mul__

    def inverse(self):
        if not self.c[0]:
            raise ZeroDivisionError("Taylor jet has zero constant coefficient")
        out=[1/self.c[0]]
        for n in range(1,len(self.c)):
            out.append(-mp.fsum(self.c[j]*out[n-j] for j in range(1,n+1))/self.c[0])
        return _Jet(out)

    def __truediv__(self,value):
        if isinstance(value,_Jet):
            return self*value.inverse()
        return _Jet([x/value for x in self.c])

    def derivative(self):
        return _Jet([(i+1)*self.c[i+1] for i in range(len(self.c)-1)]+[mp.mpf(0)])

    def value(self,offset=0,derivative=0):
        if offset == 0:
            return mp.factorial(derivative)*self.c[derivative] if derivative<len(self.c) else mp.mpf(0)
        coeff=self.c
        for _ in range(derivative):
            coeff=[(i+1)*coeff[i+1] for i in range(len(coeff)-1)]
        return mp.polyval(list(reversed(coeff)),offset) if coeff else mp.mpf(0)


def _product_coefficient(a,b,n):
    return sum((a[i]*b[n-i] for i in range(n+1)),0)


@dataclass(frozen=True)
class AxisData:
    h: float = 1e-6
    offset_j: float = .02
    regularization_sigma: float = .2
    radial_scale_Lambda: float = 64.
    amplitude_C: float = 1e4
    pressure_magnitude: float = 100.
    log_amplitude_C: str | None = None

    def __post_init__(self):
        if not 0 < self.h < .01 or not 0 < self.offset_j <= .05:
            raise ValueError("Require 0<h<.01 and 0<j<=.05 for this validation family")
        for name in ("regularization_sigma","radial_scale_Lambda","amplitude_C","pressure_magnitude"):
            if not math.isfinite(getattr(self,name)) or getattr(self,name)<=0:
                raise ValueError(f"{name} must be positive")
        if self.radial_scale_Lambda<1 or self.amplitude_C<=1:
            raise ValueError("Require Lambda>=1 and C>1")
        if self.log_amplitude_C is not None:
            value=mp.mpf(str(self.log_amplitude_C))
            if not mp.isfinite(value) or value<=0:
                raise ValueError("log_amplitude_C must be finite and positive")


class InnerSeries:
    """Local radial polynomial in Y=Lambda*X, with eta Taylor jets.

    Jet degree exceeds radial order by four, so eta derivatives needed for
    coefficient generation and center evaluation remain available. Values at
    nonzero eta offset are only a local Taylor representation, not a global fit.
    An optional axis_pressure object supplies eta_coefficients(center,order,digits);
    otherwise the analytic validation pressure in AxisData is used.
    """
    def __init__(self,data=AxisData(),eta_center=0.,order=10,digits=60,axis_pressure=None):
        if not isinstance(order,int) or not 1<=order<=24:
            raise ValueError("order must be an integer between 1 and 24")
        if digits<30 or not math.isfinite(eta_center) or abs(eta_center)>1:
            raise ValueError("Require at least 30 digits and eta_center in [-1,1]")
        self.data,self.eta_center,self.order,self.digits=data,eta_center,order,digits
        with mp.workdps(digits):
            njet=order+4
            c=lambda value:_Jet.constant(value,njet)
            eta=c(str(eta_center))
            eta.c[1]=mp.mpf(1)
            h=mp.mpf(str(data.h))
            lam=mp.mpf(str(data.radial_scale_Lambda))
            sig=mp.mpf(str(data.regularization_sigma))
            amp=mp.mpf(str(data.amplitude_C))
            A,D=mp.mpf(".5")+h,mp.mpf(".5")-h
            d=1-eta*eta
            L=1-2*h*eta*eta
            ustar=4*eta+mp.mpf(str(data.offset_j))
            Hstar=D*eta+d*ustar
            zeta=-L*Hstar/(Hstar*Hstar+sig*sig)
            xi=lam*zeta
            def zeta_scalar(w):
                H=D*w+(1-w*w)*(4*w+mp.mpf(str(data.offset_j)))
                return -(1-2*h*w*w)*H/(H*H+sig*sig)
            if data.log_amplitude_C is None:
                g0=mp.exp(lam*mp.quad(zeta_scalar,[0,mp.mpf(str(eta_center))]))/amp
                self.log_g0=mp.log(g0)
            else:
                self.log_g0,g0=logarithmic_axis_value(data,eta_center,digits)
            g=[g0]
            # g_eta = xi*g, an analytic coefficient recurrence with no numerical differentiation.
            for n in range(njet):
                g.append(mp.fsum(xi.c[j]*g[n-j] for j in range(n+1))/(n+1))
            self.g=_Jet(g)
            self.xi=xi
            if axis_pressure is None:
                pi0=-mp.mpf(str(data.pressure_magnitude))*(1+eta*eta).inverse()*(1+eta*eta).inverse()
            else:
                # An external trace supplies ordinary eta coefficients at this center.
                # Its own quadrature uncertainty is not removed by mp arithmetic.
                coefficients=axis_pressure.eta_coefficients(eta_center,njet,digits)
                if len(coefficients)!=njet+1:
                    raise ValueError("Axis pressure returned the wrong coefficient count")
                pi0=_Jet(coefficients)
            phi,U,Pi=[c(1)],[ustar],[pi0]
            g2=self.g*self.g
            for n in range(order):
                avg=[U[i]/(i+1) for i in range(n+1)]
                W=[(1 if i==0 else 0)-2*D*eta*avg[i]-d*avg[i].derivative() for i in range(n+1)]
                H=[(D*eta if i==0 else 0)+d*U[i] for i in range(n+1)]
                angular=(_product_coefficient(W,[(i+1)*phi[i] for i in range(n+1)],n)
                    +h*phi[n]-2*h*eta*_product_coefficient(U,phi,n)
                    +_product_coefficient(H,[f.derivative()+xi*f for f in phi],n))
                axial=(_product_coefficient(W,[i*U[i] for i in range(n+1)],n)
                    +A*U[n]-2*A*eta*_product_coefficient(U,U,n)
                    +_product_coefficient(H,[u.derivative() for u in U],n)
                    +d*Pi[n].derivative()-4*A*eta*Pi[n]-2*eta*n*Pi[n])
                next_phi=angular/(2*lam*(n+1)*(n+2))/L
                next_u=axial/(2*lam*(n+1)**2)/L
                next_pi=g2*_product_coefficient(phi,phi,n)/(lam*(n+1))
                phi.append(next_phi)
                U.append(next_u)
                Pi.append(next_pi)
            self.phi,self.U,self.Pi=phi,U,Pi

    def _evaluate(self,values,Y,eta_offset=0,dy=0,de=0,order=None,average=False):
        chosen=self.order if order is None else order
        coeff=[v.value(eta_offset,de)/(i+1 if average else 1) for i,v in enumerate(values[:chosen+1])]
        for _ in range(dy):
            coeff=[(i+1)*coeff[i+1] for i in range(len(coeff)-1)]
        return mp.polyval(list(reversed(coeff)),Y) if coeff else mp.mpf(0)

    def evaluate(self,Y,eta_offset=0,order=None):
        """Return high-precision profile values and direct leading-balance diagnostics."""
        if order is not None and not 0<=order<=self.order:
            raise ValueError("evaluation order must be between zero and computed order")
        with mp.workdps(self.digits):
            Y,delta=mp.mpf(str(Y)),mp.mpf(str(eta_offset))
            eta=mp.mpf(str(self.eta_center))+delta
            if not mp.isfinite(Y) or Y<0 or not mp.isfinite(eta) or abs(eta)>1:
                raise ValueError("Require finite Y>=0 and eta in [-1,1]")
            h=mp.mpf(str(self.data.h))
            lam=mp.mpf(str(self.data.radial_scale_Lambda))
            A,D=mp.mpf(".5")+h,mp.mpf(".5")-h
            d,L=1-eta*eta,1-2*h*eta*eta
            val=lambda values,dy=0,de=0,avg=False:self._evaluate(values,Y,delta,dy,de,order,avg)
            phi,py,pyy,pe=val(self.phi),val(self.phi,1),val(self.phi,2),val(self.phi,de=1)
            U,uy,uyy,ue=val(self.U),val(self.U,1),val(self.U,2),val(self.U,de=1)
            Pi,piy,pie=val(self.Pi),val(self.Pi,1),val(self.Pi,de=1)
            avg,avg_eta=val(self.U,avg=True),val(self.U,de=1,avg=True)
            W=1-2*D*eta*avg-d*avg_eta
            H=D*eta+d*U
            g=self.g.value(delta)
            geta=self.g.value(delta,1)
            angular_terms=[2*L*lam*(Y*pyy+2*py),
                -W*(phi+Y*py),-h*(1-2*eta*U)*phi,-H*(pe+geta/g*phi)]
            axial_terms=[2*L*lam*(Y*uyy+uy),-W*Y*uy,-A*(1-2*eta*U)*U,
                -H*ue,-d*pie,4*A*eta*Pi,2*eta*Y*piy]
            pressure_terms=[lam*piy,-g*g*phi*phi]
            residual=lambda terms:mp.fsum(terms)
            relative=lambda terms:abs(residual(terms))/mp.fsum(abs(x) for x in terms) if any(terms) else mp.mpf(0)
            chosen=self.order if order is None else order
            pi_increment=mp.fsum(self.Pi[i].value(delta)*Y**i for i in range(1,chosen+1))
            v0_over_X=(2*eta*U+W-1)/L
            p1=-2*Y*py/phi
            ns=-2*lam*uy
            log_abs_p2=(mp.log(Y/(2*lam))/2+mp.log(abs(ns))-mp.log(abs(g*phi))
                        if Y>0 and ns else mp.ninf)
            Sq=-2*L*lam*(Y*pyy+2*py)/phi
            return dict(Phi=phi,U=U,Pi=Pi,pressure_increment=pi_increment,F=g*phi,V0_over_X=v0_over_X,
                p1=p1,n_s=ns,log_abs_p2=log_abs_p2,Sq=Sq,
                angular_balance=residual(angular_terms),axial_balance=residual(axial_terms),
                pressure_balance=residual(pressure_terms),
                angular_relative=relative(angular_terms),axial_relative=relative(axial_terms),
                pressure_relative=relative(pressure_terms))

    def physical_field(self,r,z,t):
        """Local polynomial field for independent physical-coordinate finite differences.

        Viscosity is one. eta must stay close to the series center for this local
        eta Taylor representation. The supplied r must be strictly positive.
        """
        with mp.workdps(self.digits):
            r,z,t=mp.mpf(str(r)),mp.mpf(str(z)),mp.mpf(str(t))
            if not all(mp.isfinite(x) for x in (r,z,t)) or r<=0 or t>=1:
                raise ValueError("Require r>0 and t<1")
            h=mp.mpf(str(self.data.h))
            D,A=mp.mpf(".5")-h,mp.mpf(".5")+h
            tau=1-t
            # eta=z/q^D and q-z^2 q^(2h)=tau.
            q=mp.findroot(lambda q:q-z*z*q**(2*h)-tau,(tau+z*z,tau+z*z+mp.mpf(".01")))
            eta=z/q**D
            X=r*r/(2*q)
            s=self.evaluate(self.data.radial_scale_Lambda*X,eta_offset=eta-mp.mpf(str(self.eta_center)))
            ur=r*s["V0_over_X"]/(2*q)
            swirl=q**(-A)*mp.sqrt(2*X)*s["F"]
            uz=q**(-A)*s["U"]
            pressure=q**(-2*A)*s["Pi"]
            return ur,swirl,uz,pressure


@lru_cache(maxsize=64)
def inner_series(data=AxisData(),eta_center=0.,order=10,digits=60,axis_pressure=None):
    return InnerSeries(data,eta_center,order,digits,axis_pressure)
