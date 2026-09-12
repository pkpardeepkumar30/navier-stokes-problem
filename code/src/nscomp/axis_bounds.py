"""Logarithmic axis amplitude and a conservative complex-neighborhood bound.

The bound uses exact rational arithmetic on declared decimal inputs. It proves
the amplitude inequality B.16 alone, not the full construction's parameter
hierarchy or an inner-series error bound.
"""
from fractions import Fraction
from functools import lru_cache
import mpmath as mp


def _decimal_upper(value, places=40):
    scale = 10**places
    integer = (value.numerator*scale+value.denominator-1)//value.denominator
    whole, fractional = divmod(integer, scale)
    return f"{whole}.{fractional:0{places}d}"


def complex_amplitude_bound(data, margin=1, places=40):
    """A rectangle enclosing [-1,1] with explicit pole separation.

    Re(z) in (-1-b,1+b), |Im(z)|<b, b=min(1/1000,sigma/68).
    On this rectangle |H'|<17 and |H(z)-H(Re z)|<=sigma/4.
    Therefore |zeta|<=4*(1+2h*(1+2b)^2)/(3sigma).
    A straight integration path gives the primitive bound used for log(C).
    """
    h = Fraction(str(data.h))
    j = Fraction(str(data.offset_j))
    sigma = Fraction(str(data.regularization_sigma))
    lam = Fraction(str(data.radial_scale_Lambda))
    if not 0<h<Fraction(1,100) or not 0<j<=Fraction(1,20):
        raise ValueError("Require 0<h<.01 and 0<j<=.05")
    if sigma<=0 or lam<1 or margin<1:
        raise ValueError("Require sigma>0, Lambda>=1, and margin>=1")
    b = min(Fraction(1,1000), sigma/68)
    radius = 1+2*b
    derivative_bound = Fraction(9,2)+12*radius**2+2*j*radius
    if not derivative_bound<17:
        raise ArithmeticError("The stated derivative bound was not established")
    L_bound = 1+2*h*radius**2
    zeta_bound = 4*L_bound/(3*sigma)
    primitive_bound = radius*zeta_bound
    required_log_C = lam*primitive_bound+Fraction(str(margin))
    log_C = _decimal_upper(required_log_C, places)
    assert Fraction(log_C)>=required_log_C
    return dict(half_width=b,modulus_bound=radius,H_derivative_bound=derivative_bound,
                H_factor_modulus_lower=3*sigma/4,zeta_modulus_bound=zeta_bound,
                primitive_modulus_bound=primitive_bound,
                required_log_C=required_log_C,log_C_upper=log_C,
                margin=Fraction(str(margin)))


@lru_cache(maxsize=32)
def _poles(h_text, j_text, sigma_text, digits):
    with mp.workdps(digits):
        h,j,sigma=map(mp.mpf,(h_text,j_text,sigma_text))
        D=mp.mpf(".5")-h
        result=[]
        for sign in (-1,1):
            roots=mp.polyroots([-4,-j,D+4,j-sign*mp.j*sigma],
                               maxsteps=100,error=False,extraprec=35)
            for root in roots:
                derivative=D+4-12*root**2-2*j*root
                residue=-(1-2*h*root**2)/(2*derivative)
                result.append((root,residue))
        return tuple(result)


def axis_primitive(data, eta, digits=70):
    """Integral_0^eta zeta via rational partial fractions on the real interval.

    Each logarithm is continued from zero along a straight real path. No pole
    lies on that path. Conjugate poles make the final primitive real.
    """
    with mp.workdps(digits):
        e=mp.mpf(str(eta))
        if not mp.isfinite(e) or abs(e)>1:
            raise ValueError("Require real eta in [-1,1]")
        poles=_poles(str(data.h),str(data.offset_j),str(data.regularization_sigma),digits)
        value=mp.fsum(residue*mp.log(1-e/root) for root,residue in poles)
        if abs(mp.im(value))>mp.mpf(10)**(-digits+10):
            raise ArithmeticError("Conjugate pole sum acquired an unexpected imaginary part")
        return +mp.re(value)


def complex_axis_poles(data, digits=80):
    """Poles of zeta for numerical source auditing; bounds do not use root samples."""
    return _poles(str(data.h),str(data.offset_j),str(data.regularization_sigma),digits)


def logarithmic_axis_value(data, eta, digits=70):
    """Return log(g) and g with guard digits for the large logarithm."""
    log_C_text=data.log_amplitude_C
    if log_C_text is None:
        raise ValueError("An explicit log_amplitude_C is required")
    with mp.workdps(30):
        magnitude=max(abs(mp.mpf(log_C_text)),abs(mp.mpf(str(data.radial_scale_Lambda))),1)
        guard=max(0,int(mp.floor(mp.log10(magnitude))))+12
    with mp.workdps(digits+guard):
        primitive=axis_primitive(data,eta,digits+guard)
        log_g=mp.mpf(str(data.radial_scale_Lambda))*primitive-mp.mpf(log_C_text)
        return log_g,mp.exp(log_g)


def feature_scales(data, pressure, digits=80):
    """Locate the zero of H_star and the pressure-dependent narrower scale."""
    with mp.workdps(digits):
        h,j=mp.mpf(str(data.h)),mp.mpf(str(data.offset_j))
        A,D=mp.mpf(".5")+h,mp.mpf(".5")-h
        H=lambda e:D*e+(1-e*e)*(4*e+j)
        center=mp.findroot(H,(-j/4,mp.mpf(0)))
        derivative=D+4-12*center**2-2*j*center
        pi,pi_eta=pressure.eta_coefficients(center,1,digits)
        U=4*center+j
        Z=-A*(1-2*center*U)*U-(1-center**2)*pi_eta+4*A*center*pi
        return dict(center=center,H_derivative=derivative,Z_at_center=Z,
            chi_width=mp.mpf(str(data.regularization_sigma))/derivative,
            pressure_width=Z/(mp.mpf(str(data.radial_scale_Lambda))*derivative))


def exit_diagnostics(series,Y=4,order=None,axial_cap=3):
    """Pointwise diagnostics; a capped axial term gives an exit lower bound."""
    with mp.workdps(series.digits):
        values=series.evaluate(Y,order=order)
        e=mp.mpf(str(series.eta_center))
        h,j=mp.mpf(str(series.data.h)),mp.mpf(str(series.data.offset_j))
        lam=mp.mpf(str(series.data.radial_scale_Lambda))
        L=1-2*h*e*e
        H=(mp.mpf(".5")-h)*e+(1-e*e)*(4*e+j)
        chi=H*H/(H*H+mp.mpf(str(series.data.regularization_sigma))**2)
        target=mp.mpf("2.5")+mp.mpf(".95")*lam*L*chi
        lower=log_ratio=None
        if values["p1"]>0:
            log_ratio=2*values["log_abs_p2"]-mp.log(values["p1"])
            if chi>mp.mpf(".99"):
                lower=values["p1"]  # The angular term alone suffices if this exceeds 2.2.
            else:
                axial=mp.mpf(axial_cap) if log_ratio>=mp.log(axial_cap) else mp.exp(log_ratio)
                lower=values["p1"]+axial
        return {**values,"chi":chi,"scaled_p1":lam*values["p1"],
                "Sq_target":target,"Sq_margin":values["Sq"]-target,
                "exit_lower":lower,"log_axial_ratio":log_ratio}
