"""Momentum residuals for the paper's explicit heat exterior and a Gaussian toy.

Viscosity is dimensionless and defaults to one. The exterior is only evaluated
for positive radius; it is not extended through the axis or joined to a core.
"""

from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.special import roots_genlaguerre, roots_jacobi, gamma, poch
from .scaling import positive
from .energy import CoreEnergyScaling, gaussian_swirl


@lru_cache(maxsize=32)
def _laguerre(order, h):
    x, w = roots_genlaguerre(order, h)
    return x, w/gamma(1+h)


@lru_cache(maxsize=32)
def _pressure_rule(order, h):
    x, w = roots_jacobi(order, 0, 2*h)
    return (x+1)/2, w/2**(1+2*h)


@dataclass(frozen=True)
class HeatExterior:
    """Normalized explicit family from Eq. (4.29), not a calibrated full profile."""

    h: float = 1e-6
    amplitude: float = 1.0
    order: int = 64

    def __post_init__(self):
        if np.ndim(self.h) != 0 or not np.isfinite(self.h) or not 0 <= self.h < 0.5:
            raise ValueError("h must be scalar in [0,0.5); full proof restrictions not checked")
        positive(self.amplitude, "amplitude")
        if not isinstance(self.order, (int, np.integer)) or not 4 <= self.order <= 256:
            raise ValueError("quadrature order must be an integer from 4 to 256")

    @property
    def A(self):
        return 0.5+self.h

    def profile(self, zeta, derivative=0):
        zeta = np.asarray(zeta, dtype=float)
        if np.any(~np.isfinite(zeta)) or np.any(zeta < 0):
            raise ValueError("profile argument must be finite and nonnegative")
        if derivative not in (0, 1, 2, 3):
            raise ValueError("implemented derivative orders are 0 through 3")
        v, weights = _laguerre(self.order, self.h)
        kernel = v**derivative * (1+zeta[..., None]*v)**(-self.h-derivative)
        factor = (-1)**derivative * (poch(self.h, derivative) if derivative else 1)
        return factor*np.sum(kernel*weights, axis=-1)

    def velocity(self, radius, tau):
        r, tau = np.broadcast_arrays(positive(radius, "radius"), positive(tau, "tau"))
        radial_square = r*r/2
        return self.amplitude*radial_square**(-self.A)*self.profile(2*tau/radial_square)

    def pressure(self, radius, tau, order=64):
        """Pressure tending to zero at radial infinity, computed by integration."""
        r, tau = np.broadcast_arrays(positive(radius, "radius"), positive(tau, "tau"))
        radial_square = r*r/2
        zeta = 2*tau/radial_square
        v, weights = _pressure_rule(order, self.h)
        integral = np.sum(weights*self.profile(zeta[..., None]*v)**2, axis=-1)
        return -self.amplitude**2/2 * radial_square**(-2*self.A)*integral

    def derivatives(self, radius, tau):
        r, tau = np.broadcast_arrays(positive(radius, "radius"), positive(tau, "tau"))
        s = r*r/2
        zeta = 2*tau/s
        H, H1, H2 = [self.profile(zeta, n) for n in range(3)]
        A, c = self.A, self.amplitude
        u = c*s**(-A)*H
        us = c*s**(-A-1)*(-A*H-zeta*H1)
        uss = c*s**(-A-2)*(A*(A+1)*H+2*(A+1)*zeta*H1+zeta*zeta*H2)
        ur = r*us
        urr = us+r*r*uss
        ut = -2*c*s**(-A-1)*H1
        # Algebraically combine the radial vector Laplacian before evaluation.
        # 2*h*(1+h) avoids subtracting two nearly equal O(1) coefficients.
        lap = c*s**(-A-1)*(2*self.h*(1+self.h)*H
                            +(4*A+2)*zeta*H1+2*zeta*zeta*H2)
        raw_lap = urr+ur/r-u/(r*r)
        return dict(velocity=u, time=ut, radial=ur, radial_second=urr,
                    laplacian=lap, raw_laplacian=raw_lap)

    def terms(self, radius, tau):
        r = positive(radius, "radius")
        d = self.derivatives(r, tau)
        zero = np.zeros_like(d["velocity"])
        def vector(radial, angular, axial):
            return np.stack(np.broadcast_arrays(radial, angular, axial), axis=-1)
        time = vector(zero, d["time"], zero)
        transport = vector(-d["velocity"]**2/r, zero, zero)
        pressure_gradient = vector(d["velocity"]**2/r, zero, zero)
        viscous = vector(zero, -d["laplacian"], zero)
        residual = time+transport+pressure_gradient+viscous
        scale = abs(time)+abs(transport)+abs(pressure_gradient)+abs(viscous)
        normalized = np.divide(abs(residual), scale, out=np.zeros_like(scale), where=scale>0)
        return dict(time=time, transport=transport, pressure=pressure_gradient,
                    minus_viscosity=viscous, residual=residual,
                    component_scale=scale, component_normalized=normalized)


def gaussian_pressure(radius, z, tau, model=None):
    model = model or CoreEnergyScaling()
    a, b, U = model.widths_and_speed(tau)
    return -0.5*U*U*np.exp(1-(np.asarray(radius)/a)**2-(np.asarray(z)/b)**2)


def gaussian_terms(radius, z, tau, model=None, viscosity=1.0):
    """Exact cylindrical residual for the prescribed D2 swirl.

    Pressure balances centrifugal acceleration; its axial gradient remains.
    The returned force is manufactured and is not smooth through tau=0.
    """
    model = model or CoreEnergyScaling()
    positive(viscosity, "viscosity")
    r, z, tau = np.broadcast_arrays(positive(radius, "radius"), np.asarray(z, float),
                                    positive(tau, "tau"))
    a, b, U = model.widths_and_speed(tau)
    R, Z = r/a, z/b
    u = gaussian_swirl(r, z, a, b, U)
    ut = u/tau*(1+model.h-0.5*R*R-(0.5-model.h)*Z*Z)
    lap = u*((R*R-4)/(a*a)+(Z*Z-1)/(b*b))
    dp_r = u*u/r
    dp_z = U*U*z/(b*b)*np.exp(1-R*R-Z*Z)
    zero = np.zeros_like(u)
    time = np.stack((zero, ut, zero), axis=-1)
    transport = np.stack((-dp_r, zero, zero), axis=-1)
    pressure = np.stack((dp_r, zero, dp_z), axis=-1)
    viscous = np.stack((zero, -viscosity*lap, zero), axis=-1)
    residual = time+transport+pressure+viscous
    scale = abs(time)+abs(transport)+abs(pressure)+abs(viscous)
    normalized = np.divide(abs(residual), scale, out=np.zeros_like(scale), where=scale>0)
    return dict(time=time, transport=transport, pressure=pressure,
                minus_viscosity=viscous, residual=residual,
                component_scale=scale, component_normalized=normalized)
