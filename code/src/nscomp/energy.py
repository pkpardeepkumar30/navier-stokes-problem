"""Concentration and energy models; none is a reconstructed paper solution."""

from dataclasses import dataclass
import math
import numpy as np
from .scaling import positive


@dataclass(frozen=True)
class CoreEnergyScaling:
    """Anchored paper powers, with h=0 allowed only as a limiting example.

    a0 and b0 are radial and axial Gaussian widths in chosen consistent units.
    b0 is NOT the full axial length. No physical units are calibrated by default.
    Larger h values are accepted for mathematical counterexamples, not as
    admissible choices in the paper.
    """

    h: float = 0.0
    a0: float = 1.0
    b0: float = 2.0
    u0: float = 1.0
    density: float = 1.0

    def __post_init__(self):
        for name in ("a0", "b0", "u0", "density"):
            value = getattr(self, name)
            if np.ndim(value) != 0:
                raise ValueError(f"{name} must be scalar")
            positive(value, name)
        if np.ndim(self.h) != 0 or not np.isfinite(self.h) or not 0 <= self.h < 0.5:
            raise ValueError("h must be scalar in [0, 0.5); this does not check admissibility")

    @property
    def energy_exponent(self):
        return 0.5 - 3*self.h

    def ratios(self, time_ratio):
        s = positive(time_ratio, "time ratio")
        return {
            "radial": s**0.5,
            "axial": s**(0.5-self.h),
            "speed": s**(-0.5-self.h),
            "volume": s**(1.5-self.h),
            "energy": s**self.energy_exponent,
            "slenderness": s**self.h,  # (a/b)/(a0/b0)
        }

    def widths_and_speed(self, time_ratio):
        r = self.ratios(time_ratio)
        return self.a0*r["radial"], self.b0*r["axial"], self.u0*r["speed"]

    def toy_energy(self, time_ratio):
        """Exact integral over R^3 of density*|u|^2/2 for gaussian_swirl."""
        a, b, u = self.widths_and_speed(time_ratio)
        return self.density/2 * math.e*math.pi**1.5 * u**2*a**2*b


def gaussian_swirl(r, z, a, b, amplitude):
    """u_theta = amplitude*R*exp((1-R^2-Z^2)/2), maximum at R=1,Z=0.

    ur=uz=0. The Cartesian vector is smooth at r=0 and divergence free.
    The time-varying profile is prescribed; it does not advect its own envelope.
    """
    positive(a, "radial width")
    positive(b, "axial width")
    positive(amplitude, "amplitude")
    radial, axial = np.broadcast_arrays(np.asarray(r, float), np.asarray(z, float))
    if np.any(~np.isfinite(radial)) or np.any(radial < 0) or np.any(~np.isfinite(axial)):
        raise ValueError("r must be finite and nonnegative; z must be finite")
    R, Z = radial/a, axial/b
    return amplitude*R*np.exp((1-R**2-Z**2)/2)


def swirl_cartesian(x, y, z, a, b, amplitude):
    """Axis-regular Cartesian evaluation, including exactly x=y=0."""
    positive(a, "radial width")
    positive(b, "axial width")
    positive(amplitude, "amplitude")
    x, y, z = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float),
                                  np.asarray(z, float))
    if any(np.any(~np.isfinite(q)) for q in (x, y, z)):
        raise ValueError("coordinates must be finite")
    factor = amplitude/a*np.exp((1-(x*x+y*y)/a**2-z*z/b**2)/2)
    return np.stack((-y*factor, x*factor, np.zeros_like(factor)), axis=-1)


def cylindrical_energy_quadrature(model, time_ratio, order=64, cutoff=8.0):
    """Gauss-Legendre quadrature with the cylindrical Jacobian 2*pi*r.

    Domain: 0<=r<=cutoff*a and |z|<=cutoff*b. Compare finite-domain quadrature
    with the whole-space analytic result only after estimating the Gaussian tail.
    """
    if not isinstance(order, (int, np.integer)) or order < 2:
        raise ValueError("quadrature order must be an integer >=2")
    positive(cutoff, "cutoff")
    if np.ndim(time_ratio) != 0:
        raise ValueError("quadrature takes a scalar time ratio")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    a, b, u = model.widths_and_speed(time_ratio)
    r = a*cutoff*(nodes+1)/2
    z = b*cutoff*nodes
    wr, wz = weights*a*cutoff/2, weights*b*cutoff
    velocity = gaussian_swirl(r[:, None], z[None, :], a, b, u)
    integrand = model.density/2*velocity**2 * (2*np.pi*r[:, None])
    return float(np.sum(integrand*wr[:, None]*wz[None, :]))


def gaussian_energy_fraction(cutoff):
    """Exact fraction inside the quadrature cylinder, not a core fraction."""
    positive(cutoff, "cutoff")
    return (1-(cutoff**2+1)*math.exp(-cutoff**2))*math.erf(cutoff)


def gaussian_energy_tail_fraction(cutoff):
    """Tail fraction evaluated without subtracting nearly equal numbers."""
    positive(cutoff, "cutoff")
    return math.erfc(cutoff)+(cutoff**2+1)*math.exp(-cutoff**2)*math.erf(cutoff)


def narrow_gaussian(x, epsilon, exponent=0.5):
    """Scalar 1D toy, f_epsilon=epsilon^-exponent*exp(-x^2/(2*epsilon^2))."""
    positive(epsilon, "epsilon")
    if not np.isfinite(exponent):
        raise ValueError("exponent must be finite")
    return epsilon**(-exponent)*np.exp(-0.5*(np.asarray(x)/epsilon)**2)


def narrow_gaussian_squared_integral(epsilon, exponent=0.5):
    return math.sqrt(math.pi)*positive(epsilon, "epsilon")**(1-2*exponent)


def material_map(points, time_ratio):
    """Exact trajectories of v=(-x/(2*tau),-y/(2*tau),z/tau), tau0=1.

    Separate affine incompressible toy on 0<=t<1. It has infinite whole-space
    energy and is used ONLY to explain volume-preserving material deformation.
    """
    s = positive(time_ratio, "time ratio")
    points = np.asarray(points, dtype=float)
    if np.ndim(s) != 0 or points.ndim < 1 or points.shape[-1] != 3:
        raise ValueError("provide a scalar time ratio and points of shape (...,3)")
    if np.any(~np.isfinite(points)):
        raise ValueError("points must be finite")
    return points*np.array([np.sqrt(s), np.sqrt(s), 1/s])
