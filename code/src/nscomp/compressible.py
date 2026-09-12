"""Stationary Gaussian isentropic Euler vortex and independent radial integration.

Lengths use the peak-speed radius a; velocities use the exterior sound speed.
This is a known equilibrium family, not the OpenAI collapse or a viscous solution.
"""
from dataclasses import dataclass
import numpy as np
from .scaling import positive


def nonnegative(value, name):
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)) or np.any(out < 0):
        raise ValueError(f"{name} must be finite and nonnegative")
    return out


@dataclass(frozen=True)
class IsentropicVortex:
    peak_mach: float = .3
    gamma: float = 1.4

    def __post_init__(self):
        nonnegative(self.peak_mach, "peak Mach")
        if not np.isfinite(self.gamma) or self.gamma <= 1:
            raise ValueError("gamma must exceed one")
        if self.depression >= 1:
            raise ValueError("This vortex would have nonpositive center temperature")

    @property
    def depression(self):
        return .5*(self.gamma-1)*np.e*self.peak_mach**2

    def state(self, radius):
        r = nonnegative(radius, "dimensionless radius")
        drop = self.depression*np.exp(-r*r)
        logtheta = np.log1p(-drop)
        theta = 1-drop
        density = np.exp(logtheta/(self.gamma-1))
        pressure = np.exp(self.gamma*logtheta/(self.gamma-1))
        swirl = self.peak_mach*r*np.exp((1-r*r)/2)
        return dict(temperature_ratio=theta,density_ratio=density,pressure_ratio=pressure,
            density_deficit=-np.expm1(logtheta/(self.gamma-1)),
            pressure_deficit=-np.expm1(self.gamma*logtheta/(self.gamma-1)),
            incompressible_pressure_ratio=1-self.gamma/(self.gamma-1)*drop,
            swirl_over_exterior_sound_speed=swirl,
            local_mach=swirl/np.sqrt(theta))

    def conservative_flux(self, x, y):
        """Cartesian fluxes in units rho_inf, c_inf, a; includes smooth axis value."""
        x,y = np.broadcast_arrays(np.asarray(x,float),np.asarray(y,float))
        radius = np.hypot(x,y)
        s = self.state(radius)
        omega = self.peak_mach*np.exp((1-radius**2)/2)
        u,v = -y*omega,x*omega
        rho = s["density_ratio"]
        pressure = s["pressure_ratio"]/self.gamma
        energy = pressure/(self.gamma-1)+.5*rho*(u*u+v*v)
        fx = np.stack((rho*u,rho*u*u+pressure,rho*u*v,(energy+pressure)*u),axis=-1)
        fy = np.stack((rho*v,rho*u*v,rho*v*v+pressure,(energy+pressure)*v),axis=-1)
        return fx,fy


def density_deficit_mach(deficit, gamma=1.4):
    """Peak Mach (using exterior sound speed) giving a declared center density deficit."""
    d = positive(deficit,"density deficit")
    if np.any(d >= 1) or not np.isfinite(gamma) or gamma <= 1:
        raise ValueError("Require 0 < deficit < 1 and gamma > 1")
    return np.sqrt(-2*np.expm1((gamma-1)*np.log1p(-d))/((gamma-1)*np.e))


def integrate_density(vortex, outer_radius, intervals, boundary="exact"):
    """Classical RK4 applied inward to the nonlinear radial density equation.

    rho' = M^2 R exp(1-R^2) rho^(2-gamma).
    The analytic profile is used ONLY to set the optional exact outer datum.
    No time step, two-dimensional CFD grid, or propagating-vortex test is implied.
    """
    outer = float(positive(outer_radius,"outer radius"))
    if not isinstance(intervals,(int,np.integer)) or intervals < 2:
        raise ValueError("intervals must be an integer at least 2")
    if boundary not in ("exact","far_field"):
        raise ValueError("boundary must be exact or far_field")
    r = np.linspace(outer,0,intervals+1)
    rho = np.empty_like(r)
    rho[0] = float(vortex.state(outer)["density_ratio"]) if boundary=="exact" else 1.
    step = -outer/intervals
    def rhs(x,value):
        if value <= 0 or not np.isfinite(value):
            raise ValueError("Nonpositive intermediate density: refine the radial grid")
        return vortex.peak_mach**2*x*np.exp(1-x*x)*value**(2-vortex.gamma)
    for i,x in enumerate(r[:-1]):
        y = rho[i]
        k1 = rhs(x,y)
        k2 = rhs(x+step/2,y+step*k1/2)
        k3 = rhs(x+step/2,y+step*k2/2)
        k4 = rhs(x+step,y+step*k3)
        rho[i+1] = y+step*(k1+2*k2+2*k3+k4)/6
    if np.any(rho <= 0):
        raise ValueError("Computed density is not positive")
    return r[::-1],rho[::-1]


def local_gas_audit(vortex, radius, radius_m, gas, sutherland_temperature_K=110.4):
    """Thermodynamic gradient Kn proxies and actual shear for this declared equilibrium.

    Mean free path follows the NIST Eq.11 T,p correction, used as a local-property
    approximation. This does not validate the inviscid or calorically perfect assumptions.
    """
    a = float(positive(radius_m,"peak radius"))
    if not np.isclose(vortex.gamma,gas.heat_capacity_ratio,rtol=0,atol=1e-14):
        raise ValueError("Vortex and gas must have the same gamma")
    st = float(positive(sutherland_temperature_K,"Sutherland temperature"))
    r = nonnegative(radius,"dimensionless radius")
    s = vortex.state(r)
    theta = s["temperature_ratio"]
    pressure = s["pressure_ratio"]
    mfp = (gas.mean_free_path_m * theta/pressure
           * (1+st/gas.temperature_K)/(1+st/(gas.temperature_K*theta)))
    dlogt = 2*r*vortex.depression*np.exp(-r*r)/(a*theta)
    dlogrho = dlogt/(vortex.gamma-1)
    dlogp = vortex.gamma*dlogrho
    # |du_theta/dr - u_theta/r|: twice the magnitude of the r-theta strain entry.
    shear_rate = gas.sound_speed/a*vortex.peak_mach*r*r*np.exp((1-r*r)/2)
    thermal_speed = gas.mean_thermal_speed*np.sqrt(theta)
    return dict(mean_free_path_m=mfp,kn_temperature=mfp*dlogt,
        kn_density=mfp*dlogrho,kn_pressure=mfp*dlogp,
        shear_rate_s_inv=shear_rate,collision_shear_ratio=mfp/thermal_speed*shear_rate,
        temperature_K=gas.temperature_K*theta,pressure_Pa=gas.pressure_Pa*pressure)
