"""The finite moment solve from NS Eq. (5.16), with explicit representative bumps."""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np


@lru_cache(maxsize=1)
def bump_normalization():
    nodes,weights=np.polynomial.legendre.leggauss(128)
    return float(np.dot(weights,np.exp(-1/(1-nodes*nodes))))


@lru_cache(maxsize=16)
def quadrature(order):
    if not isinstance(order,(int,np.integer)) or order<4:
        raise ValueError("quadrature order must be an integer >=4")
    nodes,weights=np.polynomial.legendre.leggauss(order)
    raw=np.exp(-1/(1-nodes*nodes))
    return nodes,weights*raw/np.dot(weights,raw)


@dataclass(frozen=True)
class SmoothBump:
    center: float
    halfwidth: float

    def __post_init__(self):
        if not np.isfinite(self.center) or not np.isfinite(self.halfwidth) or not 0<self.halfwidth<self.center:
            raise ValueError("bump must have a positive width and support away from zero")

    def value(self,radius):
        radius=np.asarray(radius,float)
        z=(radius-self.center)/self.halfwidth
        value=np.zeros_like(radius)
        active=abs(z)<1
        # A separate high-order normalization is used for direct profile evaluation.
        value[active]=np.exp(-1/(1-z[active]**2))/(self.halfwidth*bump_normalization())
        return value

    def moment(self,power,order=96):
        nodes,weights=quadrature(order)
        return float(np.dot(weights,(self.center+self.halfwidth*nodes)**power))


def moment_matrices(lam,u_bumps,e_bumps,order=96):
    if not np.isfinite(lam) or lam<=0:
        raise ValueError("lambda must be positive; the U moment powers coincide at zero")
    if len(u_bumps)!=2 or len(e_bumps)!=3:
        raise ValueError("the source solve uses two axial and three angular bumps")
    ordered=sorted([*u_bumps,*e_bumps],key=lambda b:b.center)
    if any(a.center+a.halfwidth>=b.center-b.halfwidth for a,b in zip(ordered,ordered[1:])):
        raise ValueError("all bump supports must be disjoint")
    BU=np.array([[b.moment(power,order) for b in u_bumps] for power in (1,1-2*lam)])
    BE=np.array([[b.moment(power,order) for b in e_bumps] for power in (2,-2-2*lam,-2*lam)])
    return BU,BE


def solve_moments(lam,u_bumps,e_bumps,d_u,d_e,order=96):
    BU,BE=moment_matrices(lam,u_bumps,e_bumps,order)
    return np.linalg.solve(BU,-np.asarray(d_u,float)),np.linalg.solve(BE,-np.asarray(d_e,float))


def combine_bumps(radius,coefficients,bumps):
    return sum(coefficient*bump.value(radius) for coefficient,bump in zip(coefficients,bumps))


def original_moments(d_u,d_e,ef=1.):
    """Invert the paper's scaled discrepancy definitions; preserve moment order 1..5."""
    return np.array([d_u[0],d_e[0],2*ef*d_e[1],ef*d_u[1],-ef*d_e[2]],float)


def repaired_moments(BU,BE,alpha,beta,before,ef=1.):
    axial=BU@alpha
    angular=BE@beta
    return np.asarray(before)+np.array([axial[0],angular[0],2*ef*angular[1],ef*axial[1],-ef*angular[2]])
