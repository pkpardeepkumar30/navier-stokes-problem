"""Axis pressure from the ideal azimuthal schedule, source Lemma A.5.

The finite schedule implements Appendix A.2 with pressure-preserving angular
bumps omitted, as permitted in A.21. It does not certify the outer stress cone
or the complete existence hierarchy. Quadrature is binary64; logarithmic
weights retain very small stage contributions without underflow.
"""
from dataclasses import dataclass
from functools import lru_cache
import math
import mpmath as mp
import numpy as np
from scipy.special import expit, logsumexp


@lru_cache(maxsize=16)
def _rule(order):
    x, w = np.polynomial.legendre.leggauss(order)
    return (x+1)/2, w/2


def smooth_step(x):
    """The exact functional form in A.5, evaluated stably."""
    x = np.asarray(x, dtype=float)
    u = np.clip(x, 1e-12, 1-1e-12)
    value = expit(1/(1-u)**2-1/u**2)
    return np.where(x <= 0, 0., np.where(x >= 1, 1., value))


def step_derivative(x):
    x = np.asarray(x, dtype=float)
    u = np.clip(x, 1e-12, 1-1e-12)
    value = smooth_step(u)
    return np.where((x > 0) & (x < 1),
                    value*(1-value)*(2/u**3+2/(1-u)**3), 0.)


def integrated_step(x, order=128):
    x = np.asarray(x, dtype=float)
    nodes, weights = _rule(order)
    value = x*np.sum(smooth_step(x[..., None]*nodes)*weights, axis=-1)
    # Symmetry gives this endpoint exactly, independent of quadrature order.
    return np.where(x == 1, .5, value)


@dataclass(frozen=True)
class Schedule:
    M_d: float = 2.
    intermediate_lambda: float = .02
    h: float = 1e-9
    interpolation_length: float = 100.
    terminal_c: float = .01
    log_amplitude_margin: float = 1.

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive and finite")
        if self.M_d > 20 or self.intermediate_lambda >= .1:
            raise ValueError("This numerical implementation requires M_d<=20 and lambda<.1")
        if self.h >= min(.01, self.intermediate_lambda) or math.log(self.h) >= -self.T_d:
            raise ValueError("Require h<min(.01,lambda,exp(-T_d))")
        if self.terminal_c >= .05:
            raise ValueError("Require terminal_c<.05 in this numerical family")

    @property
    def T_d(self):
        return math.exp(self.M_d)+10

    @property
    def log_P(self):
        return self.T_d+self.log_amplitude_margin


@dataclass(frozen=True)
class Segment:
    name: str
    start_y: float
    length: float
    log_c_start: float
    ell_start: float
    ell_end: float
    theta: float
    kind: str = "power"
    rho: float = 0.

    def shape(self, y, order=128):
        """Return log(c/P_star) and theta in E/P_star=c*f^theta."""
        y = np.asarray(y, dtype=float)
        if self.kind == "slope":
            relative = ((self.ell_start-.5)*y
                        +(self.ell_end-self.ell_start)*self.length
                        *integrated_step(y/self.length, order))
            theta = np.full_like(y, self.theta)
        elif self.kind == "interpolation":
            step = smooth_step(y/self.length)
            relative = (self.ell_start-.5)*y-step*np.log(2.)
            theta = 1-step
        elif self.kind == "terminal":
            psi = 1-smooth_step((y-1)/2)
            relative = ((self.ell_start-.5)*y
                        +np.log1p(-self.rho*psi)-np.log1p(-self.rho))
            theta = np.zeros_like(y)
        else:
            relative = (self.ell_start-.5)*y
            theta = np.full_like(y, self.theta)
        return self.log_c_start+relative, theta


class PressureTrace:
    """Positive quadrature representation of -2*Pi0/P_star^2.

    Arbitrary-precision eta coefficients preserve algebraic digits and small
    contributions; their underlying quadrature weights remain binary64.
    """
    def __init__(self, schedule=Schedule(), order=128):
        if not isinstance(order, int) or not 8 <= order <= 512:
            raise ValueError("Require integer quadrature order between 8 and 512")
        self.schedule, self.order = schedule, order
        self.segments = []
        # The semi-infinite initial branch integrates to 5*f^2 exactly.
        self.components = [dict(stage="initial power", log_weight=math.log(5.), theta=1.)]
        y, logc = 0., 0.

        def add(name, length, ell0, ell1=None, theta=1., kind="power", rho=0.):
            nonlocal y, logc
            ell1 = ell0 if ell1 is None else ell1
            segment = Segment(name, y, length, logc, ell0, ell1, theta, kind, rho)
            self.segments.append(segment)
            self._integrate_segment(segment)
            logc = float(segment.shape(length, order)[0])
            y += length

        lam, h = schedule.intermediate_lambda, schedule.h
        add("initial slope transition", 1., .6, 0., kind="slope")
        add("axial decrease interval", schedule.T_d, 0.)
        add("intermediate slope transition", 1., 0., -lam, kind="slope")
        add("reserved intermediate interval", 60*math.log(1/lam), -lam)
        add("axial pulse interval", 13/lam, -lam)
        add("remove eta dependence", schedule.interpolation_length, -lam,
            kind="interpolation")
        add("uniform angular correction interval", 30*math.log(1/lam), -lam, theta=0.)

        self.release_start_y = y
        Q = (lam-h)/(1-lam)
        Q = self._advance_Q(Q, -lam, -1.)
        add("release to steep slope", 1., -lam, -1., theta=0., kind="slope")
        hold = 4*math.log(1/h)
        Q += (1-h)*hold  # ell=-1 makes the homogeneous coefficient vanish.
        add("steep power interval", hold, -1., theta=0.)
        Q = self._advance_Q(Q, -1., -h)
        add("release to heat exponent", 1., -1., -h, theta=0., kind="slope")

        nodes, weights = _rule(order)
        rho = schedule.terminal_c*h
        # In A.13 the f factors cancel in the integrating factor.
        v = 1+2*nodes
        Qp = rho/(1-rho)*np.sum(weights*np.exp((1-h)*v)*step_derivative(nodes))
        waiting_length = math.log(Q/Qp)/(1-h)
        if not waiting_length > 0:
            raise ValueError("The source's preterminal waiting interval must be positive")
        self.Q_before_wait, self.Q_target = float(Q), float(Qp)
        self.waiting_length = waiting_length
        add("wait for terminal Q target", waiting_length, -h, theta=0.)
        add("terminal smooth factor", 3., -h, theta=0., kind="terminal", rho=rho)
        self.tail_start_y, self.tail_log_c = y, logc
        self.components.append(dict(stage="infinite exterior power",
                                    log_weight=2*logc-math.log(1+2*h), theta=0.))
        self.log_weights = np.array([row["log_weight"] for row in self.components])
        self.thetas = np.array([row["theta"] for row in self.components])

    def _advance_Q(self, initial, ell0, ell1):
        nodes, weights = _rule(self.order)
        ell = ell0+(ell1-ell0)*smooth_step(nodes)
        G = (1+ell0)*nodes+(ell1-ell0)*integrated_step(nodes, self.order)
        G_end = 1+(ell0+ell1)/2
        return math.exp(-G_end)*(initial+np.sum(weights*np.exp(G)*(-ell-self.schedule.h)))

    def _integrate_segment(self, segment):
        if segment.kind == "power":
            decay = 1-2*segment.ell_start
            local_integral = -math.expm1(-decay*segment.length)/decay
            self.components.append(dict(stage=segment.name,
                log_weight=2*segment.log_c_start+math.log(local_integral),
                theta=segment.theta))
            return
        nodes, weights = _rule(self.order)
        # Split the terminal factor at its exactly constant first unit.
        intervals = ((0., 1.), (1., 3.)) if segment.kind == "terminal" else ((0., segment.length),)
        for left, right in intervals:
            y = left+(right-left)*nodes
            logc, theta = segment.shape(y, self.order)
            logw = np.log((right-left)*weights)+2*logc
            if segment.kind == "interpolation":
                self.components.extend(dict(stage=segment.name, log_weight=float(w), theta=float(t))
                                       for w, t in zip(logw, theta))
            else:
                self.components.append(dict(stage=segment.name,
                    log_weight=float(logsumexp(logw)), theta=float(theta[0])))

    def normalized(self, eta, derivative=0):
        """Pi0/P_star^2 and its first two eta derivatives."""
        eta = np.asarray(eta, dtype=float)
        if np.any(~np.isfinite(eta)) or np.any(np.abs(eta)>1):
            raise ValueError("Require finite eta in [-1,1]")
        if derivative not in (0, 1, 2):
            raise ValueError("Implemented derivative orders are zero through two")
        e = eta[..., None]
        weights = np.exp(self.log_weights-2*self.thetas*np.log1p(e*e))
        if derivative:
            a = -4*self.thetas*e/(1+e*e)
            factor = a if derivative == 1 else a*a-4*self.thetas*(1-e*e)/(1+e*e)**2
            weights = weights*factor
        return -.5*np.sum(weights, axis=-1)

    def stage_log_magnitudes(self, eta):
        """Log absolute normalized pressure contribution, without tiny-tail underflow."""
        if not math.isfinite(eta) or abs(eta)>1:
            raise ValueError("Require finite eta in [-1,1]")
        names = list(dict.fromkeys(row["stage"] for row in self.components))
        return {name:float(logsumexp([row["log_weight"]
                    -2*row["theta"]*math.log1p(eta*eta)
                    for row in self.components if row["stage"] == name])-math.log(2.))
                for name in names}

    def eta_coefficients(self, eta, order, digits=60, normalized=False):
        """Ordinary Taylor coefficients, with binary64 quadrature accuracy.

        The recurrence (1+eta^2)*f'=-2*beta*eta*f generates f=(1+eta^2)^(-beta).
        The small stage integrals remain represented even when their addition
        to the total pressure cannot affect the selected working digits.
        """
        with mp.workdps(digits):
            e = mp.mpf(str(eta))
            total = [mp.mpf(0)]*(order+1)
            scale = mp.mpf(0) if normalized else 2*mp.mpf(str(self.schedule.log_P))
            for row in self.components:
                beta = 2*mp.mpf(str(row["theta"]))
                weight = -mp.exp(mp.mpf(str(row["log_weight"]))+scale)/2
                coefficients = [(1+e*e)**(-beta)]
                for n in range(order):
                    previous = coefficients[n-1] if n else 0
                    value = (-2*e*(n+beta)*coefficients[n]
                             -(n-1+2*beta)*previous)/((1+e*e)*(n+1))
                    coefficients.append(value)
                for n, value in enumerate(coefficients):
                    total[n] += weight*value
            return tuple(total)

    def axis_condition_bound(self, j=.02, sigma=1e-4, delta=1., digits=60):
        """An algebraic sufficient bound for B.2, using only A.22 and its sign.

        If |Z_star|<=delta, the lower pressure bound confines eta to a band
        where H_star is bounded away from zero. This assesses B.2 alone.
        """
        if not 0<j<=.05 or not sigma>0 or not delta>0:
            raise ValueError("Require 0<j<=.05, sigma>0 and delta>0")
        with mp.workdps(digits):
            j, sigma, delta = map(lambda x:mp.mpf(str(x)), (j,sigma,delta))
            h = mp.mpf(str(self.schedule.h))
            A, D = mp.mpf(".5")+h, mp.mpf(".5")-h
            umax = 4+j
            geometry_bound = A*(1+2*umax)*umax+4*(D+umax)
            eta_bound = ((delta+geometry_bound)/(mp.mpf("2.5")*A)
                         *mp.exp(-2*mp.mpf(str(self.schedule.log_P))))
            H_lower = j-(D+4)*eta_bound-j*eta_bound**2
            chi_lower = H_lower**2/(H_lower**2+sigma**2) if H_lower>0 else mp.mpf(0)
            return dict(geometry_bound=geometry_bound,eta_band_bound=eta_bound,
                        H_lower=H_lower,chi_lower=chi_lower,
                        condition_sufficient=bool(eta_bound<1 and H_lower>0
                                                  and chi_lower>mp.mpf(".99")))


@lru_cache(maxsize=16)
def pressure_trace(schedule=Schedule(), order=128):
    return PressureTrace(schedule, order)
