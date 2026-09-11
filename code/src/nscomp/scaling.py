"""Leading-order scaling surrogates, with SI inputs and explicit assumptions.

These functions evaluate power laws; they do not solve Navier-Stokes. The paper's
construction has 0 < h < 0.01 and additional smallness conditions. The endpoints
accepted here are limiting sensitivity cases, not constructed solutions.
"""

from dataclasses import dataclass
import numpy as np


def positive(value, name):
    """Validate scalar or array inputs without allowing silent NaNs or zeros."""
    array = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(array)) or np.any(array <= 0):
        raise ValueError(f"{name} must contain only finite positive values")
    return array


@dataclass(frozen=True)
class CoreScaling:
    """An anchored characteristic-radius / characteristic-speed surrogate.

    r0_m and u0_m_s refer to the SAME nonzero-speed reference stage. A different
    normalization is represented by a different reference pair, not a multiplier
    that silently breaks the anchor. h=0 is the inverse-law limiting proxy.
    """

    r0_m: float = 1e-3
    u0_m_s: float = 1.0
    h: float = 0.0

    def __post_init__(self):
        positive(self.r0_m, "r0_m")
        positive(self.u0_m_s, "u0_m_s")
        if not np.isfinite(self.h) or not 0 <= self.h <= 0.01:
            raise ValueError("h must lie in [0, 0.01]; endpoints are limit cases")

    @property
    def exponent(self):
        return 1.0 / (1.0 + 2.0 * self.h)

    def radius(self, speed_m_s):
        speed = positive(speed_m_s, "speed_m_s")
        return self.r0_m * np.exp(-self.exponent * np.log(speed / self.u0_m_s))

    def speed(self, radius_m):
        radius = positive(radius_m, "radius_m")
        return self.u0_m_s * np.exp(-np.log(radius / self.r0_m) / self.exponent)

    def remaining_time_ratio(self, speed_m_s):
        """tau/tau0 only; no physical collapse time is calibrated here."""
        ratio = positive(speed_m_s, "speed_m_s") / self.u0_m_s
        return np.exp(-np.log(ratio) / (0.5 + self.h))


def mach_number(speed_m_s, sound_speed_m_s):
    return positive(speed_m_s, "speed") / positive(sound_speed_m_s, "sound speed")


def knudsen_number(mean_free_path_m, gradient_length_m):
    return positive(mean_free_path_m, "mean free path") / positive(
        gradient_length_m, "gradient length"
    )


def number_density(density_kg_m3, molar_mass_kg_mol, avogadro_mol_inv):
    return (positive(density_kg_m3, "density")
            / positive(molar_mass_kg_mol, "molar mass")
            * positive(avogadro_mol_inv, "Avogadro constant"))


def mean_spacing_m(density_kg_m3, molar_mass_kg_mol, avogadro_mol_inv):
    """n^(-1/3), a cubic mean-spacing estimate, NOT a molecular diameter."""
    return number_density(density_kg_m3, molar_mass_kg_mol, avogadro_mol_inv) ** (-1 / 3)


def dimensionless_viscosity(nu_m2_s, length_unit_m, time_unit_s):
    return (positive(nu_m2_s, "kinematic viscosity")
            * positive(time_unit_s, "time unit")
            / positive(length_unit_m, "length unit") ** 2)


def exponent_envelope(speed_m_s, r0_m=1e-3, u0_m_s=1.0, h_max=0.01):
    """Endpoint envelope for the stated exponent interval at a FIXED anchor.

    This is not a statistical confidence interval or a bound on missing profile
    prefactors. min/max also handle speed below the reference speed correctly.
    """
    first = CoreScaling(r0_m, u0_m_s, 0).radius(speed_m_s)
    last = CoreScaling(r0_m, u0_m_s, h_max).radius(speed_m_s)
    return np.minimum(first, last), np.maximum(first, last)
