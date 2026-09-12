"""Explicit physical screening proxies; no reconstructed fluid field or CFD grid."""

from dataclasses import dataclass
import numpy as np
from .scaling import CoreScaling, positive


@dataclass(frozen=True)
class GasReference:
    """Frozen ideal-gas reference with an independently specified mean-free-path convention."""

    temperature_K: float
    pressure_Pa: float
    gas_constant_J_kg_K: float
    heat_capacity_ratio: float
    dynamic_viscosity_Pa_s: float
    mean_free_path_m: float

    def __post_init__(self):
        for name, value in vars(self).items():
            positive(value, name)
        if self.heat_capacity_ratio <= 1:
            raise ValueError("heat_capacity_ratio must exceed one")

    @property
    def density(self):
        return self.pressure_Pa / (self.gas_constant_J_kg_K * self.temperature_K)

    @property
    def sound_speed(self):
        return np.sqrt(self.heat_capacity_ratio * self.gas_constant_J_kg_K * self.temperature_K)

    @property
    def mean_thermal_speed(self):
        return np.sqrt(8 * self.gas_constant_J_kg_K * self.temperature_K / np.pi)

    @property
    def collision_time_proxy(self):
        return self.mean_free_path_m / self.mean_thermal_speed

    @property
    def kinematic_viscosity(self):
        return self.dynamic_viscosity_Pa_s / self.density


def gas_diagnostics(speed_m_s, length_m, gas):
    """U/L is a variation-rate proxy, not the measured symmetric strain tensor."""
    speed = positive(speed_m_s, "speed")
    length = positive(length_m, "variation length")
    return {
        "mach": speed / gas.sound_speed,
        "knudsen_proxy": gas.mean_free_path_m / length,
        "reynolds_proxy": speed * length / gas.kinematic_viscosity,
        "advection_time_proxy_s": length / speed,
        "sound_crossing_time_s": length / gas.sound_speed,
        "diffusion_time_proxy_s": length**2 / gas.kinematic_viscosity,
        "collision_rate_proxy": gas.collision_time_proxy * speed / length,
    }


def kn_crossing_speed(model, beta, mean_free_path_m, kn_screen):
    """Formal crossing; it may precede the nonzero-speed reference stage."""
    radius = (positive(mean_free_path_m, "mean free path")
              / (positive(kn_screen, "Kn screen") * positive(beta, "beta")))
    return model.speed(radius)


def rate_crossing_speed(model, beta, collision_time_s, rate_marker):
    """Solve tau_collision U/(beta*a(U)) = marker in logarithmic form."""
    log_ratio = (np.log(positive(rate_marker, "rate marker"))
                 + np.log(positive(beta, "beta") * model.r0_m)
                 - np.log(positive(collision_time_s, "collision time") * model.u0_m_s))
    return model.u0_m_s * np.exp(log_ratio / (1 + model.exponent))


def carrier_power_length(model, speed_m_s, carrier_over_core_at_reference):
    """Anchored bare Q power only: ell ~ Q^((1+h)/2), U ~ Q^(-1/2-h).

    Missing phase-gradient, envelope, and profile factors are not supplied here.
    This is a scale convention; no factor 2*pi is silently applied.
    """
    ratio = positive(speed_m_s, "speed") / model.u0_m_s
    prefactor = positive(carrier_over_core_at_reference, "carrier/core anchor")
    exponent = (1 + model.h) / (1 + 2 * model.h)
    return model.r0_m * prefactor * np.exp(-exponent * np.log(ratio))


def averaging_cube_count(variation_length_m, mean_spacing_m, fraction):
    """Mean particle count in a cube of side fraction*L; no universal cutoff."""
    frac = positive(fraction, "averaging fraction")
    if np.any(frac >= 1):
        raise ValueError("averaging fraction must be smaller than one")
    return (frac * positive(variation_length_m, "variation length")
            / positive(mean_spacing_m, "mean spacing"))**3
