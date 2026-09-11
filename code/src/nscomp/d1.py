"""D1 data, figures and tabular presentation from one parameter file."""

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .scaling import CoreScaling, exponent_envelope, mean_spacing_m, number_density

ROOT = Path(__file__).resolve().parents[3]
STUDY = ROOT / "code" / "d1_core_size"
WRITING = ROOT / "writing" / "d1_core_size"
TEAL, GOLD, INK, PURPLE = "#147d83", "#bc681b", "#1d2b3a", "#7558a4"


def load_parameters(path=None):
    return json.loads(Path(path or STUDY / "parameters.json").read_text(encoding="utf-8"))


def derived(p):
    a, w, c = p["air"], p["water"], p["constants"]
    density_air = a["pressure_Pa"] / (a["gas_constant_J_kg_K"] * a["temperature_K"])
    sound_air = np.sqrt(a["heat_capacity_ratio"] * a["gas_constant_J_kg_K"] * a["temperature_K"])
    thermal_speed = np.sqrt(8*a["gas_constant_J_kg_K"]*a["temperature_K"]/np.pi)
    return {
        "air_density_kg_m3": density_air,
        "air_sound_speed_m_s": float(sound_air),
        "air_kinematic_viscosity_m2_s": a["dynamic_viscosity_Pa_s"] / density_air,
        "air_mean_thermal_speed_m_s": float(thermal_speed),
        "air_collision_time_proxy_s": a["mean_free_path_m"] / float(thermal_speed),
        "air_mean_spacing_m": (c["boltzmann_J_K"]*a["temperature_K"]/a["pressure_Pa"])**(1/3),
        "water_kinematic_viscosity_m2_s": w["dynamic_viscosity_Pa_s"] / w["density_kg_m3"],
        "water_mean_spacing_m": float(mean_spacing_m(w["density_kg_m3"], w["molar_mass_kg_mol"], c["avogadro_mol_inv"])),
        "water_number_density_m3": float(number_density(w["density_kg_m3"], w["molar_mass_kg_mol"], c["avogadro_mol_inv"])),
        "necessary_h_ceiling": float(np.exp(p["sensitivity"]["necessary_h_ceiling_exponent"])),
    }


def reference_model(p):
    r = p["reference"]
    return CoreScaling(r["radius_m"], r["speed_m_s"], r["h"])


def speed_rows(p):
    model, d = reference_model(p), derived(p)
    c = p["constants"]["speed_of_light_m_s"]
    targets = [
        ("Air Ma=0.3 screen", 0.3*d["air_sound_speed_m_s"], "air"),
        ("Air Ma=1", d["air_sound_speed_m_s"], "air"),
        ("Air Ma=5 marker", 5*d["air_sound_speed_m_s"], "air"),
        ("Water Ma=0.3 screen", 0.3*p["water"]["sound_speed_m_s"], "water"),
        ("Water Ma=1 marker", p["water"]["sound_speed_m_s"], "water"),
        ("1 km/s", 1000.0, "none"),
        ("0.01c (formal)", 0.01*c, "none"),
        ("0.1c (formal)", 0.1*c, "none"),
    ]
    rows = []
    for label, speed, fluid in targets:
        lo, hi = exponent_envelope(speed, model.r0_m, model.u0_m_s, d["necessary_h_ceiling"])
        broad_lo, broad_hi = exponent_envelope(speed, model.r0_m, model.u0_m_s)
        radius = float(model.radius(speed))
        rows.append(dict(
            marker=label, fluid=fluid, speed_m_s=speed, radius_m=radius,
            radius_nm=radius*1e9, diameter_m=2*radius,
            necessary_h_envelope_lower_m=float(lo), necessary_h_envelope_upper_m=float(hi),
            introductory_h_envelope_upper_m=float(broad_hi),
            relative_time_remaining=float(model.remaining_time_ratio(speed)),
            air_mach=speed/d["air_sound_speed_m_s"],
            water_mach=speed/p["water"]["sound_speed_m_s"],
            air_kn_proxy=p["air"]["mean_free_path_m"]/(p["reference"]["gradient_length_over_radius"]*radius),
            water_spacing_ratio=radius/d["water_mean_spacing_m"],
            interpretation="anchored power-law extrapolation; frozen reference properties",
        ))
    return rows


def length_rows(p):
    model, d = reference_model(p), derived(p)
    beta = p["reference"]["gradient_length_over_radius"]
    markers = [
        ("Air Kn=0.01 radius proxy", p["air"]["mean_free_path_m"]/(p["screening"]["knudsen"]*beta), "gas screen"),
        ("Air mean free path", p["air"]["mean_free_path_m"], "length ruler"),
        ("Air mean spacing", d["air_mean_spacing_m"], "length ruler"),
        ("10 water mean spacings", 10*d["water_mean_spacing_m"], "geometric marker, not threshold"),
        ("3 water mean spacings", 3*d["water_mean_spacing_m"], "geometric marker, not threshold"),
        ("1 water mean spacing", d["water_mean_spacing_m"], "geometric marker, not threshold"),
        ("1 Bohr radius", p["constants"]["bohr_radius_m"], "atomic reference, not atom diameter"),
    ]
    return [dict(marker=label, radius_m=radius, radius_nm=radius*1e9,
                 formal_speed_m_s=float(model.speed(radius)),
                 fraction_of_c=float(model.speed(radius))/p["constants"]["speed_of_light_m_s"],
                 interpretation=kind) for label, radius, kind in markers]


def gradient_rows(p):
    model, d = reference_model(p), derived(p)
    ma_speed = p["screening"]["mach"] * d["air_sound_speed_m_s"]
    rows = []
    for beta in p["sensitivity"]["gradient_factors"]:
        radius = p["air"]["mean_free_path_m"] / (p["screening"]["knudsen"]*beta)
        kn_speed = float(model.speed(radius))
        rows.append(dict(beta=beta, kn_screen_radius_m=radius, kn_screen_speed_m_s=kn_speed,
                         mach_screen_speed_m_s=ma_speed,
                         first_of_two_screens="Kn" if kn_speed < ma_speed else "Mach",
                         note="formal crossing; if below U0 the screen is already crossed at the reference"))
    return rows


def sensitivity_rows(p):
    base, d = reference_model(p), derived(p)
    target = .01*p["constants"]["speed_of_light_m_s"]
    rows = []
    groups = [("h", p["sensitivity"]["h_values"] + [d["necessary_h_ceiling"]]),
              ("radius_factor", p["sensitivity"]["radius_factors"]),
              ("speed_factor", p["sensitivity"]["speed_factors"]),
              ("normalization_factor", p["sensitivity"]["normalization_factors"])]
    for parameter, values in groups:
        for value in values:
            r0, u0, h = base.r0_m, base.u0_m_s, base.h
            if parameter == "h": h = value
            elif parameter in ("radius_factor", "normalization_factor"): r0 *= value
            else: u0 *= value
            radius = float(CoreScaling(r0, u0, h).radius(target))
            rows.append(dict(parameter=parameter, value=value, effective_r0_m=r0, u0_m_s=u0,
                             h=h, speed_m_s=target, radius_m=radius, radius_nm=radius*1e9,
                             radius_over_baseline=radius/float(base.radius(target)),
                             note="one factor at a time; normalization_factor is degenerate with radius_factor"))
    return rows


def markdown_table(headers, rows):
    def cell(x):
        return str(x).replace("|", "\\|").replace("\n", " ")
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join("---" for _ in headers) + " |"]
                     + ["| " + " | ".join(cell(x) for x in row) + " |" for row in rows])


def format_length(m):
    for factor, unit in [(1e-3, "mm"), (1e-6, "µm"), (1e-9, "nm"), (1e-12, "pm")]:
        if m >= factor*.999:
            return f"{m/factor:.3g} {unit}"
    return f"{m:.3g} m"


def speed_table(p):
    return markdown_table(["Speed marker", "Speed (m/s)", "Core radius, inverse-law limit"],
                          [(r["marker"], f'{r["speed_m_s"]:.4g}', format_length(r["radius_m"])) for r in speed_rows(p)])


def length_table(p):
    return markdown_table(["Length used as radius", "Length", "Formal speed (m/s)"],
                          [(r["marker"], format_length(r["radius_m"]), f'{r["formal_speed_m_s"]:.3g}') for r in length_rows(p)])


def fluid_table(p):
    d = derived(p)
    return markdown_table(["Reference property", "Air", "Liquid water"], [
        ("Temperature / pressure", "296.15 K / 101.3 kPa", "298.15 K / 100 kPa"),
        ("Density (kg/m³)", f'{d["air_density_kg_m3"]:.4g} (ideal gas)', f'{p["water"]["density_kg_m3"]:.6g}'),
        ("Sound speed (m/s)", f'{d["air_sound_speed_m_s"]:.4g} (ideal gas)', f'{p["water"]["sound_speed_m_s"]:.6g}'),
        ("Kinematic viscosity (m²/s)", f'{d["air_kinematic_viscosity_m2_s"]:.4g}', f'{d["water_kinematic_viscosity_m2_s"]:.4g}'),
        ("Mean spacing n^(-1/3)", format_length(d["air_mean_spacing_m"]), format_length(d["water_mean_spacing_m"])),
        ("Mean free path", format_length(p["air"]["mean_free_path_m"]), "Not assigned a dilute-gas value"),
    ])


def save_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": INK, "text.color": INK,
                         "axes.titleweight": "bold", "axes.titlepad": 13,
                         "svg.fonttype": "none", "svg.hashsalt": "nscomp-d1",
                         "savefig.facecolor": "white"})


def save_figure(fig, name, directory):
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / f"{name}.png", dpi=190, bbox_inches="tight")
    fig.savefig(directory / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def make_figures(p, directory):
    style()
    d, base = derived(p), reference_model(p)
    h_cap = d["necessary_h_ceiling"]
    ratios = np.geomspace(1, 1e8, 350)
    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.5), layout="constrained")
    for h, color, label in [(0, INK, "Inverse-law limit"), (0.01, GOLD, "Introductory ceiling h → 0.01")]:
        model = CoreScaling(1, 1, h)
        ax[0].loglog(ratios, model.radius(ratios), color=color, lw=2, label=label)
        ax[1].semilogx(ratios, model.radius(ratios)*ratios, color=color, lw=2)
    ax[1].semilogx(ratios, CoreScaling(1, 1, h_cap).radius(ratios)*ratios, color=TEAL, lw=2,
                   label="Necessary ceiling from Lemma 4.8")
    ax[0].set(title="A. Radius falls almost inversely with speed", xlabel="Speed ratio U / U₀", ylabel="Radius ratio r / r₀")
    ax[1].set(title="B. How much can the exponent change it?", xlabel="Speed ratio U / U₀", ylabel="Radius / inverse-law estimate")
    ax[0].legend(loc="lower left", fontsize=8)
    ax[1].legend(loc="upper left", fontsize=8)
    for a in ax: a.grid(which="major", alpha=.18)
    save_figure(fig, "01_dimensionless_scaling", directory)

    speeds = np.geomspace(1, .1*p["constants"]["speed_of_light_m_s"], 450)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 6.1), sharey=True, layout="constrained")
    for ax, fluid, sound in zip(axes, ["air", "water"], [d["air_sound_speed_m_s"], p["water"]["sound_speed_m_s"]]):
        screen = p["screening"]["mach"]*sound
        ax.axvspan(screen, speeds[-1]*1.2, color="#f1f2f4", zorder=0)
        ax.loglog(speeds, base.radius(speeds), color=INK, ls="--", lw=2)
        before = speeds <= screen
        if fluid == "air":
            kn_radius = p["air"]["mean_free_path_m"]/(p["screening"]["knudsen"]*p["reference"]["gradient_length_over_radius"])
            before &= base.radius(speeds) >= kn_radius
            ax.axhline(kn_radius, color=GOLD, lw=1.4)
            ax.text(.98, kn_radius*1.35, "Kn=0.01, Lgrad=r", ha="right", color=GOLD, fontsize=8,
                    transform=ax.get_yaxis_transform())
        ax.loglog(speeds[before], base.radius(speeds[before]), color=TEAL, lw=3)
        spacing = d["water_mean_spacing_m"]
        atom = p["constants"]["bohr_radius_m"]
        ax.axhspan(spacing, 3*spacing, color=TEAL, alpha=.12)
        ax.axhspan(atom, 2*atom, color=PURPLE, alpha=.13)
        ax.text(.03, 4.2*spacing, "1–3 water mean spacings (reference ruler)", fontsize=8, color=TEAL, transform=ax.get_yaxis_transform())
        ax.text(.03, atom/1.8, "1–2 Bohr radii (atomic reference ruler)", fontsize=8, color=PURPLE, transform=ax.get_yaxis_transform())
        ax.axvline(screen, color=GOLD, lw=1, ls=":")
        ax.text(screen*1.3, 2.4e-3, "Beyond Ma=0.3 screen", fontsize=8, color="#606b79")
        for frac, label in [(.01, "0.01c"), (.1, "0.1c")]:
            u = frac*p["constants"]["speed_of_light_m_s"]
            r = float(base.radius(u))
            ax.scatter([u], [r], s=28, color=INK, zorder=5)
            offset, alignment = ((8, 12), "left") if frac == .01 else ((-10, -20), "right")
            ax.annotate(label, (u, r), xytext=offset, textcoords="offset points", ha=alignment, fontsize=8)
        ax.set(xlabel="Characteristic speed U (m/s)", xlim=(1, speeds[-1]*1.2), ylim=(1e-11, 1e-2))
        state = "23 °C, 101.3 kPa" if fluid == "air" else "25 °C, 100 kPa"
        ax.set_title(f"{fluid.capitalize()} reference: {state}")
        ax.grid(which="major", alpha=.15)
    axes[0].set_ylabel("Characteristic core radius r (m)")
    fig.supxlabel("Illustrative anchor: r₀=1 mm, U₀=1 m/s. Dashed curves are formal extrapolations; shaded length bands are not universal cutoffs.", fontsize=8)
    save_figure(fig, "02_core_size_and_physical_scales", directory)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), layout="constrained")
    factors = np.geomspace(.1, 10, 120)
    target = .01*p["constants"]["speed_of_light_m_s"]
    r_target_nm = float(base.radius(target))*1e9
    axes[0].loglog(factors, factors*r_target_nm, color=TEAL, lw=2, label="Either reference factor, inverse law")
    speed_sweep = [float(CoreScaling(base.r0_m, base.u0_m_s*f, .01).radius(target))*1e9 for f in factors]
    axes[0].loglog(factors, speed_sweep, color=GOLD, ls="--", lw=1.8,
                   label="Speed factor, introductory h→0.01")
    axes[0].scatter([1], [r_target_nm], color=INK, zorder=5)
    axes[0].axhline(d["water_mean_spacing_m"]*1e9, color=PURPLE, ls=":", label="Water mean spacing")
    axes[0].set(title="A. Calibration dominates the length estimate", xlabel="Factor multiplying reference radius or speed", ylabel="Formal radius at 0.01c (nm)")
    axes[0].legend(fontsize=8, loc="upper left")
    beta = np.geomspace(.1, 10, 120)
    radii = p["air"]["mean_free_path_m"]/(p["screening"]["knudsen"]*beta)
    axes[1].loglog(beta, base.speed(radii), color=TEAL, lw=2, label="Air Kn=0.01 crossing")
    axes[1].axhline(.3*d["air_sound_speed_m_s"], color=GOLD, lw=2, label="Air Ma=0.3 crossing")
    axes[1].set(title="B. The first screen depends on gradient length", xlabel="Assumed β in Lgrad = βr", ylabel="Formal crossing speed (m/s)")
    axes[1].legend(fontsize=8, loc="upper left")
    for ax in axes: ax.grid(which="major", alpha=.18)
    save_figure(fig, "03_calibration_and_gradient_sensitivity", directory)


def generate(parameter_path=None):
    p = load_parameters(parameter_path)
    model, d = reference_model(p), derived(p)
    data_dir, figure_dir = STUDY / "data", WRITING / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    WRITING.mkdir(parents=True, exist_ok=True)
    tables = {"speed_thresholds": speed_rows(p), "length_thresholds": length_rows(p),
              "sensitivity": sensitivity_rows(p), "gradient_sensitivity": gradient_rows(p)}
    for name, rows in tables.items(): save_csv(data_dir / f"{name}.csv", rows)
    grid = np.geomspace(1, .1*p["constants"]["speed_of_light_m_s"], 450)
    curves = []
    for h in sorted(set(p["sensitivity"]["h_values"] + [d["necessary_h_ceiling"]])):
        m = CoreScaling(model.r0_m, model.u0_m_s, h)
        for u, r in zip(grid, m.radius(grid)):
            curves.append(dict(h=h, speed_m_s=float(u), radius_m=float(r), r0_m=m.r0_m, u0_m_s=m.u0_m_s))
    save_csv(data_dir / "radius_speed_curves.csv", curves)
    # Also export the exact plotting abscissas for dimensionless and sensitivity panels.
    q = np.geomspace(1, 1e8, 350)
    save_csv(data_dir / "dimensionless_curves.csv", [
        dict(speed_ratio=float(x), radius_ratio=float(CoreScaling(1, 1, h).radius(x)), h=h)
        for h in [0, .01, d["necessary_h_ceiling"]] for x in q])
    factors = np.geomspace(.1, 10, 120)
    target = .01*p["constants"]["speed_of_light_m_s"]
    target_r = float(model.radius(target))
    save_csv(data_dir / "continuous_sensitivity.csv", [
        dict(factor=float(f), radius_vary_r0_m=float(f*target_r),
             radius_vary_u0_introductory_ceiling_m=float(CoreScaling(model.r0_m, model.u0_m_s*f, .01).radius(target)),
             air_kn_crossing_speed_m_s=float(model.speed(p["air"]["mean_free_path_m"]/(p["screening"]["knudsen"]*f))))
        for f in factors])
    make_figures(p, figure_dir)
    summary = dict(d)
    summary.update({
        "radius_at_0_01c_m": target_r,
        "radius_at_0_1c_m": float(model.radius(.1*p["constants"]["speed_of_light_m_s"])),
        "exponent_ceiling_relative_effect_at_0_01c": float(CoreScaling(model.r0_m, model.u0_m_s, d["necessary_h_ceiling"]).radius(target)/target_r-1),
        "air_reference_reynolds": model.r0_m*model.u0_m_s/d["air_kinematic_viscosity_m2_s"],
        "water_reference_reynolds": model.r0_m*model.u0_m_s/d["water_kinematic_viscosity_m2_s"],
        "inverse_law_proxy_not_constructed_solution": True,
        "profile_viscosity_compatibility_established": False,
        "liquid_relaxation_evaluated": False,
        "source_sha256": next(s["sha256"] for s in json.loads((STUDY/"source-manifest.json").read_text(encoding="utf-8"))["sources"] if s["id"] == "openai-ns"),
        "parameters_sha256": hashlib.sha256(Path(parameter_path or STUDY/"parameters.json").read_bytes()).hexdigest(),
    })
    (data_dir / "summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    return summary


if __name__ == "__main__":
    print(json.dumps(generate(), indent=2))
