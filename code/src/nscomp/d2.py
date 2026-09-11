"""D2 figures, data, and an offline animation with explicit normalizations."""

import argparse
import base64
import csv
import hashlib
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np
from PIL import Image

from .energy import (
    CoreEnergyScaling, cylindrical_energy_quadrature, gaussian_energy_tail_fraction,
    gaussian_swirl, narrow_gaussian, narrow_gaussian_squared_integral,
)

ROOT = Path(__file__).resolve().parents[3]
STUDY = ROOT / "code/d2_finite_energy"
WRITING = ROOT / "writing/d2_finite_energy"
TEAL, GOLD, INK, PURPLE = "#147d83", "#bc681b", "#1d2b3a", "#7558a4"


def load_parameters():
    return json.loads((STUDY/"parameters.json").read_text(encoding="utf-8"))


def reference_model(p=None):
    values = (p or load_parameters())["reference"]
    return CoreEnergyScaling(**{k: values[k] for k in ("h", "a0", "b0", "u0", "density")})


def write_csv(name, rows):
    with (STUDY/"data"/name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig, name):
    fig.savefig(WRITING/"figures"/f"{name}.png", dpi=180, bbox_inches="tight",
                metadata={"Software": "Navier-Stokes companion"})
    fig.savefig(WRITING/"figures"/f"{name}.svg", bbox_inches="tight",
                metadata={"Date": None})
    plt.close(fig)


def ratio_table(p=None):
    model = reference_model(p)
    lines = ["| Time remaining / reference | Speed / reference | Radius / reference | Axial width / reference | Volume / reference | Energy / reference |",
             "| --- | --- | --- | --- | --- | --- |"]
    for s in (1, 1e-2, 1e-4, 1e-6):
        r = model.ratios(s)
        values = [s, r["speed"], r["radial"], r["axial"], r["volume"], r["energy"]]
        lines.append("| " + " | ".join(f"{float(v):.3g}" for v in values) + " |")
    return "\n".join(lines)


def data_and_summary(p):
    model = reference_model(p)
    times = np.geomspace(1, p["scaling"]["time_ratio_min"], p["scaling"]["samples"])
    rows = []
    for label, h in [("h_zero_limit", 0), ("necessary_ceiling_endpoint", p["scaling"]["necessary_h_ceiling"])]:
        values = CoreEnergyScaling(h=h).ratios(times)
        for index, s in enumerate(times):
            rows.append(dict(case=label, h=h, time_ratio=s,
                             **{key+"_ratio": float(value[index]) for key, value in values.items()}))
    write_csv("core_scaling.csv", rows)
    toy_rows = []
    for exponent in p["one_dimensional_toy"]["exponents"]:
        for width in np.geomspace(1, 1e-6, 121):
            toy_rows.append(dict(exponent=exponent, width=width,
                                 peak=width**(-exponent),
                                 squared_integral=float(narrow_gaussian_squared_integral(width, exponent))))
    write_csv("one_dimensional_energy.csv", toy_rows)
    checks = []
    cutoff = p["quadrature"]["cutoff_in_widths"]
    for s in (1, 1e-4, 1e-8):
        exact = float(model.toy_energy(s))
        for order in p["quadrature"]["orders"]:
            numeric = cylindrical_energy_quadrature(model, s, order, cutoff)
            checks.append(dict(time_ratio=s, quadrature_order=order,
                               numeric_energy=numeric, whole_space_energy=exact,
                               relative_error=abs(numeric/exact-1),
                               truncation_tail_fraction=gaussian_energy_tail_fraction(cutoff)))
    write_csv("quadrature_convergence.csv", checks)
    material_rows = []
    for s in np.geomspace(1, 1e-2, 121):
        r = model.ratios(s)
        material_rows.append(dict(time_ratio=s, parcel_radial_ratio=np.sqrt(s),
                                  parcel_axial_ratio=1/s, parcel_volume_ratio=1,
                                  core_radial_ratio=float(r["radial"]),
                                  core_axial_ratio=float(r["axial"]),
                                  core_volume_ratio=float(r["volume"])))
    write_csv("material_and_core.csv", material_rows)
    summary = {
        "h_baseline": model.h,
        "baseline_is_limiting_toy": True,
        "radial_exponent": 0.5,
        "axial_exponent": 0.5-model.h,
        "speed_exponent": -0.5-model.h,
        "volume_exponent": 1.5-model.h,
        "energy_exponent": model.energy_exponent,
        "gaussian_energy_constant": float(np.e*np.pi**1.5),
        "gaussian_reference_energy": float(model.toy_energy(1)),
        "ratios_at_time_ratio_1e_6": {k: float(v) for k, v in model.ratios(1e-6).items()},
        "max_order64_relative_quadrature_error": max(row["relative_error"] for row in checks if row["quadrature_order"] == 64),
        "quadrature_tail_fraction": gaussian_energy_tail_fraction(cutoff),
        "necessary_ceiling_slenderness_ratio_at_1e_12": float((1e-12)**p["scaling"]["necessary_h_ceiling"]),
        "paper_profile_evaluated": False,
        "total_paper_energy_evaluated": False,
        "parameters_sha256": hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest(),
    }
    (STUDY/"data/summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    return summary


def static_figures(p):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none",
                         "text.color": INK, "axes.labelcolor": INK})
    model = reference_model(p)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    x = np.linspace(-3, 3, 2401)
    for epsilon, color in zip(p["one_dimensional_toy"]["widths"], (TEAL, GOLD, PURPLE)):
        f = narrow_gaussian(x, epsilon, 0.5)
        axes[0].plot(x, f, color=color, label=f"width = {epsilon:g}")
        axes[1].plot(x, f*f, color=color)
    axes[0].set(xlabel="x (arbitrary units)", ylabel="f(x)", title="A taller, narrower function")
    axes[0].legend(fontsize=8)
    axes[1].set(xlabel="x (arbitrary units)", ylabel="f(x) squared", title="Same integrated squared value")
    widths = np.geomspace(1, 1e-6, 100)
    for exponent, color in zip(p["one_dimensional_toy"]["exponents"], (TEAL, GOLD, PURPLE)):
        normalized = widths**(1-2*exponent)
        axes[2].loglog(widths, normalized, color=color, label=f"p = {exponent:g}")
    axes[2].invert_xaxis()
    axes[2].set(xlabel="Width epsilon (shrinking rightward)", ylabel="Integral / initial integral",
                title="Different growth rates, different outcomes")
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.14)
    save_figure(fig, "01_tall_and_narrow")

    times = np.geomspace(1, 1e-12, 241)
    ratios = model.ratios(times)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4), constrained_layout=True)
    for key, label, color, style in [
        ("speed", "Peak speed", GOLD, "-"), ("volume", "Core volume scale", PURPLE, "-"),
        ("energy", "Toy energy", TEAL, "-"), ("radial", "Radius (overlaps energy here)", INK, ":"),
        ("axial", "Axial width (overlaps radius here)", "#7e8c9b", "--"),
    ]:
        axes[0].loglog(times, ratios[key], label=label, color=color, ls=style)
    axes[0].set(title="Peak grows while volume and energy fall",
                xlabel="Time remaining / reference (toward collapse)", ylabel="Value / reference")
    axes[0].invert_xaxis()
    axes[0].legend(fontsize=8, loc="lower left")
    ceiling = p["scaling"]["necessary_h_ceiling"]
    axes[1].plot(-np.log10(times), 100*(1-times**ceiling), color=TEAL, lw=2)
    axes[1].set(title="Slenderness changes very slowly at the necessary ceiling",
                xlabel="Decades closer to the singular time",
                ylabel="Decrease in (radius / axial width), percent")
    axes[1].text(0.03, 0.90, "h = exp(-10): excluded upper endpoint\nActual permitted h must be smaller",
                 transform=axes[1].transAxes, fontsize=9, va="top")
    for ax in axes:
        ax.grid(alpha=0.15)
    fig.suptitle("Paper powers with fixed toy normalization; left panel uses the h -> 0+ limit", fontsize=12)
    save_figure(fig, "02_core_scaling")

    times = np.geomspace(1, 1e-2, 121)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.2), constrained_layout=True)
    for ax, title, values in [
        (axes[0], "Same material parcel: radial squeezing, axial stretching",
         {"Radius": np.sqrt(times), "Axial width": 1/times, "Volume": np.ones_like(times)}),
        (axes[1], "Intense-flow envelope: both widths shrink",
         {"Radius": np.sqrt(times), "Axial width": times**(0.5-model.h), "Volume": times**(1.5-model.h)}),
    ]:
        for (label, value), color, style in zip(values.items(), (TEAL, GOLD, PURPLE), ("-", "--", "-")):
            ax.loglog(times, value, label=label, color=color, ls=style, lw=2)
        ax.invert_xaxis()
        ax.set(title=title, xlabel="Time remaining / reference", ylabel="Value / reference")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.15)
    fig.suptitle("Two distinct objects and two explicitly separate toy models", fontsize=12)
    save_figure(fig, "03_material_and_core")


def animation(p):
    """Render scientific PNG frames, then package them in a local HTML player."""
    model = reference_model(p)
    settings = p["animation"]
    times = np.geomspace(1, settings["time_ratio_min"], settings["frames"])
    count = settings["grid_points"]
    x = model.a0*np.linspace(-3, 3, count)
    z = model.a0*np.linspace(-6, 6, count)
    R, Z = np.linspace(-3, 3, count), np.linspace(-3, 3, count)
    normalized = gaussian_swirl(np.abs(R)[None, :], Z[:, None], 1, 1, 1)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.5), layout="constrained")
    maximum = float(model.ratios(times[-1])["speed"])
    first = np.maximum(gaussian_swirl(np.abs(x)[None, :], z[:, None], model.a0, model.b0, model.u0), 1e-4)
    left = axes[0].imshow(first/model.u0, extent=[-3, 3, -6, 6], origin="lower", aspect="auto",
                          cmap="viridis", norm=LogNorm(1e-4, maximum))
    right = axes[1].imshow(normalized, extent=[-3, 3, -3, 3], origin="lower", aspect="auto",
                           cmap="viridis", norm=Normalize(0, 1))
    axes[0].set(xlabel="x / initial radial width", ylabel="z / initial radial width",
                title="Fixed spatial coordinates")
    axes[1].set(xlabel="x / current radial width", ylabel="z / current axial width",
                title="Similarity coordinates: fixed profile")
    fig.colorbar(left, ax=axes[0], label="Speed / initial peak (fixed log colors)")
    fig.colorbar(right, ax=axes[1], label="Speed / current peak (fixed linear colors)")
    title = fig.suptitle("", fontsize=11)
    frames, rows = [], []
    for index, s in enumerate(times):
        a, b, U = model.widths_and_speed(s)
        field = gaussian_swirl(np.abs(x)[None, :], z[:, None], a, b, U)
        left.set_data(np.maximum(field/model.u0, 1e-4))
        ratios = model.ratios(s)
        title.set_text(f"Prescribed Gaussian swirl, h=0 limit | remaining time ratio {s:.4f}\n"
                       f"Peak / initial = {float(ratios['speed']):.2f}     Energy / initial = {float(ratios['energy']):.3f}")
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", dpi=105)
        payload = buffer.getvalue()
        compressed = io.BytesIO()
        Image.open(io.BytesIO(payload)).save(compressed, format="WEBP", lossless=True, method=4)
        frames.append("data:image/webp;base64,"+base64.b64encode(compressed.getvalue()).decode("ascii"))
        if index in (0, len(times)//2, len(times)-1):
            suffix = {0: "start", len(times)//2: "middle", len(times)-1: "end"}[index]
            (WRITING/"animations"/f"frame_{suffix}.png").write_bytes(payload)
        if index == len(times)-1:
            fig.savefig(WRITING/"figures/04_animation_static.png", dpi=180, bbox_inches="tight")
            fig.savefig(WRITING/"figures/04_animation_static.svg", bbox_inches="tight", metadata={"Date": None})
        rows.append(dict(frame=index, time_ratio=float(s), radial_width=float(a),
                         axial_width=float(b), peak=float(U), energy_ratio=float(ratios["energy"])))
    plt.close(fig)
    write_csv("animation_frames.csv", rows)
    template = (STUDY/"animation-player.html").read_text(encoding="utf-8")
    html = template.replace("__FRAMES_JSON__", json.dumps(frames)).replace("__FPS__", str(settings["fps"]))
    (WRITING/"animations/core_collapse.html").write_text(html, encoding="utf-8")
    (STUDY/"animation-manifest.json").write_text(json.dumps({
        "frames": len(frames), "fps": settings["fps"],
        "time_ratio_range": [float(times[0]), float(times[-1])],
        "frame_sampling": "uniform in log(time ratio)",
        "embedded_lossless_webp_frames": True, "autoplay": False,
        "paper_profile": False, "particle_tracking": False,
        "axis_aspect": "panels stretch axes to fit; read labeled coordinate units",
        "left_color_range": [1e-4, maximum], "right_color_range": [0, 1],
    }, indent=2)+"\n", encoding="utf-8")


def generate(include_animation=False):
    p = load_parameters()
    for directory in (STUDY/"data", WRITING/"figures", WRITING/"animations"):
        directory.mkdir(parents=True, exist_ok=True)
    summary = data_and_summary(p)
    static_figures(p)
    if include_animation:
        animation(p)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--animation", action="store_true")
    args = parser.parse_args()
    print(json.dumps(generate(include_animation=args.animation), indent=2))
