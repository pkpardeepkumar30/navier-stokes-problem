"""Simple reader illustrations; no new Navier–Stokes experiment.

Run from the repository root with the established Python environment.
Only the pulse curve, pressure share, and gas densities reuse numerical
study results. Other graphics are explicitly labelled arithmetic examples
or conceptual diagrams in their article captions.
"""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INK, TEAL, ORANGE, BLUE, PALE = "#203044", "#087e86", "#ad5815", "#6359a6", "#edf5f5"


def setup():
    plt.rcParams.update({"font.size": 12, "axes.titlesize": 14,
                         "axes.labelsize": 12, "svg.fonttype": "none",
                         "axes.spines.top": False, "axes.spines.right": False,
                         "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": INK, "ytick.color": INK,
                         "svg.hashsalt": "navier-stokes-reader-series"})


def save(fig, folder, records):
    target = ROOT / "writing" / folder / "figures"
    target.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "svg"):
        options = {"metadata": {"Date": None}} if suffix == "svg" else {}
        fig.savefig(target / f"reader.{suffix}", dpi=170, bbox_inches="tight", **options)
        if suffix == "svg":
            svg = target / "reader.svg"
            svg.write_text("\n".join(line.rstrip() for line in
                svg.read_text(encoding="utf-8").splitlines()) + "\n", encoding="utf-8")
    plt.close(fig)
    records.append(folder)


def box(ax, x, y, width, height, title, detail=None, color=TEAL):
    ax.add_patch(FancyBboxPatch((x, y), width, height,
        boxstyle="round,pad=0.02,rounding_size=0.05", fc=PALE, ec=color, lw=1.5))
    ax.text(x+width/2, y+height*(.63 if detail else .5), title,
            ha="center", va="center", fontsize=12, weight="bold", color=color)
    if detail:
        ax.text(x+width/2, y+height*.28, detail, ha="center", va="center",
                fontsize=10.5, linespacing=1.4)


def arrow(ax, start, end, color=INK):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                                mutation_scale=15, lw=1.5, color=color))


def load_summary(folder):
    return json.loads((ROOT / "code" / folder / "data/summary.json").read_text(encoding="utf-8"))


def generate():
    setup()
    records = []
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    ax.set(xlim=(-.15, 11.4), ylim=(-.05, 5.4))
    ax.axis("off")
    positions = [(0, 3.25), (3.9, 3.25), (7.8, 3.25),
                 (7.8, .5), (3.9, .5), (0, .5)]
    labels = [
        ("Concentrate a vortex", "Smaller fast region;\nbounded energy"),
        ("Join the regions", "Smooth core, annulus,\nand exterior"),
        ("Find the missing push", "Calculate the remaining\nmomentum balance"),
        ("Prepare shear + pulses", "Fine ripples and waves\nsupply useful transport"),
        ("Repair the new errors", "Waves, mean motion,\npressure, and totals"),
        ("Finish the construction", "Confine the flow;\nkeep the force smooth")]
    for xy, pair in zip(positions, labels):
        box(ax, *xy, 3.25, 1.65, *pair)
    for start, end in [((3.3, 4.08), (3.85, 4.08)), ((7.2, 4.08), (7.75, 4.08)),
                       ((9.42, 3.2), (9.42, 2.2)), ((7.75, 1.33), (7.2, 1.33)),
                       ((3.85, 1.33), (3.3, 1.33))]:
        arrow(ax, start, end)
    ax.text(5.53, .08, "Repeat repairs with controlled errors", ha="center", fontsize=11)
    save(fig, "paper-map", records)

    fig, ax = plt.subplots(figsize=(10.8, 5.7))
    ax.set(xlim=(0, 10), ylim=(0, 6))
    ax.axis("off")
    box(ax, .3, 3.9, 3.9, 1.45, "Calculate the residual", "Include all newly created terms")
    box(ax, 5.8, 3.9, 3.9, 1.45, "Repair the waves", "Adjust oscillations and\naveraged transport")
    box(ax, 5.8, .6, 3.9, 1.45, "Repair the mean field", "Pressure, mean motion,\nand weighted totals")
    box(ax, .3, .6, 3.9, 1.45, "Control the improvement", "Check the full remaining error")
    for start, end in [((4.3, 4.63), (5.7, 4.63)), ((7.75, 3.8), (7.75, 2.15)),
                       ((5.7, 1.32), (4.3, 1.32)), ((2.25, 2.15), (2.25, 3.8))]:
        arrow(ax, start, end)
    ax.text(5, 2.98, "A repair can create\nanother error", ha="center", va="center", fontsize=13)
    save(fig, "completing-the-flow", records)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    x = np.linspace(0, 1, 4097)
    for count, color in [(2, ORANGE), (8, BLUE), (32, TEAL)]:
        phase = (count*x) % 1
        height = np.minimum(phase, 1-phase)/count
        axes[0].plot(x, height, color=color, lw=1.8, label=f"{count} up/down pairs")
        assert abs(height.max()-1/(2*count)) < 1e-12
    axes[0].axhline(0, color=INK, lw=.8)
    axes[0].set(xlabel="Distance across the page", ylabel="Height above the flat line",
                title="Smaller height, the same slopes")
    axes[0].legend(fontsize=10)
    axes[1].bar(["Upward\nslope", "Downward\nslope", "Equal-share\naverage"],
                [1, -1, 0], color=[ORANGE, BLUE, TEAL], width=.58)
    axes[1].plot([2], [0], "o", color=TEAL, ms=9)
    axes[1].axhline(0, color=INK, lw=.8)
    axes[1].set(ylim=(-1.3, 1.3), yticks=[-1, 0, 1], ylabel="Slope",
                title="The average can be zero")
    save(fig, "convex-integration", records)

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.3), constrained_layout=True)
    stages = np.arange(4)
    relative_speed = 2.**stages
    mass = 8.**(-stages)
    relative_energy = mass*relative_speed**2
    assert np.allclose(relative_energy, 2.**(-stages))
    for ax, values, label, color in [(axes[0], relative_speed, "Speed compared with the first stage", ORANGE),
                                    (axes[1], relative_energy, "Energy compared with the first stage", TEAL)]:
        ax.bar(stages, values, color=color, width=.6)
        for i, v in enumerate(values):
            ax.text(i, v+.03*max(values), f"{v:g}", ha="center", fontsize=12)
        ax.set(xticks=stages, xticklabels=["1", "2", "3", "4"], xlabel="Illustrative stage",
               ylabel=label, ylim=(0, max(values)*1.2))
    axes[0].set_title("Speed rises")
    axes[1].set_title("Energy in the shrinking region falls")
    save(fig, "d2_finite_energy", records)

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.5), constrained_layout=True)
    for radius in [.3, .6, .9]:
        axes[0].add_patch(Circle((0, 0), radius, fill=False, ec=TEAL, lw=1.8))
        arrow(axes[0], (radius, 0), (radius, .4*radius), ORANGE)
    axes[0].plot(0, 0, "o", color=INK, ms=5)
    axes[0].text(0, -.12, "Centre", ha="center")
    axes[0].set(xlim=(-1.1, 1.1), ylim=(-1.1, 1.1), aspect="equal",
                title="A rigid rotating disk")
    axes[0].axis("off")
    axes[1].plot([0, 1], [0, 1], color=TEAL, lw=3)
    axes[1].scatter([0, .5, 1], [0, .5, 1], color=ORANGE, s=50, zorder=3)
    axes[1].set(xlabel="Distance from centre (relative)", ylabel="Circular speed (relative)",
                title="Same turn time, shorter path near centre",
                xlim=(-.05, 1.05), ylim=(-.05, 1.1))
    save(fig, "d6_inner_profile", records)

    pressure = load_summary("d6_axis_pressure")
    share = pressure["pressure_initial_branch_fraction"]*100
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.5), constrained_layout=True)
    axes[0].step(range(4), [10, 8, 7, 6.5], where="mid", color=TEAL, lw=3)
    axes[0].scatter(range(4), [10, 8, 7, 6.5], color=ORANGE, s=55, zorder=3)
    axes[0].set(xticks=range(4), xticklabels=["Outside", "After\nring 1", "After\nring 2", "Inside"],
                ylabel="Illustrative pressure", title="Add the pressure drops", ylim=(5.8, 10.6))
    axes[1].barh(["Source\ncalculation"], [share], color=TEAL, label="Initial branch")
    axes[1].barh(["Source\ncalculation"], [100-share], left=[share], color=ORANGE, label="Remaining intervals")
    axes[1].text(share/2, 0, f"{share:.0f}%", color="white", va="center", ha="center", weight="bold", fontsize=16)
    axes[1].text(share+(100-share)/2, 0, f"{100-share:.0f}%", color="white", va="center", ha="center", weight="bold", fontsize=16)
    axes[1].set(xlim=(0, 100), ylim=(-.9, .9), xlabel="Share of centre-pressure magnitude (%)",
                title="The remaining contribution matters")
    axes[1].legend(loc="upper center", bbox_to_anchor=(.5, -.16), fontsize=10)
    save(fig, "d6_axis_pressure", records)

    fig, ax = plt.subplots(figsize=(9.6, 4.3), constrained_layout=True)
    ax.barh(["Missing", "Already accounted for", "Required"], [1, 2, 3],
            color=[ORANGE, BLUE, TEAL], height=.55)
    for i, value in enumerate([1, 2, 3]):
        ax.text(value+.06, i, f"{value} N", va="center", weight="bold")
    ax.set(xlim=(0, 3.5), xlabel="Force on the one-kilogram cart",
           title="The residual is the missing part of the balance")
    save(fig, "d3_residual", records)

    fig, ax = plt.subplots(figsize=(10.5, 4.2), constrained_layout=True)
    ax.axis("off")
    table = ax.table(cellText=[["A", "+2", "+3", "+6"],
                              ["B", "−2", "−3", "+6"],
                              ["Average", "0", "0", "+6"]],
        colLabels=["Portion", "Upward fluctuation", "Horizontal fluctuation", "Product"],
        cellLoc="center", loc="center", colWidths=[.16, .28, .30, .19])
    table.auto_set_font_size(False)
    table.set_fontsize(13)
    table.scale(1, 2.7)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_facecolor(TEAL if r == 0 else (PALE if c != 3 else "#f7eadf"))
        if r == 0:
            cell.get_text().set_color("white")
            cell.get_text().set_fontsize(11)
            cell.get_text().set_weight("bold")
    ax.set_title("The velocities cancel; their products reinforce", pad=15)
    save(fig, "d4_mean_stress", records)

    with (ROOT / "code/d5_shearing_wave/data/time_history.csv").open(encoding="utf-8", newline="") as stream:
        history = list(csv.DictReader(stream))
    # Column names are checked explicitly against the existing study's export.
    time = np.array([float(row["time"]) for row in history])
    amplitude = np.array([float(row["speed"]) for row in history])
    amplitude /= amplitude[0]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5), constrained_layout=True)
    yy = np.linspace(-1, 1, 50)
    for xx in np.linspace(-.7, .7, 5):
        axes[0].plot(np.full_like(yy, xx), yy, color=BLUE, alpha=.55, lw=2)
        axes[0].plot(xx+.6*yy, yy, color=TEAL, lw=2)
    axes[0].plot([], [], color=BLUE, label="Before")
    axes[0].plot([], [], color=TEAL, label="After shearing")
    axes[0].set(xlim=(-1, 1), ylim=(-1, 1), xlabel="Along the layers",
                ylabel="Across the layers", title="Different speeds tilt the pattern")
    axes[0].legend(fontsize=10)
    axes[1].plot(time, amplitude, color=TEAL, lw=2.5)
    axes[1].axhline(1, color=INK, ls="--", lw=1, label="Initial speed")
    axes[1].set(xlabel="Time in the simple shear example",
                ylabel="Disturbance speed / initial speed", title="Temporary growth, then decay")
    axes[1].legend(fontsize=10)
    save(fig, "d5_shearing_wave", records)

    fig, ax = plt.subplots(figsize=(9.8, 4.5), constrained_layout=True)
    ax.bar([1, 3], [-1, 1], color=[ORANGE, TEAL], width=.6)
    ax.axhline(0, color=INK, lw=1)
    ax.text(1, -1.22, "Decrease by 1", ha="center")
    ax.text(3, 1.1, "Increase by 1", ha="center")
    ax.text(4.1, .4, "Total change:\n−1 + 1 = 0\n\nWeighted change:\n1 × (−1) + 3 × 1 = 2",
            ha="left", va="center", fontsize=13, linespacing=1.6,
            bbox={"facecolor":"white", "edgecolor":"none", "pad":7})
    ax.set(xlim=(.2, 6.5), ylim=(-1.5, 1.5), xticks=[1, 3], xlabel="Position",
           ylabel="Correction amount", title="Preserve one total while changing another")
    assert -1+1 == 0 and -1+3 == 2
    save(fig, "d6_moment_correction", records)

    fig, ax = plt.subplots(figsize=(9.8, 4.5), constrained_layout=True)
    x = np.linspace(0, 1, 3001)
    road = lambda v: 1-1.4*np.exp(-((v-.43)/.015)**2)
    sample = np.linspace(0, 1, 5)
    ax.plot(x, road(x), color=TEAL, lw=2.5, label="Illustrative curve")
    ax.scatter(sample, road(sample), color=ORANGE, s=70, zorder=3, label="Coarse samples")
    ax.axhline(0, color=INK, ls="--")
    ax.annotate("A missed dip", xy=(.43, -.4), xytext=(.61, -.33),
                arrowprops={"arrowstyle":"->", "color":INK}, fontsize=12)
    ax.set(xlabel="Position", ylabel="Value being checked", ylim=(-.62, 1.3),
           title="Passing at sample points does not cover the gaps")
    ax.legend(loc="upper right", fontsize=10)
    assert all(road(sample)>0) and road(.43)<0
    save(fig, "d6_global_axis", records)

    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.set(xlim=(0, 10), ylim=(0, 4))
    ax.axis("off")
    box(ax, .2, 2.15, 5.9, 1.1, "1.2344 − 1.2343", "Keep the original information")
    box(ax, .2, .35, 5.9, 1.1, "1.234 − 1.234", "Round first to four significant digits")
    arrow(ax, (6.3, 2.7), (7.25, 2.7))
    arrow(ax, (6.3, .9), (7.25, .9))
    ax.text(8.35, 2.7, "0.0001", fontsize=23, ha="center", va="center", color=TEAL)
    ax.text(8.35, .9, "0", fontsize=23, ha="center", va="center", color=ORANGE)
    save(fig, "d7_conditioning", records)

    fig, ax = plt.subplots(figsize=(11.2, 4.2))
    ax.set(xlim=(-.1, 11.3), ylim=(0, 4))
    ax.axis("off")
    for i, (speed, size) in enumerate([("1 m/s", "1 millimetre"), ("1,000 m/s", "1 micrometre"),
                                      ("1,000,000 m/s", "1 nanometre")]):
        box(ax, i*3.9, .95, 3.25, 2, speed, size)
    arrow(ax, (3.3, 1.95), (3.8, 1.95))
    arrow(ax, (7.2, 1.95), (7.7, 1.95))
    ax.text(5.45, .28, "Every 1,000-fold speed increase makes the radius 1,000 times smaller",
            ha="center", fontsize=12)
    save(fig, "d1_core_size", records)

    fig, ax = plt.subplots(figsize=(10.7, 4.5))
    ax.set(xlim=(0, 11), ylim=(0, 5))
    ax.axis("off")
    ax.add_patch(Rectangle((.7, 1.2), 2.7, 2.7, fc=PALE, ec=TEAL, lw=2))
    for xx in np.linspace(.88, 3.2, 10):
        ax.scatter(np.full(10, xx), np.linspace(1.38, 3.7, 10), s=8, color=TEAL)
    ax.text(2.05, .67, "1,000 molecules on average\nin the larger cube", ha="center", fontsize=12)
    ax.add_patch(Rectangle((8.1, 2.42), .27, .27, fc=PALE, ec=TEAL, lw=2))
    ax.scatter([8.235], [2.555], color=TEAL, s=14)
    ax.text(8.23, 1.5, "About 1 molecule\nin the smaller cube", ha="center", fontsize=12)
    arrow(ax, (4, 2.55), (7.3, 2.55))
    ax.text(5.65, 3.38, "Each side ÷ 10\nVolume ÷ 1,000", ha="center", fontsize=13, linespacing=1.7)
    ax.text(2.05, 4.3, "Side view of a cube", ha="center", fontsize=11)
    ax.text(5.5, .1, "Same density • The large square uses dots as a texture, not a literal count",
            ha="center", fontsize=10, color=INK)
    save(fig, "d9_regime_map", records)

    gas = load_summary("d10_compressible_vortex")
    states = gas["representative_states"][:3]
    fig, ax = plt.subplots(figsize=(9.6, 4.6), constrained_layout=True)
    fractions = [100*s["center_density_ratio"] for s in states]
    ax.bar(range(3), fractions, color=[TEAL, BLUE, ORANGE], width=.58)
    for i, v in enumerate(fractions):
        ax.text(i, v+2, f"{v:.0f}%", ha="center", fontsize=15, weight="bold")
    speed_labels = [f'{states[0]["speed_m_s"]:.1f} m/s'] + [f'{s["speed_m_s"]:.0f} m/s' for s in states[1:]]
    ax.set(xticks=range(3), xticklabels=speed_labels,
           xlabel="Maximum circular speed in the stationary example",
           ylabel="Centre density / outer density (%)", ylim=(0, 113),
           title="The gas centre becomes less dense at higher rotation")
    save(fig, "d10_compressible_vortex", records)

    expected = [a["folder"] for a in json.loads((HERE/"manifest.json").read_text(encoding="utf-8"))["articles"]]
    assert sorted(records) == sorted(expected)
    inputs = ["code/d5_shearing_wave/data/time_history.csv",
              "code/d6_axis_pressure/data/summary.json",
              "code/d10_compressible_vortex/data/summary.json"]
    data = {"classification":"Arithmetic illustrations, conceptual diagrams, or redisplays of existing study results",
            "energy":{"relative_speed":relative_speed.tolist(),
                      "relative_mass":mass.tolist(),"relative_energy":relative_energy.tolist()},
            "pressure_initial_branch_percent":share,
            "gas_centre_density_percent":fractions,
            "source_data_sha256":{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
            "generated_folders":records}
    (HERE/"figure-data.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print(f"Generated {len(records)} reader figures in PNG and SVG.")
    return records


if __name__ == "__main__":
    generate()
