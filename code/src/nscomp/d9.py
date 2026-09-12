"""D9: anchored scale scenarios, frozen-property screens, and explicit missing physics."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import numpy as np
from .scaling import CoreScaling, mean_spacing_m
from .regimes import (GasReference, gas_diagnostics, kn_crossing_speed,
    rate_crossing_speed, carrier_power_length, averaging_cube_count)
from .d1 import markdown_table

ROOT = Path(__file__).resolve().parents[3]
STUDY = ROOT / "code/d9_regime_map"
WRITING = ROOT / "writing/d9_regime_map"
TEAL, GOLD, INK, PURPLE = "#147d83", "#bc681b", "#1d2b3a", "#7558a4"
COLORS = [TEAL, GOLD, PURPLE, "#727b83"]


def load_parameters():
    return json.loads((STUDY/"parameters.json").read_text(encoding="utf-8"))


def reference(p):
    r = p["reference"]
    model = CoreScaling(r["radius_m"], r["speed_m_s"], r["h"])
    gas = GasReference(**{key: p["air"][key] for key in GasReference.__dataclass_fields__})
    water = p["water"]
    spacing = float(mean_spacing_m(water["density_kg_m3"], water["molar_mass_kg_mol"],
                                  p["constants"]["avogadro_mol_inv"]))
    return model, gas, spacing


def crossing_rows(p):
    model, gas, _ = reference(p)
    mach = p["screens"]["mach"]*gas.sound_speed
    rows = []
    for beta in p["gradient_factors"]:
        kn = float(kn_crossing_speed(model,beta,gas.mean_free_path_m,p["screens"]["gas_knudsen"]))
        rows.append(dict(beta=beta,kn_crossing_speed_m_s=kn,
            core_radius_at_kn_m=float(model.radius(kn)),
            gradient_length_at_kn_m=float(beta*model.radius(kn)),
            mach_at_kn=kn/gas.sound_speed,mach_crossing_speed_m_s=mach,
            first_of_two_screens="Kn" if kn < mach else "Mach",
            kn_already_exceeded_at_reference=kn < model.u0_m_s,
            collision_rate_at_kn=float(gas_diagnostics(kn,beta*model.radius(kn),gas)["collision_rate_proxy"]),
            speed_at_rate_marker_m_s=float(rate_crossing_speed(model,beta,gas.collision_time_proxy,
                                                              p["screens"]["collision_rate_marker"]))))
    return rows


def crossing_table(p):
    return markdown_table(["Assumed L / core radius", "Air Kn screen (m/s)",
        "Air Mach screen (m/s)", "First of these two"], [
        (f'{r["beta"]:g}', f'{r["kn_crossing_speed_m_s"]:.3g}',
         f'{r["mach_crossing_speed_m_s"]:.3g}', r["first_of_two_screens"]+
         ("; already crossed at reference" if r["kn_already_exceeded_at_reference"] else ""))
        for r in crossing_rows(p)])


def write_csv(name, rows):
    with (STUDY/"data"/name).open("w",encoding="utf-8",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig, name):
    for ext in ("png","svg"):
        options = {"metadata":{"Date":None}} if ext=="svg" else {}
        fig.savefig(WRITING/"figures"/f"{name}.{ext}",dpi=180,bbox_inches="tight",**options)
    plt.close(fig)


def generate():
    p = load_parameters()
    model, gas, spacing = reference(p)
    (STUDY/"data").mkdir(parents=True,exist_ok=True)
    (WRITING/"figures").mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    speed = np.geomspace(p["speed_range_m_s"][0],p["speed_range_m_s"][1],p["samples"])
    radius = model.radius(speed)
    crosses = crossing_rows(p)
    write_csv("air_crossings.csv",crosses)
    rows = []
    for beta in p["gradient_factors"]:
        length = beta*radius
        gas_values = gas_diagnostics(speed,length,gas)
        for i,u in enumerate(speed):
            rows.append(dict(speed_m_s=u,beta=beta,core_radius_m=radius[i],
                assumed_variation_length_m=length[i],
                air_mach=gas_values["mach"][i],air_knudsen_proxy=gas_values["knudsen_proxy"][i],
                air_reynolds_proxy=gas_values["reynolds_proxy"][i],
                air_collision_rate_proxy=gas_values["collision_rate_proxy"][i],
                water_mach=u/p["water"]["sound_speed_m_s"],
                water_length_over_mean_spacing=length[i]/spacing,
                water_molecules_in_assumed_averaging_cube=float(averaging_cube_count(length[i],spacing,.1))))
    write_csv("scale_trajectories.csv",rows)
    times = []
    time_speed = np.geomspace(1,1e4,301)
    for beta in p["gradient_factors"]:
        vals = gas_diagnostics(time_speed,beta*model.radius(time_speed),gas)
        for i,u in enumerate(time_speed):
            times.append(dict(speed_m_s=u,beta=beta,collision_time_proxy_s=gas.collision_time_proxy,
                              **{key:values[i] for key,values in vals.items()}))
    write_csv("air_time_scales.csv",times)

    # Each panel uses its own material information; liquid bands are geometric only.
    fig,axes = plt.subplots(1,2,figsize=(11.6,4.9),sharey=True,constrained_layout=True)
    for ax in axes:
        for beta,color in zip(p["gradient_factors"],COLORS):
            ax.loglog(speed,beta*radius,color=color,label=f"L = {beta:g} a")
        ax.set(xlim=(speed[0],speed[-1]),ylim=(1e-11,.03),xlabel="Characteristic speed U (m/s)")
        ax.set_xticks([1,1e2,1e4,1e6])
        ax.set_yticks([1e-10,1e-8,1e-6,1e-4,1e-2])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.grid(alpha=.12)
    axes[0].set(ylabel="Assumed flow-variation length L (m)",title="Air: compressibility and rarefaction screens")
    axes[1].set_title("Water: compressibility screen and molecular rulers")
    for ax,cs in zip(axes,(gas.sound_speed,p["water"]["sound_speed_m_s"])):
        ax.axvspan(p["screens"]["mach"]*cs,speed[-1],color=GOLD,alpha=.075)
        ax.axvline(p["screens"]["mach"]*cs,color=INK,ls="--",lw=1)
        ax.text(.97,.96,"Right of dashed line:\nMach > 0.3",ha="right",va="top",
                transform=ax.transAxes,fontsize=8)
    air_length = gas.mean_free_path_m/p["screens"]["gas_knudsen"]
    axes[0].axhspan(1e-11,air_length,color=PURPLE,alpha=.07)
    axes[0].axhline(air_length,color=PURPLE,ls=":",lw=1.5)
    axes[0].text(.97,.57,"Below dotted line: Kn proxy > 0.01",ha="right",
                 transform=axes[0].transAxes,fontsize=8,color=PURPLE)
    axes[0].legend(loc="lower left",fontsize=8)
    for mult in (1,10):
        axes[1].axhline(mult*spacing,color=PURPLE,ls=":",lw=1)
        axes[1].text(.035,mult*spacing*1.25,f"{mult} water spacing"+("s" if mult>1 else ""),
                     transform=axes[1].get_yaxis_transform(),fontsize=8,color=PURPLE)
    axes[1].axhspan(spacing,10*spacing,color=PURPLE,alpha=.1)
    bohr = p["constants"]["bohr_radius_m"]
    axes[1].axhline(bohr,color=INK,ls=":",lw=.8)
    axes[1].text(.035,bohr*1.3,"Bohr radius (reference only)",
                 transform=axes[1].get_yaxis_transform(),fontsize=8)
    for ax in axes: ax.title.set_fontsize(10)
    save_figure(fig,"01_physical_regime_map")

    fig,axes = plt.subplots(1,2,figsize=(11.6,4.3),constrained_layout=True)
    betas = np.geomspace(.001,1,201)
    axes[0].loglog(betas,kn_crossing_speed(model,betas,gas.mean_free_path_m,p["screens"]["gas_knudsen"]),
                  color=TEAL,label="Air Kn = 0.01 crossing")
    axes[0].axhline(p["screens"]["mach"]*gas.sound_speed,color=GOLD,ls="--",label="Air Mach = 0.3 crossing")
    axes[0].axhspan(.1,model.u0_m_s,color=INK,alpha=.06,label="Earlier than the reference stage")
    axes[0].set(xlabel="Assumed variation length / core radius",ylabel="Formal crossing speed (m/s)",
                ylim=(.1,300),xlim=(.001,1),title="A smaller scale can reverse the ordering")
    axes[0].set_xticks([.001,.01,.1,1],["0.001","0.01","0.1","1"])
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].legend(fontsize=7,loc="lower right")
    for beta,color in zip((1.,.1,.01),COLORS):
        vals = gas_diagnostics(time_speed,beta*model.radius(time_speed),gas)
        axes[1].loglog(time_speed,vals["advection_time_proxy_s"],color=color,label=f"L/U, L = {beta:g} a")
    axes[1].axhline(gas.collision_time_proxy,color=INK,ls="--",label="Mean-free-path / thermal-speed proxy")
    axes[1].axvline(p["screens"]["mach"]*gas.sound_speed,color=GOLD,ls=":",label="Mach = 0.3")
    axes[1].set(xlabel="Characteristic speed U (m/s)",ylabel="Time (s)",
                title="Molecular time and flow time answer another question")
    axes[1].legend(fontsize=7,loc="lower left")
    for ax in axes:
        ax.title.set_fontsize(9)
        ax.grid(alpha=.15)
    save_figure(fig,"02_screen_order_and_times")

    pulse = []
    pulse_speed = np.geomspace(model.u0_m_s,p["pulse"]["maximum_speed_m_s"],301)
    h_values = [0.,p["pulse"]["example_h"],p["pulse"]["necessary_h_ceiling"]]
    for h in h_values:
        m = CoreScaling(model.r0_m,model.u0_m_s,h)
        for alpha in p["pulse"]["carrier_core_anchor_factors"]:
            ell = carrier_power_length(m,pulse_speed,alpha)
            for i,u in enumerate(pulse_speed):
                pulse.append(dict(h=h,carrier_core_anchor=alpha,speed_m_s=u,
                    core_power_length_m=float(m.radius(u)),carrier_power_length_m=ell[i],
                    carrier_over_core=ell[i]/m.radius(u),
                    relative_to_reference_ratio=ell[i]/(alpha*m.radius(u)),
                    status="bare powers only; phase and profile prefactors uncomputed"))
    write_csv("pulse_power_comparison.csv",pulse)
    fig,axes = plt.subplots(1,2,figsize=(11.6,4.2),constrained_layout=True)
    for h,label,color in zip(h_values,("h = 0 limit","h = 1e-6 scenario","h = exp(-10) excluded ceiling"),COLORS):
        m = CoreScaling(model.r0_m,model.u0_m_s,h)
        normalized = carrier_power_length(m,pulse_speed,1)/m.radius(pulse_speed)
        axes[0].semilogx(pulse_speed/model.u0_m_s,100*(1-normalized),color=color,label=label)
    axes[0].set(xlabel="Speed ratio U/U0 (U0 = 1 m/s)",ylabel="Decrease in carrier/core ratio (%)",
                title="Bare exponents give little extra separation here")
    axes[0].legend(fontsize=8)
    m = CoreScaling(model.r0_m,model.u0_m_s,p["pulse"]["example_h"])
    for alpha,color in zip(p["pulse"]["carrier_core_anchor_factors"],COLORS):
        axes[1].loglog(pulse_speed/model.u0_m_s,carrier_power_length(m,pulse_speed,alpha)/m.radius(pulse_speed),
                       color=color,label=f"Assumed initial ratio {alpha:g}")
    axes[1].set(xlabel="Speed ratio U/U0",ylabel="Bare carrier scale / core radius",
                title="Unspecified initial factors still control the ratio",ylim=(5e-4,2))
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.set_xticks([1,1e2,1e4,1e6])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.title.set_fontsize(9)
        ax.grid(alpha=.15)
    save_figure(fig,"03_carrier_scale_audit")

    windows = []
    separation = np.geomspace(1,1000,301)
    fig,ax = plt.subplots(figsize=(7.6,4.3),constrained_layout=True)
    for fraction,color in zip(p["averaging_fractions"],COLORS):
        counts = averaging_cube_count(separation,1.,fraction)
        ax.loglog(separation,counts,color=color,label=f"Cube side = {fraction:g} L")
        for ratio,count in zip(separation,counts):
            windows.append(dict(variation_length_over_spacing=ratio,averaging_fraction=fraction,
                                expected_count_in_cube=count))
    write_csv("averaging_windows.csv",windows)
    ax.axhline(1,color=INK,ls=":",label="One molecule on average")
    ax.axhline(p["illustrative_molecule_count"],color=INK,ls="--",label="1000 molecules: illustration, not a cutoff")
    ax.set(xlabel="Flow-variation length / mean molecular spacing",ylabel="Mean molecules in averaging cube",
           title="The averaging volume must fit inside the changing flow",ylim=(1e-5,1e7))
    ax.legend(fontsize=8,loc="upper left")
    ax.grid(alpha=.15)
    save_figure(fig,"04_averaging_scale_separation")

    ceiling_model = CoreScaling(model.r0_m,model.u0_m_s,p["pulse"]["necessary_h_ceiling"])
    endpoint = p["pulse"]["maximum_speed_m_s"]
    summary = dict(air_sound_speed_m_s=float(gas.sound_speed),
        air_mean_thermal_speed_m_s=float(gas.mean_thermal_speed),
        air_collision_time_proxy_s=float(gas.collision_time_proxy),
        water_mean_spacing_m=spacing,
        air_mach_screen_speed_m_s=float(p["screens"]["mach"]*gas.sound_speed),
        water_mach_screen_speed_m_s=float(p["screens"]["mach"]*p["water"]["sound_speed_m_s"]),
        beta_for_simultaneous_air_screens=float(air_length/model.radius(p["screens"]["mach"]*gas.sound_speed)),
        air_core_radius_at_mach_screen_m=float(model.radius(p["screens"]["mach"]*gas.sound_speed)),
        water_core_radius_at_mach_screen_m=float(model.radius(p["screens"]["mach"]*p["water"]["sound_speed_m_s"])),
        ceiling_bare_carrier_ratio_decrease_percent=float(100*(1-carrier_power_length(ceiling_model,endpoint,1)/ceiling_model.radius(endpoint))),
        formal_endpoint_speed_m_s=endpoint,
        actual_paper_gradient_evaluated=False,actual_paper_pulse_wavelength_evaluated=False,
        liquid_relaxation_evaluated=False,full_dimensional_profile_calibration=False,
        parameters_sha256=hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest())
    (STUDY/"data/summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary


if __name__ == "__main__":
    print(json.dumps(generate(),indent=2))
