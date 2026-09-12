"""D6 stage C artifacts: pressure trace, scale audit and inner-solver handoff."""
import csv
import hashlib
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from .pressure import Schedule, pressure_trace
from .inner import AxisData, inner_series

ROOT=Path(__file__).resolve().parents[3]
STUDY=ROOT/"code/d6_axis_pressure"
WRITING=ROOT/"writing/d6_axis_pressure"
TEAL,GOLD,INK,PURPLE="#147d83","#bc681b","#1d2b3a","#7558a4"


def load_parameters():
    return json.loads((STUDY/"parameters.json").read_text(encoding="utf-8"))


def write_csv(name,rows):
    with (STUDY/"data"/name).open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig,name):
    for ext in ("png","svg"):
        kwargs={"metadata":{"Date":None}} if ext=="svg" else {}
        fig.savefig(WRITING/"figures"/f"{name}.{ext}",dpi=180,bbox_inches="tight",**kwargs)
    plt.close(fig)


def generate():
    p=load_parameters()
    trace=pressure_trace(Schedule(**p["schedule"]),p["quadrature_order"])
    reference=pressure_trace(trace.schedule,p["reference_order"])
    (STUDY/"data").mkdir(parents=True,exist_ok=True)
    (WRITING/"figures").mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})

    eta=np.linspace(-1,1,p["eta_samples"])
    profiles=[dict(eta=float(e),Pi0_over_P_squared=float(trace.normalized(e)),
                   first_eta_derivative=float(trace.normalized(e,1)),
                   second_eta_derivative=float(trace.normalized(e,2)),
                   initial_branch_only=-2.5/(1+e*e)**2) for e in eta]
    write_csv("axis_pressure.csv",profiles)
    stage_logs=trace.stage_log_magnitudes(0.)
    contributions=[]
    with mp.workdps(45):
        for name,value in stage_logs.items():
            contributions.append(dict(stage=name,eta=0,
                log10_normalized_pressure_magnitude=value/math.log(10),
                normalized_magnitude_decimal=mp.nstr(mp.exp(mp.mpf(str(value))),20)))
    write_csv("pressure_contributions.csv",contributions)
    schedule_rows=[]
    for s in trace.segments:
        schedule_rows.append(dict(stage=s.name,start_log_X_over_XR=s.start_y,
            end_log_X_over_XR=s.start_y+s.length,log_c_over_P_start=s.log_c_start,
            log_c_over_P_end=float(s.shape(s.length,trace.order)[0]),
            kind=s.kind))
    write_csv("radial_schedule.csv",schedule_rows)
    convergence=[]
    for n in p["quadrature_scan"]:
        t=pressure_trace(trace.schedule,n)
        convergence.append(dict(quadrature_order=n,
            max_pressure_difference=float(np.max(np.abs(t.normalized(eta)-reference.normalized(eta)))),
            relative_Q_target_difference=abs(t.Q_target/reference.Q_target-1),
            waiting_length_difference=abs(t.waiting_length-reference.waiting_length),
            max_stage_log_difference=max(abs(t.stage_log_magnitudes(0)[name]-stage_logs_ref)
                for name,stage_logs_ref in reference.stage_log_magnitudes(0).items())))
    write_csv("quadrature_convergence.csv",convergence)
    tail_log=stage_logs["infinite exterior power"]
    tails=[dict(extra_log_radius=L,
                log10_normalized_pressure_remainder=(tail_log-(1+2*trace.schedule.h)*L)/math.log(10))
           for L in p["tail_cut_offsets"]]
    write_csv("exterior_tail.csv",tails)

    hand=p["inner_handoff"]
    data=AxisData(h=trace.schedule.h,offset_j=hand["offset_j"],
        regularization_sigma=hand["regularization_sigma"],amplitude_C=hand["amplitude_C"],
        radial_scale_Lambda=hand["Lambda_over_P_squared"]*math.exp(2*trace.schedule.log_P))
    series=inner_series(data,hand["eta"],max(hand["degrees"]),hand["digits"],axis_pressure=trace)
    handoff=[]
    for degree in hand["degrees"]:
        v=series.evaluate(hand["Y"],order=degree)
        handoff.append(dict(radial_degree=degree,eta=hand["eta"],Y=hand["Y"],
            Phi=float(v["Phi"]),U=float(v["U"]),
            angular_relative=float(v["angular_relative"]),axial_relative=float(v["axial_relative"]),
            pressure_relative=float(v["pressure_relative"])))
    write_csv("inner_handoff.csv",handoff)
    with mp.workdps(hand["digits"]):
        base_value=series.evaluate(hand["Y"])
        reference_value=inner_series(data,hand["eta"],max(hand["degrees"]),hand["digits"],
                                     axis_pressure=reference).evaluate(hand["Y"])
        propagated=[dict(quantity=key,
            absolute_change=float(abs(base_value[key]-reference_value[key])),
            relative_change=float(abs((base_value[key]-reference_value[key])/reference_value[key])))
            for key in ("Phi","U","Pi")]
    write_csv("handoff_quadrature_change.csv",propagated)
    bound=trace.axis_condition_bound(hand["offset_j"],hand["regularization_sigma"],p["B2_delta"])

    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    axes[0].plot(eta,[r["Pi0_over_P_squared"] for r in profiles],color=TEAL,label="Complete ideal-schedule integral")
    axes[0].plot(eta,[r["initial_branch_only"] for r in profiles],color=GOLD,ls="--",label="Initial branch alone")
    axes[0].set(xlabel="Axial similarity parameter eta",ylabel="Axis pressure / P_star squared",
                title="The outer schedule fixes the axis pressure")
    axes[0].legend(fontsize=8)
    axes[1].plot(eta,[r["first_eta_derivative"] for r in profiles],color=TEAL)
    axes[1].axhline(0,color=INK,lw=.7)
    axes[1].set(xlabel="Axial similarity parameter eta",ylabel="Derivative of normalized pressure",
                title="The derivative has the sign required by the source")
    for ax in axes:
        ax.title.set_fontsize(10)
        ax.grid(alpha=.15)
    save_figure(fig,"01_source_axis_pressure")

    fig,ax=plt.subplots(figsize=(10,6.4),constrained_layout=True)
    vals=[r["log10_normalized_pressure_magnitude"] for r in contributions]
    labels=[r["stage"] for r in contributions]
    ax.scatter(vals,range(len(vals)),color=TEAL,s=30)
    ax.set_yticks(range(len(vals)),labels,fontsize=8)
    ax.invert_yaxis()
    ax.set(xlabel="log10 of absolute pressure contribution / P_star squared",
           title="Small late contributions are retained as logarithms")
    ax.grid(axis="x",alpha=.15)
    ax.axvline(math.log10(np.finfo(float).smallest_subnormal),color=GOLD,ls="--",
               label="Smallest positive binary64 value")
    ax.legend(fontsize=8,loc="upper left")
    save_figure(fig,"02_pressure_contribution_scales")

    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    orders=[r["quadrature_order"] for r in convergence]
    for key,label,color in (("max_pressure_difference","Normalized pressure",TEAL),
                            ("relative_Q_target_difference","Terminal datum (relative)",GOLD)):
        axes[0].loglog(orders,[r[key] for r in convergence],"o-",color=color,label=label)
    axes[0].set(title="Quadrature converges to binary64 accuracy",
                xlabel="Gauss-Legendre nodes per interval",ylabel="Difference from 256-node result")
    axes[0].set_xticks(orders,[str(n) for n in orders])
    axes[0].minorticks_off()
    for key,label,color in (("angular_relative","Angular leading balance",TEAL),
                            ("axial_relative","Axial leading balance",GOLD),
                            ("pressure_relative","Radial pressure balance",PURPLE)):
        axes[1].semilogy(hand["degrees"],[r[key] for r in handoff],"o-",color=color,label=label)
    axes[1].set(title="The computed pressure enters the inner recurrence",
                xlabel="Inner radial degree",ylabel="Relative leading balance defect")
    axes[1].set_xticks(hand["degrees"])
    for ax in axes:
        ax.legend(fontsize=8)
        ax.title.set_fontsize(9)
        ax.grid(alpha=.15)
    save_figure(fig,"03_pressure_and_handoff_convergence")

    summary=dict(pressure_center_normalized=float(trace.normalized(0)),
        pressure_initial_branch_fraction=2.5/abs(float(trace.normalized(0))),
        log_P_star=trace.schedule.log_P,P_star=math.exp(trace.schedule.log_P),
        tail_start_log_radius=trace.tail_start_y,
        tail_start_log10_radius=trace.tail_start_y/math.log(10),
        final_tail_log10_contribution=tail_log/math.log(10),
        Q_target=trace.Q_target,waiting_length=trace.waiting_length,
        baseline_max_pressure_difference=convergence[-1]["max_pressure_difference"],
        baseline_relative_Q_difference=convergence[-1]["relative_Q_target_difference"],
        inner_radial_Lambda=data.radial_scale_Lambda,
        final_inner_handoff=handoff[-1],
        handoff_quadrature_change=propagated,
        B2_sufficient=bound["condition_sufficient"],
        B2_eta_band_bound=float(bound["eta_band_bound"]),
        B2_H_lower=float(bound["H_lower"]),B2_chi_lower=float(bound["chi_lower"]),
        full_outer_hierarchy_certified=False,outer_moments_and_stress_cone_computed=False,
        global_inner_amplitude_and_exit_certified=False,joined_profile=False,
        first_q_power_correction_computed=False,
        parameters_sha256=hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest())
    (STUDY/"data/summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary


if __name__=="__main__":
    print(json.dumps(generate(),indent=2))
