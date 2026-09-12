"""D6 stage D: a bounded complex amplitude and targeted inner-exit sampling."""
import csv
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from .axis_bounds import (complex_amplitude_bound,axis_primitive,complex_axis_poles,
                          feature_scales,exit_diagnostics)
from .inner import AxisData,inner_series
from .pressure import Schedule,pressure_trace

ROOT=Path(__file__).resolve().parents[3]
STUDY=ROOT/"code/d6_global_axis"
WRITING=ROOT/"writing/d6_global_axis"
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


def make_data(p,trace,ratio):
    data=AxisData(h=trace.schedule.h,offset_j=p["offset_j"],
        regularization_sigma=p["regularization_sigma"],
        radial_scale_Lambda=ratio*math.exp(2*trace.schedule.log_P))
    bound=complex_amplitude_bound(data,p["log_C_margin"],p["log_C_decimal_places"])
    return replace(data,log_amplitude_C=bound["log_C_upper"]),bound


def generate():
    p=load_parameters()
    trace=pressure_trace(Schedule(**p["pressure_schedule"]),p["pressure_quadrature_order"])
    data,bound=make_data(p,trace,p["Lambda_over_P_squared"])
    (STUDY/"data").mkdir(parents=True,exist_ok=True)
    (WRITING/"figures").mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    digits=p["decimal_digits"]
    with mp.workdps(p["coordinate_digits"]):
        feature=feature_scales(data,trace,p["coordinate_digits"])
        center,width=feature["center"],feature["pressure_width"]
        coordinates=[("global",mp.mpf(str(e))) for e in np.linspace(-1,1,p["global_eta_samples"])]
        coordinates += [("chi",center+mp.mpf(str(c))*feature["chi_width"]) for c in p["chi_offsets"]]
        coordinates += [("critical",center+mp.mpf(str(c))*width) for c in p["critical_offsets"]]
        profiles,exits,amplitudes=[],[],[]
        for group,e in coordinates:
            degree=p["critical_radial_degree"] if group=="critical" else p["radial_degree"]
            series=inner_series(data,e,degree,digits,axis_pressure=trace)
            amplitude=axis_primitive(data,e,digits)
            amplitudes.append(dict(group=group,eta_decimal=mp.nstr(e,80),
                eta_float=float(e),primitive=float(amplitude),
                log_g_over_Lambda=float(amplitude-mp.mpf(data.log_amplitude_C)/data.radial_scale_Lambda)))
            for Y in p["radial_samples"]:
                d=exit_diagnostics(series,Y)
                profiles.append(dict(group=group,eta_decimal=mp.nstr(e,80),Y=Y,radial_degree=degree,
                    Phi=float(d["Phi"]),U=float(d["U"]),p1=float(d["p1"]),
                    Sq_margin=float(d["Sq_margin"]),
                    angular_relative=float(d["angular_relative"]),axial_relative=float(d["axial_relative"]),
                    pressure_relative=float(d["pressure_relative"])))
            d=exit_diagnostics(series,p["exit_Y"],axial_cap=p["axial_lower_bound_cap"])
            exits.append(dict(group=group,eta_decimal=mp.nstr(e,80),
                chi_offset=float((e-center)/feature["chi_width"]),chi=float(d["chi"]),
                p1=float(d["p1"]),scaled_p1=float(d["scaled_p1"]),
                exit_lower=float(d["exit_lower"]) if d["exit_lower"] is not None else None,
                Sq_margin=float(d["Sq_margin"]),
                log10_axial_ratio=float(d["log_axial_ratio"]/mp.log(10)) if d["log_axial_ratio"] is not None else None))
        write_csv("profile_checks.csv",profiles)
        write_csv("exit_checks.csv",exits)
        write_csv("amplitude_real_axis.csv",amplitudes)

        poles=[]
        for root,residue in complex_axis_poles(data,digits):
            poles.append(dict(real_decimal=mp.nstr(mp.re(root),50),imag_decimal=mp.nstr(mp.im(root),50),
                local_real=float((mp.re(root)-center)/feature["chi_width"]),
                local_imag=float(mp.im(root)/feature["chi_width"])))
        write_csv("complex_poles.csv",poles)
        critical=[]
        for ratio in (p["comparison_Lambda_over_P_squared"],p["Lambda_over_P_squared"]):
            datum,_=make_data(p,trace,ratio)
            local=feature_scales(datum,trace,p["coordinate_digits"])
            for c in p["critical_offsets"]:
                e=local["center"]+mp.mpf(str(c))*local["pressure_width"]
                s=inner_series(datum,e,p["critical_radial_degree"],digits,axis_pressure=trace)
                d=exit_diagnostics(s,p["exit_Y"])
                critical.append(dict(Lambda_over_P_squared=ratio,scaled_offset=c,eta_decimal=mp.nstr(e,80),
                    eta_offset_decimal=mp.nstr(e-local["center"],50),
                    scaled_p1=float(d["scaled_p1"]),Sq_margin=float(d["Sq_margin"]),
                    leading_relative_max=float(max(d[k] for k in ("angular_relative","axial_relative","pressure_relative")))))
        write_csv("critical_scale_scan.csv",critical)
        reference=inner_series(data,mp.mpf(".2"),p["degree_reference"],digits,axis_pressure=trace)
        ref=reference.evaluate("4.1")
        convergence=[]
        for degree in p["degree_scan"]:
            d=reference.evaluate("4.1",order=degree)
            convergence.append(dict(radial_degree=degree,Phi_error=float(abs(d["Phi"]-ref["Phi"])),
                angular_relative=float(d["angular_relative"]),pressure_relative=float(d["pressure_relative"])))
        write_csv("degree_convergence.csv",convergence)
        precision=[]
        for ratio,c in ((p["comparison_Lambda_over_P_squared"],.5),(p["Lambda_over_P_squared"],20)):
            datum,_=make_data(p,trace,ratio)
            f=feature_scales(datum,trace,p["coordinate_digits"])
            e=f["center"]+mp.mpf(str(c))*f["pressure_width"]
            a=exit_diagnostics(inner_series(datum,e,12,65,axis_pressure=trace))
            b=exit_diagnostics(inner_series(datum,e,18,90,axis_pressure=trace))
            precision.append(dict(Lambda_over_P_squared=ratio,scaled_offset=c,
                low_degree=12,high_degree=18,low_digits=65,high_digits=90,
                scaled_p1_low=float(a["scaled_p1"]),scaled_p1_high=float(b["scaled_p1"]),
                scaled_p1_relative_change=float(abs((a["scaled_p1"]-b["scaled_p1"])/b["scaled_p1"])),
                Sq_margin_absolute_change=float(abs(a["Sq_margin"]-b["Sq_margin"]))))
        write_csv("precision_checks.csv",precision)
        e=center+20*width
        a=exit_diagnostics(inner_series(data,e,12,digits,axis_pressure=trace))
        b=exit_diagnostics(inner_series(data,e,12,digits,axis_pressure=pressure_trace(
            trace.schedule,p["pressure_reference_order"])))
        sensitivity=[dict(quantity=k,baseline=float(a[k]),reference=float(b[k]),
                          absolute_change=float(abs(a[k]-b[k])))
                     for k in ("scaled_p1","Sq_margin")]
        write_csv("pressure_sensitivity.csv",sensitivity)

    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    rows=sorted(amplitudes,key=lambda r:r["eta_float"])
    axes[0].plot([r["eta_float"] for r in rows],[r["log_g_over_Lambda"] for r in rows],color=TEAL)
    axes[0].set(xlabel="Axial parameter eta",ylabel="log(g) / Lambda",
                title="The amplitude stays in logarithmic form")
    axes[0].ticklabel_format(axis="y",style="plain",useOffset=False)
    strip=float(bound["half_width"])/float(feature["chi_width"])
    axes[1].axhspan(-strip,strip,color=TEAL,alpha=.18,label="Complex neighborhood (local view)")
    near=[r for r in poles if abs(r["local_real"])<2]
    axes[1].scatter([r["local_real"] for r in near],[r["local_imag"] for r in near],
                    color=GOLD,marker="x",s=65,label="Nearest poles of the axis integrand")
    axes[1].set(xlim=(-2,2),ylim=(-1.3,1.3),xlabel="Real displacement / (sigma / H_star')",
                ylabel="Imaginary part / (sigma / H_star')",
                title="The analytic neighborhood avoids the poles")
    axes[1].legend(fontsize=8,loc="center right",bbox_to_anchor=(1,.27))
    for ax in axes:
        ax.title.set_fontsize(10)
        ax.grid(alpha=.15)
    save_figure(fig,"01_logarithmic_amplitude_bound")

    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    for ratio,color in ((p["comparison_Lambda_over_P_squared"],GOLD),(p["Lambda_over_P_squared"],TEAL)):
        rows=[r for r in critical if r["Lambda_over_P_squared"]==ratio]
        for ax,key in zip(axes,("scaled_p1","Sq_margin")):
            ax.plot([r["scaled_offset"] for r in rows],[r[key] for r in rows],"o-",ms=3,color=color,
                    label=f"Lambda / P_star^2 = {ratio:.0e}")
    axes[0].set(ylabel="Lambda times angular shear p1",title="A narrow band reveals negative shear")
    axes[1].set(ylabel="Angular source minus its required lower bound",
                title="Increasing Lambda recovers the sampled sign")
    for ax in axes:
        ax.set_yscale("symlog",linthresh=1)
        ax.axhline(0,color=INK,lw=.8)
        ax.set_xlabel("Displacement from H_star zero / pressure-dependent width")
        ax.title.set_fontsize(9)
        ax.legend(fontsize=8)
        ax.grid(alpha=.15)
    save_figure(fig,"02_hidden_parameter_band")

    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    rows=sorted([r for r in exits if r["group"]=="chi"],key=lambda r:r["chi_offset"])
    axes[0].plot([r["chi_offset"] for r in rows],[r["p1"] for r in rows],"o-",color=GOLD,ms=3,label="Angular term p1")
    axes[0].plot([r["chi_offset"] for r in rows],[r["exit_lower"] for r in rows],"o-",color=TEAL,ms=3,
                 label="Pointwise exit lower bound")
    axes[0].axhline(p["exit_threshold"],color=INK,ls="--",label="Chosen exit threshold 2.2")
    axes[0].set(xlabel="Displacement from H_star zero / (sigma / H_star')",
                ylabel="Exit expression or a stated lower bound",title="Complementary terms pass the sampled exit checks")
    for key,label,color in (("Phi_error","Phi vs degree 24",TEAL),
                            ("angular_relative","Angular leading defect",GOLD),
                            ("pressure_relative","Pressure balance defect",PURPLE)):
        axes[1].semilogy([r["radial_degree"] for r in convergence],[r[key] for r in convergence],"o-",color=color,label=label)
    axes[1].set(xlabel="Radial degree at eta=0.2, Y=4.1",ylabel="Difference or normalized equation defect",
                title="Radial convergence is checked separately")
    axes[1].set_xticks(p["degree_scan"])
    for ax in axes:
        ax.title.set_fontsize(9)
        ax.legend(fontsize=8)
        ax.grid(alpha=.15)
    save_figure(fig,"03_sampled_exit_and_convergence")

    good=[r for r in critical if r["Lambda_over_P_squared"]==p["Lambda_over_P_squared"]]
    bad=[r for r in critical if r["Lambda_over_P_squared"]==p["comparison_Lambda_over_P_squared"]]
    with mp.workdps(p["coordinate_digits"]):
        critical_floats={float(feature["center"]+mp.mpf(str(c))*feature["pressure_width"]) for c in p["critical_offsets"]}
        summary=dict(log_C_upper=bound["log_C_upper"],complex_half_width=float(bound["half_width"]),
            primitive_bound=float(bound["primitive_modulus_bound"]),B16_exact_rational_envelope=True,
            H_zero_decimal=mp.nstr(feature["center"],80),
            pressure_width_decimal=mp.nstr(feature["pressure_width"],50),
            pressure_width=float(feature["pressure_width"]),
            binary64_distinct_critical_coordinates=len(critical_floats),
            critical_coordinate_count=len(p["critical_offsets"]),
            min_sampled_Phi=min(r["Phi"] for r in profiles),
            min_sampled_Sq_margin=min(r["Sq_margin"] for r in profiles),
            min_exit_lower=min(r["exit_lower"] for r in exits if r["exit_lower"] is not None),
            all_sampled_exit_checks_pass=all(r["p1"]>0 and r["exit_lower"] is not None
                and r["exit_lower"]>p["exit_threshold"] for r in exits),
            bad_min_scaled_p1=min(r["scaled_p1"] for r in bad),
            good_min_scaled_p1=min(r["scaled_p1"] for r in good),
            good_min_critical_Sq_margin=min(r["Sq_margin"] for r in good),
            pressure_Sq_margin_change=sensitivity[1]["absolute_change"],
            global_eta_coordinate_count=len(coordinates),profile_sample_count=len(profiles),
            uniform_inner_remainder_or_between_sample_bound=False,
            complete_outer_hierarchy_certified=False,joined_profile=False,first_q_correction=False,
            parameters_sha256=hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest())
    (STUDY/"data/summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    (STUDY/"data/complex-bound.json").write_text(json.dumps({k:str(v) for k,v in bound.items()},indent=2)+"\n",encoding="utf-8")
    return summary


if __name__=="__main__":
    print(json.dumps(generate(),indent=2))
