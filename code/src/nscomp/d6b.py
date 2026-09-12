"""D6 stage B: compute source-defined local inner equations with validation axis data."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.optimize import brentq
from .inner import AxisData,inner_series

ROOT=Path(__file__).resolve().parents[3]
STUDY=ROOT/"code/d6_inner_profile"
WRITING=ROOT/"writing/d6_inner_profile"
TEAL,GOLD,INK,PURPLE="#147d83","#bc681b","#1d2b3a","#7558a4"
COLORS=[TEAL,GOLD,PURPLE,"#727b83"]


def load_parameters():
    return json.loads((STUDY/"parameters.json").read_text(encoding="utf-8"))


def axis_values(data,eta):
    h=mp.mpf(str(data.h)); A=mp.mpf(".5")+h; D=mp.mpf(".5")-h
    d,L=1-eta*eta,1-2*h*eta*eta
    U=4*eta+mp.mpf(str(data.offset_j))
    H=D*eta+d*U
    pi=-data.pressure_magnitude/(1+eta*eta)**2
    dpi=4*data.pressure_magnitude*eta/(1+eta*eta)**3
    Z=-A*(1-2*eta*U)*U-4*H-d*dpi+4*A*eta*pi
    chi=H*H/(H*H+mp.mpf(str(data.regularization_sigma))**2)
    return U,H,Z,chi,L


def write_csv(name,rows):
    with (STUDY/"data"/name).open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def save_figure(fig,name):
    for ext in ("png","svg"):
        options={"metadata":{"Date":None}} if ext=="svg" else {}
        fig.savefig(WRITING/"figures"/f"{name}.{ext}",dpi=180,bbox_inches="tight",**options)
    plt.close(fig)


def generate():
    p=load_parameters(); data=AxisData(**p["axis_data"])
    digits,degree=p["decimal_digits"],p["radial_degree"]
    (STUDY/"data").mkdir(parents=True,exist_ok=True)
    (WRITING/"figures").mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    ys=np.linspace(0,p["maximum_Y"],p["radial_samples"])
    profiles=[]
    with mp.workdps(digits):
        for eta in p["eta_centers"]:
            series=inner_series(data,eta,degree,digits)
            for Y in ys:
                v=series.evaluate(Y)
                profiles.append(dict(eta=eta,Y=Y,X=Y/data.radial_scale_Lambda,
                    Phi=float(v["Phi"]),U=float(v["U"]),Pi=float(v["Pi"]),
                    F_decimal=mp.nstr(v["F"],25),pressure_increment_decimal=mp.nstr(v["pressure_increment"],25),
                    log10_F=float(mp.log10(abs(v["F"]))),V0_over_X=float(v["V0_over_X"]),
                    angular_relative=float(v["angular_relative"]),
                    axial_relative=float(v["axial_relative"]),pressure_relative=float(v["pressure_relative"])))
        write_csv("inner_profiles.csv",profiles)
        eta,Y=p["degree_scan"]["eta"],p["degree_scan"]["Y"]
        reference=inner_series(data,eta,p["reference_degree"],digits)
        target=reference.evaluate(Y)
        convergence=[]
        for n in p["degree_scan"]["degrees"]:
            v=reference.evaluate(Y,order=n)
            convergence.append(dict(radial_degree=n,eta=eta,Y=Y,
                Phi_error_vs_reference=float(abs(v["Phi"]-target["Phi"])),
                U_error_vs_reference=float(abs(v["U"]-target["U"])),
                angular_relative=float(v["angular_relative"]),axial_relative=float(v["axial_relative"]),
                pressure_relative=float(v["pressure_relative"])))
        write_csv("degree_convergence.csv",convergence)
        comparisons=[]
        e=mp.mpf(str(p["lambda_scan"]["eta"])); y=mp.mpf(str(p["lambda_scan"]["Y"]))
        for lam in p["lambda_scan"]["values"]:
            datum=AxisData(**{**p["axis_data"],"radial_scale_Lambda":lam})
            v=inner_series(datum,float(e),p["lambda_scan"]["degree"],digits).evaluate(y)
            U,H,Z,chi,L=axis_values(datum,e)
            phi_limit=mp.hyper([],[2],-y*chi/2)
            axial_limit=U-y*Z/(2*lam*L)
            comparisons.append(dict(Lambda=lam,eta=float(e),Y=float(y),
                Phi=float(v["Phi"]),Phi_limit=float(phi_limit),Phi_comparison_error=float(abs(v["Phi"]-phi_limit)),
                U=float(v["U"]),U_first_comparison=float(axial_limit),U_comparison_error=float(abs(v["U"]-axial_limit)),
                max_relative_balance=float(max(v[k] for k in ("angular_relative","axial_relative","pressure_relative")))))
        write_csv("large_lambda_comparison.csv",comparisons)
        storage=[]
        for eta in p["eta_centers"]:
            series=inner_series(data,eta,degree,digits)
            v=series.evaluate(p["storage_Y"])
            axis=series.evaluate(0)
            recovered=v["Pi"]-axis["Pi"]
            inc=v["pressure_increment"]
            storage.append(dict(eta=eta,Y=p["storage_Y"],working_digits=digits,
                separate_increment_decimal=mp.nstr(inc,25),
                increment_recovered_from_total_decimal=mp.nstr(recovered,25),
                relative_increment_loss=float(abs((recovered-inc)/inc))))
        write_csv("pressure_storage.csv",storage)
        root=brentq(lambda eta:float(axis_values(data,mp.mpf(str(eta)))[2]),-.01,.01)
        chi_at_root=float(axis_values(data,mp.mpf(str(root)))[3])

    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    for eta,color in zip(p["eta_centers"],COLORS):
        rows=[r for r in profiles if r["eta"]==eta]
        axes[0].plot([r["Y"] for r in rows],[r["Phi"] for r in rows],color=color,label=f"eta = {eta:g}")
        axes[1].plot([r["Y"] for r in rows],[r["U"] for r in rows],color=color,label=f"eta = {eta:g}")
    axes[0].set(title="The normalized swirl stays regular at the axis",ylabel="Phi = F / axis value of F",ylim=(0,1.05))
    axes[1].set(title="The axial profile follows its coupled equation",ylabel="Axial similarity profile U")
    for ax in axes:
        ax.set_xlabel("Rescaled radial coordinate Y = Lambda X")
        ax.legend(fontsize=8)
        ax.title.set_fontsize(10)
        ax.grid(alpha=.15)
    save_figure(fig,"01_regular_inner_profiles")

    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    for key,label,color in (("angular_relative","Angular leading balance",TEAL),
                            ("axial_relative","Axial leading balance",GOLD),
                            ("pressure_relative","Radial pressure balance",PURPLE)):
        axes[0].semilogy([r["radial_degree"] for r in convergence],
                         [max(r[key],1e-70) for r in convergence],"o-",color=color,label=label)
    axes[0].set(ylabel="Absolute balance defect / sum of term magnitudes",
                title="Increasing radial degree reduces the leading defects")
    for key,label,color in (("Phi_error_vs_reference","Normalized swirl vs degree 24",TEAL),
                            ("U_error_vs_reference","Axial velocity vs degree 24",GOLD)):
        axes[1].semilogy([r["radial_degree"] for r in convergence],
                         [max(r[key],1e-70) for r in convergence],"o-",color=color,label=label)
    axes[1].set(ylabel="Absolute difference from higher-degree profile",
                title="Profile values converge as well")
    for ax in axes:
        ax.set_xlabel("Radial polynomial degree (not correction order)")
        ax.set_xticks(p["degree_scan"]["degrees"])
        ax.legend(fontsize=8)
        ax.title.set_fontsize(9)
        ax.grid(alpha=.15)
    save_figure(fig,"02_radial_degree_convergence")

    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    lambdas=np.array([r["Lambda"] for r in comparisons])
    for ax,key,power,title,color in ((axes[0],"Phi_comparison_error",-1,"Swirl approaches the source comparison profile",TEAL),
                                     (axes[1],"U_comparison_error",-2,"Including the first axial term reduces the difference",GOLD)):
        errors=np.array([r[key] for r in comparisons])
        ax.loglog(lambdas,errors,"o-",color=color,label="Computed difference")
        ax.loglog(lambdas,errors[0]*(lambdas/lambdas[0])**power,color=INK,ls="--",label=f"Lambda^{power} reference")
        ax.set(title=title,xlabel="Radial scale parameter Lambda",ylabel="Absolute difference from comparison")
        ax.set_xticks(lambdas,[str(x) for x in lambdas])
        ax.minorticks_off()
        ax.legend(fontsize=8)
        ax.title.set_fontsize(9)
        ax.grid(alpha=.15)
    save_figure(fig,"03_large_lambda_comparison")
    summary=dict(radial_degree=degree,decimal_digits=digits,
        minimum_sampled_Phi=min(r["Phi"] for r in profiles),
        maximum_sampled_angular_relative=max(r["angular_relative"] for r in profiles),
        maximum_sampled_axial_relative=max(r["axial_relative"] for r in profiles),
        maximum_sampled_pressure_relative=max(r["pressure_relative"] for r in profiles),
        highest_degree_scan_defects={k:convergence[-1][k] for k in ("angular_relative","axial_relative","pressure_relative")},
        eta_of_validation_Z_zero=root,chi_at_validation_Z_zero=chi_at_root,
        condition_B2_satisfied=False if chi_at_root<.99 else None,
        pressure_storage_max_relative_increment_loss=max(r["relative_increment_loss"] for r in storage),
        axis_pressure_from_prepared_outer_profile=False,joined_inner_annular_profile=False,
        first_q_power_correction_computed=False,
        scope="source leading inner equations with declared validation axis data; no fully joined paper realization",
        parameters_sha256=hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest())
    (STUDY/"data/summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary


if __name__=="__main__":
    print(json.dumps(generate(),indent=2))
