"""D3 datasets, figures, and small article tables."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from .energy import CoreEnergyScaling
from .residual import HeatExterior, gaussian_terms, gaussian_force_l2

ROOT=Path(__file__).resolve().parents[3]
STUDY=ROOT/"code/d3_residual"
WRITING=ROOT/"writing/d3_residual"
TEAL,GOLD,INK,PURPLE="#147d83","#bc681b","#1d2b3a","#7558a4"


def load_parameters():
    return json.loads((STUDY/"parameters.json").read_text(encoding="utf-8"))


def models(p=None):
    p=p or load_parameters()
    e=p["heat_exterior"]
    return HeatExterior(e["h"],e["amplitude"],e["quadrature_order"]), CoreEnergyScaling(**p["gaussian"])


def write_csv(name,rows):
    with (STUDY/"data"/name).open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig,name):
    for ext in ("png","svg"):
        options={"metadata":{"Date":None}} if ext=="svg" else {}
        fig.savefig(WRITING/"figures"/f"{name}.{ext}",dpi=180,bbox_inches="tight",**options)
    plt.close(fig)


def reference_terms(p=None):
    p=p or load_parameters()
    exterior,toy=models(p)
    tau=p["reference"]["tau"]
    radius=np.sqrt(2*tau*p["reference"]["xi"])
    a,b,U=toy.widths_and_speed(tau)
    return exterior.terms(radius,tau),gaussian_terms(p["reference"]["gaussian_R"]*a,
               p["reference"]["gaussian_Z"]*b,tau,toy)


def term_table(p=None):
    heat,toy=reference_terms(p)
    lines=["| Term | Heat exterior, radial | Heat exterior, angular | Gaussian toy, angular | Gaussian toy, axial |",
           "| --- | --- | --- | --- | --- |"]
    for key,label in [("time","Local acceleration"),("transport","Transport"),
                      ("pressure","Pressure gradient"),("minus_viscosity","Minus viscosity"),
                      ("residual","Sum / required force")]:
        values=[heat[key][0],heat[key][1],toy[key][1],toy[key][2]]
        lines.append("| "+label+" | "+" | ".join(f"{float(v):.4g}" for v in values)+" |")
    return "\n".join(lines)


def generate():
    p=load_parameters()
    exterior,toy=models(p)
    (STUDY/"data").mkdir(parents=True,exist_ok=True)
    (WRITING/"figures").mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    times=np.geomspace(1,p["scan"]["time_min"],p["scan"]["samples"])
    time_rows=[]
    for tau in times:
        radius=np.sqrt(2*tau*p["reference"]["xi"])
        heat=exterior.terms(radius,tau)
        a,b,U=toy.widths_and_speed(tau)
        gaussian=gaussian_terms(a,p["reference"]["gaussian_Z"]*b,tau,toy)
        d=exterior.derivatives(radius,tau)
        time_rows.append(dict(tau=tau,heat_radius=radius,heat_speed=float(d["velocity"]),
            heat_time_angular=float(heat["time"][1]),
            heat_minus_viscosity_angular=float(heat["minus_viscosity"][1]),
            heat_residual_magnitude=float(np.linalg.norm(heat["residual"])),
            heat_normalized_angular=float(heat["component_normalized"][1]),
            heat_laplacian_raw_relative_error=float(abs(d["raw_laplacian"]/d["laplacian"]-1)),
            gaussian_force_magnitude=float(np.linalg.norm(gaussian["residual"])),
            gaussian_force_l2=gaussian_force_l2(tau,toy),
            gaussian_angular_normalized=float(gaussian["component_normalized"][1])))
    write_csv("time_scan.csv",time_rows)
    radius_rows=[]
    tau=p["reference"]["tau"]
    for xi in np.geomspace(p["scan"]["xi_min"],p["scan"]["xi_max"],81):
        r=np.sqrt(2*tau*xi)
        term=exterior.terms(r,tau)
        radius_rows.append(dict(xi=xi,radius=r,tau=tau,velocity=float(exterior.velocity(r,tau)),
            pressure=float(exterior.pressure(r,tau)),
            radial_transport=float(term["transport"][0]),radial_pressure_gradient=float(term["pressure"][0]),
            angular_time=float(term["time"][1]),angular_minus_viscosity=float(term["minus_viscosity"][1]),
            angular_normalized_residual=float(term["component_normalized"][1])))
    write_csv("exterior_radial_scan.csv",radius_rows)
    zeta=p["quadrature"]["test_zeta"]
    with mp.workdps(p["quadrature"]["reference_digits"]):
        h=mp.mpf(str(exterior.h))
        z=mp.mpf(str(zeta))
        H_ref=mp.quad(lambda v: mp.exp(-v)*v**h*(1+z*v)**(-h),[0,1,mp.inf])/mp.gamma(1+h)
        H1_ref=-h*mp.quad(lambda v: mp.exp(-v)*v**(h+1)*(1+z*v)**(-h-1),[0,1,mp.inf])/mp.gamma(1+h)
    convergence=[]
    for order in p["quadrature"]["orders"]:
        chosen=HeatExterior(exterior.h,exterior.amplitude,order)
        terms=chosen.terms(np.sqrt(4/zeta),1)
        convergence.append(dict(order=order,
            profile_relative_error=abs(float(chosen.profile(zeta))/float(H_ref)-1),
            profile_derivative_relative_error=abs(float(chosen.profile(zeta,1))/float(H1_ref)-1),
            angular_normalized_residual=float(terms["component_normalized"][1])))
    write_csv("quadrature_convergence.csv",convergence)
    finite=[]
    r,tau=p["finite_difference"]["radius"],p["finite_difference"]["tau"]
    exact=exterior.derivatives(r,tau)
    for relative in np.geomspace(p["finite_difference"]["relative_step_max"],
                                 p["finite_difference"]["relative_step_min"],45):
        dr,dt=relative*r,relative*tau
        radial=(exterior.velocity(r+dr,tau)-exterior.velocity(r-dr,tau))/(2*dr)
        temporal=-(exterior.velocity(r,tau+dt)-exterior.velocity(r,tau-dt))/(2*dt)
        finite.append(dict(relative_step=relative,
            radial_derivative_relative_error=abs(float(radial/exact["radial"])-1),
            time_derivative_relative_error=abs(float(temporal/exact["time"])-1)))
    write_csv("finite_difference_check.csv",finite)

    heat,gaussian=reference_terms(p)
    keys=["time","transport","pressure","minus_viscosity"]
    labels=["Local\nacceleration","Transport","Pressure\ngradient","Minus\nviscosity"]
    colors=[GOLD,INK,TEAL,PURPLE]
    fig,axes=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
    for ax,terms,component,title in [
        (axes[0],heat,0,"Heat exterior: radial balance"),
        (axes[1],heat,1,"Heat exterior: angular balance"),
        (axes[2],gaussian,1,"Gaussian toy: angular force remains"),
    ]:
        values=[float(terms[k][component]) for k in keys]
        ax.bar(labels,values,color=colors)
        ax.axhline(0,color="#7e8c9b",lw=0.8)
        ax.set(title=title,ylabel="Signed acceleration (dimensionless)")
        ax.tick_params(axis="x",labelsize=8)
        ax.grid(axis="y",alpha=0.12)
        ax.ticklabel_format(axis="y",style="sci",scilimits=(-2,2))
    fig.suptitle("Separate component scales reveal which terms cancel",fontsize=12)
    save_figure(fig,"01_momentum_terms")

    fig,axes=plt.subplots(1,2,figsize=(11.5,4.3),constrained_layout=True)
    point=np.array([row["gaussian_force_magnitude"] for row in time_rows])
    l2=np.array([row["gaussian_force_l2"] for row in time_rows])
    axes[0].loglog(times,point/point[0],color=GOLD,label="Pointwise force: tau^(-3/2)")
    axes[0].loglog(times,l2/l2[0],color=PURPLE,label="L2 force: tau^(-3/4)")
    axes[0].set(title="The toy needs a growing force",xlabel="Time remaining (toward collapse)",
                ylabel="Required force / reference")
    axes[0].invert_xaxis()
    axes[0].legend(fontsize=9)
    heat_relative=np.array([row["heat_normalized_angular"] for row in time_rows])
    axes[1].loglog(times,np.maximum(heat_relative,1e-17),color=TEAL,label="Heat exterior: numerical residual")
    axes[1].loglog(times,[row["gaussian_angular_normalized"] for row in time_rows],
                   color=GOLD,label="Gaussian toy: nonzero residual")
    axes[1].set(title="Angular residual relative to its own terms",xlabel="Time remaining",
                ylabel="Absolute sum / sum of absolute terms")
    axes[1].invert_xaxis()
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=0.15)
    save_figure(fig,"02_required_force")

    fig,axes=plt.subplots(1,2,figsize=(11.5,4.3),constrained_layout=True)
    for key,label,color in [
        ("profile_derivative_relative_error","H' vs high-precision integral",GOLD),
        ("angular_normalized_residual","Angular cancellation",TEAL)]:
        axes[0].semilogy([q["order"] for q in convergence],
                         np.maximum([q[key] for q in convergence],1e-17),
                         "o-",color=color,label=label)
    axes[0].set(title="Quadrature must converge\nbefore interpreting cancellation",
                xlabel="Generalized Laguerre order",ylabel="Relative error / normalized residual")
    axes[0].legend(fontsize=8)
    for key,label,color in [
        ("radial_derivative_relative_error","Radial derivative",TEAL),
        ("time_derivative_relative_error","Small time derivative",GOLD)]:
        axes[1].loglog([q["relative_step"] for q in finite],
                       np.maximum([q[key] for q in finite],1e-17),color=color,label=label)
    axes[1].invert_xaxis()
    axes[1].set(title="A smaller difference step\neventually loses accuracy",
                xlabel="Relative finite-difference step (decreasing)",ylabel="Relative derivative error")
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=0.15)
    save_figure(fig,"03_numerical_checks")
    summary=dict(
        heat_max_angular_normalized_residual=max(row["heat_normalized_angular"] for row in time_rows),
        heat_max_sampled_residual_magnitude=max(row["heat_residual_magnitude"] for row in time_rows),
        heat_max_raw_laplacian_relative_difference=max(row["heat_laplacian_raw_relative_error"] for row in time_rows),
        gaussian_reference_force_magnitude=float(point[0]),
        gaussian_reference_force_l2=float(l2[0]),
        gaussian_force_growth_factor=float(point[-1]/point[0]),
        gaussian_l2_growth_factor=float(l2[-1]/l2[0]),
        best_time_finite_difference_error=min(row["time_derivative_relative_error"] for row in finite),
        last_quadrature_derivative_relative_error=convergence[-1]["profile_derivative_relative_error"],
        full_background_assembled=False,
        parameters_sha256=hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest())
    (STUDY/"data/summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary


if __name__=="__main__":
    print(json.dumps(generate(),indent=2))
