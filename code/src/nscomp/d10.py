"""D10 stationary compressible/incompressible vortex comparison."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import numpy as np
from scipy.optimize import minimize_scalar
from .compressible import IsentropicVortex, integrate_density, density_deficit_mach, local_gas_audit
from .regimes import GasReference
from .d1 import markdown_table

ROOT = Path(__file__).resolve().parents[3]
STUDY = ROOT/"code/d10_compressible_vortex"
WRITING = ROOT/"writing/d10_compressible_vortex"
TEAL,GOLD,INK,PURPLE = "#147d83","#bc681b","#1d2b3a","#7558a4"
COLORS = [TEAL,GOLD,PURPLE,"#727b83"]


def load_parameters():
    return json.loads((STUDY/"parameters.json").read_text(encoding="utf-8"))


def reference_gas(p):
    return GasReference(**{k:p["air"][k] for k in GasReference.__dataclass_fields__})


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


def observable_rows(p):
    gas=reference_gas(p)
    rows=[]
    for m in p["profile_mach_values"]:
        v=IsentropicVortex(m,gas.heat_capacity_ratio)
        s=v.state(0)
        result=minimize_scalar(lambda r:-float(v.state(r)["local_mach"]),bounds=(0,3),method="bounded")
        rows.append(dict(peak_mach_using_exterior_sound=m,speed_m_s=float(m*gas.sound_speed),
            center_density_ratio=float(s["density_ratio"]),
            center_temperature_K=float(gas.temperature_K*s["temperature_ratio"]),
            center_compressible_pressure_ratio=float(s["pressure_ratio"]),
            center_incompressible_pressure_ratio=float(s["incompressible_pressure_ratio"]),
            maximum_local_mach=float(-result.fun)))
    return rows


def observable_table(p):
    return markdown_table(["U / exterior sound speed","Center density / exterior density",
        "Center gas pressure / exterior pressure","Constant-density pressure / exterior pressure"],
        [(f'{r["peak_mach_using_exterior_sound"]:g}',f'{r["center_density_ratio"]:.4f}',
          f'{r["center_compressible_pressure_ratio"]:.4f}',f'{r["center_incompressible_pressure_ratio"]:.4f}')
         for r in observable_rows(p)])


def conservative_residual(v,x,y,step):
    derivatives=[]
    for axis in (0,1):
        values=[]
        for shift in (-2,-1,1,2):
            xx,yy=(x+shift*step,y) if axis==0 else (x,y+shift*step)
            values.append(v.conservative_flux(xx,yy)[axis])
        derivatives.append((values[0]-8*values[1]+8*values[2]-values[3])/(12*step))
    return derivatives[0]+derivatives[1]


def generate():
    p=load_parameters()
    gas=reference_gas(p)
    gamma=gas.heat_capacity_ratio
    (STUDY/"data").mkdir(parents=True,exist_ok=True)
    (WRITING/"figures").mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    radius=np.linspace(0,p["profile_outer_radius_over_a"],p["profile_samples"])
    profiles=[]
    for m in p["profile_mach_values"]:
        v=IsentropicVortex(m,gamma)
        state=v.state(radius)
        audit=local_gas_audit(v,radius,p["radius_m"],gas,p["sutherland_temperature_K"])
        for i,r in enumerate(radius):
            profiles.append(dict(peak_mach_exterior=m,radius_over_a=r,radius_m=r*p["radius_m"],
                **{key:float(value[i]) for key,value in state.items()},
                **{key:float(value[i]) for key,value in audit.items()}))
    write_csv("radial_profiles.csv",profiles)
    mach=np.geomspace(*p["mach_scan_limits"],p["mach_scan_samples"])
    scan=[]
    for m in mach:
        v=IsentropicVortex(m,gamma)
        s=v.state(0)
        # Resolve maxima independently of the plotted radial grid.
        maximum_mach=minimize_scalar(lambda r:-float(v.state(r)["local_mach"]),
                                     bounds=(0,3),method="bounded")
        maximum_kn=minimize_scalar(lambda r:-float(local_gas_audit(v,r,p["radius_m"],gas,p["sutherland_temperature_K"])["kn_pressure"]),
                                   bounds=(0,4),method="bounded")
        maximum_rate=minimize_scalar(lambda r:-float(local_gas_audit(v,r,p["radius_m"],gas,p["sutherland_temperature_K"])["collision_shear_ratio"]),
                                     bounds=(0,4),method="bounded")
        scan.append(dict(peak_mach_exterior=m,center_density_deficit=float(s["density_deficit"]),
            center_density_ratio=float(s["density_ratio"]),center_temperature_ratio=float(s["temperature_ratio"]),
            center_pressure_ratio=float(s["pressure_ratio"]),
            center_incompressible_pressure_ratio=float(s["incompressible_pressure_ratio"]),
            center_pressure_difference_over_exterior=float(s["pressure_ratio"]-s["incompressible_pressure_ratio"]),
            max_local_mach=-maximum_mach.fun,max_thermodynamic_kn_proxy=-maximum_kn.fun,
            max_collision_shear_ratio=-maximum_rate.fun))
    write_csv("mach_scan.csv",scan)

    v=IsentropicVortex(p["solver"]["peak_mach"],gamma)
    refinements=[]
    outer=p["solver"]["outer_radius_over_a"]
    for n in p["solver"]["interval_counts"]:
        for boundary in ("exact","far_field"):
            r,rho=integrate_density(v,outer,n,boundary)
            exact=v.state(r)["density_ratio"]
            refinements.append(dict(intervals=n,step_over_a=outer/n,outer_radius_over_a=outer,
                boundary=boundary,max_density_error_over_exterior=float(np.max(np.abs(rho-exact))),
                center_density_error_over_exterior=float(rho[0]-exact[0])))
    write_csv("radial_refinement.csv",refinements)
    domains=[]
    for outer_value in p["solver"]["outer_radius_scan"]:
        n=int(np.ceil(outer_value/p["solver"]["domain_scan_step"]))
        r,rho=integrate_density(v,outer_value,n,"far_field")
        expected=(1-v.depression*(1-np.exp(-outer_value**2)))**(1/(gamma-1))
        true=float(v.state(0)["density_ratio"])
        boundary_error=true*np.expm1(np.log1p(v.depression*np.exp(-outer_value**2)/(1-v.depression))/(gamma-1))
        domains.append(dict(outer_radius_over_a=outer_value,intervals=n,radial_step_over_a=outer_value/n,
            center_error_vs_infinite_domain=float(rho[0]-true),
            exact_finite_boundary_error=float(boundary_error),
            radial_integration_error=float(rho[0]-expected)))
    write_csv("outer_boundary_scan.csv",domains)
    residuals=[]
    coords=np.linspace(-2.7,2.7,19)
    x,y=np.meshgrid(coords,coords)
    for step in p["residual_steps_over_a"]:
        res=conservative_residual(v,x,y,step)
        maxima=np.max(np.abs(res),axis=(0,1))
        residuals.append(dict(step_over_a=step,mass_residual=maxima[0],
            x_momentum_residual=maxima[1],y_momentum_residual=maxima[2],energy_residual=maxima[3]))
    write_csv("conservative_residual_scan.csv",residuals)

    fig,axes=plt.subplots(1,2,figsize=(11.4,4.3),constrained_layout=True)
    for m,color in zip(p["profile_mach_values"],COLORS):
        s=IsentropicVortex(m,gamma).state(radius)
        axes[0].plot(radius,s["density_ratio"],color=color,label=f"U/c_inf = {m:g}")
    axes[0].axhline(1,color=INK,ls=":",label="Constant-density reference")
    axes[0].set(xlim=(0,3.5),ylim=(.25,1.06),xlabel="Radius / peak-speed radius",
                ylabel="Density / exterior density",title="Faster rotation requires a deeper density depression")
    axes[0].legend(fontsize=8,loc="lower right")
    for m,color in zip(p["pressure_profile_mach_values"],(TEAL,GOLD)):
        s=IsentropicVortex(m,gamma).state(radius)
        axes[1].plot(radius,s["pressure_ratio"],color=color,label=f"Compressible, U/c_inf = {m:g}")
        axes[1].plot(radius,s["incompressible_pressure_ratio"],color=color,ls="--",label=f"Constant density, U/c_inf = {m:g}")
    axes[1].axhline(0,color=INK,lw=.8)
    axes[1].set(xlim=(0,3.5),xlabel="Radius / peak-speed radius",ylabel="Pressure / exterior pressure",
                title="The same exterior pressure does not give the same core")
    axes[1].legend(fontsize=7,loc="lower right")
    for ax in axes:
        ax.title.set_fontsize(9)
        ax.grid(alpha=.15)
    save_figure(fig,"01_density_and_pressure")

    fig,axes=plt.subplots(1,2,figsize=(11.4,4.2),constrained_layout=True)
    deficits=np.array([row["center_density_deficit"] for row in scan])
    marker=float(density_deficit_mach(p["density_deficit_marker"],gamma))
    axes[0].plot(mach,100*deficits,color=TEAL,label="Exact center density deficit")
    axes[0].plot(mach,100*.5*np.e*mach**2,color=INK,ls="--",label="Leading small-Mach term")
    axes[0].axhline(100*p["density_deficit_marker"],color=GOLD,ls=":")
    axes[0].axvline(marker,color=GOLD,ls=":",label=f"5% density change at U/c_inf = {marker:.3f}")
    axes[0].set(xlim=(0,p["mach_scan_limits"][1]),xlabel="Peak speed / exterior sound speed",
                ylabel="Center density deficit (%)",title="A stated density tolerance gives a profile-specific marker")
    axes[0].legend(fontsize=7,loc="upper left")
    difference=[row["center_pressure_difference_over_exterior"] for row in scan]
    axes[1].loglog(mach,difference,color=TEAL,label="Compressible minus constant-density pressure")
    axes[1].loglog(mach,gamma*np.e**2*mach**4/8,color=INK,ls="--",label="Leading fourth-power term")
    axes[1].set(xlabel="Peak speed / exterior sound speed",ylabel="Center pressure difference / exterior pressure",
                title="The pressure discrepancy is fourth order at small Mach")
    axes[1].legend(fontsize=7,loc="upper left")
    for ax in axes:
        ax.title.set_fontsize(9)
        ax.grid(alpha=.15)
    save_figure(fig,"02_departure_from_constant_density")

    fig,axes=plt.subplots(1,2,figsize=(11.4,4.2),constrained_layout=True)
    for boundary,color,label in (("exact",TEAL,"Exact outer density"),("far_field",GOLD,"Set outer density to exterior value")):
        selected=[r for r in refinements if r["boundary"]==boundary]
        axes[0].loglog([r["intervals"] for r in selected],
                      [r["max_density_error_over_exterior"] for r in selected],"o-",color=color,label=label)
    first=next(r for r in refinements if r["boundary"]=="exact")
    ns=np.array(p["solver"]["interval_counts"])
    axes[0].loglog(ns,first["max_density_error_over_exterior"]*(ns/ns[0])**-4.,
                  color=INK,ls=":",label="Fourth-order reference")
    axes[0].set(xlabel="Number of radial intervals",ylabel="Maximum density error / exterior density",
                title=f"Refining the grid cannot fix an outer datum at R = {outer:g}")
    axes[0].set_xticks(ns,[str(n) for n in ns])
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].legend(fontsize=7)
    axes[1].semilogy([r["outer_radius_over_a"] for r in domains],
                     [abs(r["center_error_vs_infinite_domain"]) for r in domains],"o-",color=GOLD,label="Total center error")
    axes[1].semilogy([r["outer_radius_over_a"] for r in domains],
                     [abs(r["radial_integration_error"]) for r in domains],"o-",color=TEAL,label="Radial integration error only")
    axes[1].semilogy([r["outer_radius_over_a"] for r in domains],
                     [r["exact_finite_boundary_error"] for r in domains],ls=":",color=INK,label="Exact outer-boundary contribution")
    axes[1].set(xlabel="Outer radius / peak-speed radius",ylabel="Center density error / exterior density",
                title="Separate outer-boundary and integration errors")
    axes[1].legend(fontsize=7)
    for ax in axes:
        ax.title.set_fontsize(9)
        ax.grid(alpha=.15)
    save_figure(fig,"03_radial_and_boundary_accuracy")

    fig,axes=plt.subplots(1,2,figsize=(11.4,4.2),constrained_layout=True)
    for m,color in zip(p["audit_mach_values"],COLORS):
        audit=local_gas_audit(IsentropicVortex(m,gamma),radius,p["radius_m"],gas,p["sutherland_temperature_K"])
        axes[0].semilogy(radius,np.maximum(audit["kn_pressure"],1e-12),color=color,label=f"U/c_inf = {m:g}")
        axes[1].semilogy(radius,np.maximum(audit["collision_shear_ratio"],1e-12),color=color,label=f"U/c_inf = {m:g}")
    axes[0].axhline(.01,color=INK,ls="--",label="Kn = 0.01 screening convention")
    axes[0].set(ylabel="Mean free path times |gradient of log pressure|",
                title="Thermodynamic variations remain well above gas transport scales",ylim=(1e-8,.03))
    axes[1].set(ylabel="Collision-time proxy times swirl shear rate",
                title="Collisions remain fast compared with the computed shear",ylim=(1e-8,.03))
    for ax in axes:
        ax.set(xlim=(0,4),xlabel="Radius / peak-speed radius")
        ax.title.set_fontsize(8)
        ax.legend(fontsize=7,loc="lower left")
        ax.grid(alpha=.15)
    save_figure(fig,"04_local_gas_scale_check")

    selected_errors=[r["max_density_error_over_exterior"] for r in refinements if r["boundary"]=="exact"]
    examples=observable_rows(p)
    summary=dict(five_percent_density_marker_mach=marker,
        density_marker_speed_m_s=float(marker*gas.sound_speed),
        mach_point_three_density_deficit=float(IsentropicVortex(.3,gamma).state(0)["density_deficit"]),
        constant_density_zero_center_pressure_mach=float(np.sqrt(2/(gamma*np.e))),
        positive_temperature_mach_limit=float(np.sqrt(2/((gamma-1)*np.e))),
        max_thermodynamic_kn_proxy=float(max(row["max_thermodynamic_kn_proxy"] for row in scan)),
        max_collision_shear_ratio=float(max(row["max_collision_shear_ratio"] for row in scan)),
        max_local_mach=float(max(row["max_local_mach"] for row in scan)),
        minimum_temperature_K=float(gas.temperature_K*min(row["center_temperature_ratio"] for row in scan)),
        finest_radial_density_error=float(selected_errors[-1]),
        radial_refinement_error_ratios=[float(x/y) for x,y in zip(selected_errors[:-1],selected_errors[1:])],
        finest_conservative_residual_max=float(max(list(residuals[-1].values())[1:])),
        representative_states=examples,
        scope="steady inviscid Euler equilibrium; no time evolution or blowup arrest",
        research_classification="known benchmark family; educational comparison; no novelty claim",
        parameters_sha256=hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest())
    (STUDY/"data/summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary


if __name__=="__main__":
    print(json.dumps(generate(),indent=2))
