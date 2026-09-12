"""D5 time histories, snapshots, and independent amplitude refinement."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import minimize_scalar
from .shear import KelvinWave, integrate_amplitude

ROOT=Path(__file__).resolve().parents[3]
STUDY=ROOT/"code/d5_shearing_wave"
WRITING=ROOT/"writing/d5_shearing_wave"
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
        options={"metadata":{"Date":None}} if ext=="svg" else {}
        fig.savefig(WRITING/"figures"/f"{name}.{ext}",dpi=180,bbox_inches="tight",**options)
    plt.close(fig)


def generate():
    p=load_parameters()
    wave=KelvinWave(**p["wave"])
    (STUDY/"data").mkdir(parents=True,exist_ok=True)
    (WRITING/"figures").mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    times=np.linspace(0,p["time"]["end"],p["time"]["samples"])
    state=wave.state(times)
    inviscid=KelvinWave(**{**p["wave"],"viscosity":0}).state(times)
    peak=minimize_scalar(lambda t:-float(wave.state(t)["speed"]),
        bounds=(0,p["time"]["end"]),method="bounded",options={"xatol":1e-12})
    if not peak.success: raise RuntimeError("speed-peak search did not converge")
    peak_time=float(peak.x)
    peak_speed=float(-peak.fun)
    # Check endpoints too; the restricted example has a single interior peak,
    # while a changed parameter set can have its largest value at an endpoint.
    candidates=[(0,float(state["speed"][0])),(peak_time,peak_speed),
                (float(times[-1]),float(state["speed"][-1]))]
    peak_time,peak_speed=max(candidates,key=lambda pair:pair[1])
    rows=[]
    for i,t in enumerate(times):
        row=dict(time=t,kx=wave.kx,ky=state["ky"][i],
                 Bx=state["velocity_amplitude"][i,0],By=state["velocity_amplitude"][i,1])
        row.update({key:values[i] for key,values in state.items() if key not in ("ky","velocity_amplitude")})
        row["inviscid_speed"]=inviscid["speed"][i]
        rows.append(row)
    write_csv("time_history.csv",rows)
    refinement=[]
    for step in p["refinement_steps"]:
        tt,bb=integrate_amplitude(wave,p["time"]["end"],step)
        exact=wave.state(tt)
        errors=np.linalg.norm(bb-exact["velocity_amplitude"],axis=1)
        qxy=bb[:,0]*bb[:,1]/2
        refinement.append(dict(step=step,steps=len(tt)-1,
            max_velocity_error=float(np.max(errors)),
            relative_to_initial_speed=float(np.max(errors)/wave.initial_speed),
            max_covariance_error=float(np.max(abs(qxy-exact["qxy"]))),
            max_constraint_defect=float(np.max(abs(wave.kx*bb[:,0]+exact["ky"]*bb[:,1]))),
            sampled_peak_speed=float(np.max(np.linalg.norm(bb,axis=1))),
            energy_budget_error=float(abs(simpson(exact["energy_rate"],x=tt)-
                    (exact["energy"][-1]-exact["energy"][0])))))
    write_csv("time_refinement.csv",refinement)

    snap=p["snapshots"]
    points=np.linspace(0,snap["domain_length"],snap["grid_size"],endpoint=False)
    xx,yy=np.meshgrid(points,points,indexing="ij")
    zeta0=wave.initial_speed*np.hypot(wave.kx,wave.ky0)
    snapshot_rows=[]
    snapshot_fields=[]
    for t in snap["times"]:
        s=wave.state(t)
        vorticity=s["vorticity_amplitude"]/zeta0*np.cos(wave.kx*xx+s["ky"]*yy)
        snapshot_fields.append(vorticity)
        snapshot_rows.append(dict(time=t,kx=wave.kx,ky=float(s["ky"]),speed=float(s["speed"]),
            vorticity_ratio=float(s["vorticity_amplitude"]/zeta0),wavelength=float(s["wavelength"]),
            max_coordinate_wavenumber=max(abs(wave.kx),abs(float(s["ky"]))),
            grid_size=snap["grid_size"],spacing=snap["domain_length"]/snap["grid_size"],
            minimum_axis_samples_per_wavelength=2*np.pi/(max(abs(wave.kx),abs(float(s["ky"])))*
                                                       snap["domain_length"]/snap["grid_size"])))
    write_csv("snapshot_parameters.csv",snapshot_rows)
    # These exact phase samples provide compact reproducible profiles for each snapshot.
    phase=2*np.pi*np.arange(256)/256
    write_csv("snapshot_phase_profiles.csv",[dict(time=t,phase=ph,
              vorticity_ratio=wave.state(t)["vorticity_amplitude"]/zeta0*np.cos(ph),
              wx=wave.state(t)["velocity_amplitude"][0]*np.sin(ph),
              wy=wave.state(t)["velocity_amplitude"][1]*np.sin(ph))
              for t in snap["times"] for ph in phase])

    fig,axes=plt.subplots(2,2,figsize=(11,7.8),constrained_layout=True)
    axes[0,0].plot(times,state["speed"]/wave.initial_speed,color=TEAL,label="Shear + viscosity")
    axes[0,0].plot(times,inviscid["speed"]/wave.initial_speed,color=PURPLE,ls="--",label="Shear, zero viscosity")
    axes[0,0].plot(times,state["heat_reference_speed"]/wave.initial_speed,color=GOLD,label="Heat reference: no shear")
    axes[0,0].scatter([peak_time],[peak_speed/wave.initial_speed],color=TEAL,s=25)
    axes[0,0].set(title="The wave grows before it decays",ylabel="Peak wave speed / initial speed")
    axes[0,0].legend(fontsize=8)
    axes[0,1].plot(times,state["wavelength"],color=TEAL)
    axes[0,1].set(title="Wavelength first lengthens, then shortens",ylabel="Normal wavelength / reference length")
    axes[1,0].plot(times,state["production"],color=TEAL,label="From background shear: -S Qxy")
    axes[1,0].plot(times,-state["dissipation"],color=GOLD,label="Viscous loss")
    axes[1,0].plot(times,state["energy_rate"],color=INK,ls="--",label="Net energy change")
    axes[1,0].axhline(0,color="#8c98a1",lw=0.8)
    axes[1,0].set(title="The energy budget identifies the source",ylabel="Energy transfer per unit time")
    axes[1,0].legend(fontsize=8)
    axes[1,1].plot(times,state["qxy"],color=TEAL,label="Qxy = mean(wx wy)")
    axes[1,1].axhline(0,color="#8c98a1",lw=0.8)
    axes[1,1].set(title="The mean momentum flux changes sign",ylabel="Momentum covariance")
    axes[1,1].legend(fontsize=8)
    for ax in axes.flat:
        ax.set_xlabel("Dimensionless time")
        ax.grid(alpha=0.15)
        ax.title.set_fontsize(10)
    save_figure(fig,"01_wave_lifecycle")

    fig,axes=plt.subplots(1,len(snap["times"]),figsize=(11.7,4),constrained_layout=True,squeeze=False)
    for ax,t,values in zip(axes[0],snap["times"],snapshot_fields):
        s=wave.state(t)
        im=ax.pcolormesh(points,points,values.T,cmap="RdBu_r",vmin=-1,vmax=1,
                         shading="nearest",rasterized=True)
        ax.set(aspect="equal",xlabel="x / reference length",ylabel="y / reference length",
               title=f"t={t:g}; ky={float(s['ky']):g}\nSpeed / initial = {float(s['speed']/wave.initial_speed):.2f}")
        ax.title.set_fontsize(10)
    fig.colorbar(im,ax=list(axes[0]),label="Perturbation vorticity / initial amplitude",shrink=0.8)
    save_figure(fig,"02_wave_orientation")

    fig,axes=plt.subplots(1,2,figsize=(10.8,4),constrained_layout=True)
    steps=np.array([r["step"] for r in refinement])
    errors=np.array([r["max_velocity_error"] for r in refinement])
    axes[0].loglog(steps,errors,"o-",color=TEAL,label="RK4 vs analytic velocity amplitude")
    axes[0].loglog(steps,errors[0]*(steps/steps[0])**4,ls="--",color=INK,label="Fourth-order reference")
    axes[0].set(title="An independent time integrator converges",xlabel="Time step (decreasing)",
                ylabel="Maximum vector amplitude error")
    axes[0].invert_xaxis()
    axes[0].set_xticks(steps,[f"{step:g}" for step in steps])
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].legend(fontsize=8)
    axes[1].semilogy(times,np.exp(-wave.viscosity*state["damping_integral"]),color=TEAL,label="Sheared wave: exp(-nu integral |k|^2 dt)")
    axes[1].semilogy(times,state["heat_reference_speed"]/wave.initial_speed,color=GOLD,label="Fixed wavevector: exp(-nu |k0|^2 t)")
    axes[1].set(title="Changing wavelength changes accumulated damping",xlabel="Dimensionless time",
                ylabel="Vorticity amplitude / initial amplitude")
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.title.set_fontsize(10)
        ax.grid(alpha=0.15)
    save_figure(fig,"03_time_accuracy_and_damping")

    summary=dict(peak_time=peak_time,peak_speed_ratio=peak_speed/wave.initial_speed,
        peak_energy_ratio=(peak_speed/wave.initial_speed)**2,
        turn_time=wave.ky0/(wave.shear*wave.kx),
        shortest_wavelength=float(np.min(state["wavelength"])),
        longest_wavelength=float(np.max(state["wavelength"])),
        final_speed_ratio=float(state["speed"][-1]/wave.initial_speed),
        final_vorticity_ratio=float(state["vorticity_amplitude"][-1]/zeta0),
        time_sample_spacing=float(times[1]-times[0]),
        finest_rk4_amplitude_error=refinement[-1]["max_velocity_error"],
        finest_rk4_constraint_defect=refinement[-1]["max_constraint_defect"],
        refinement_error_ratios=list(errors[:-1]/errors[1:]),
        integrated_energy_budget_error=float(abs(simpson(state["energy_rate"],x=times)-
                                                    (state["energy"][-1]-state["energy"][0]))),
        snapshot_minimum_axis_samples_per_wavelength=min(r["minimum_axis_samples_per_wavelength"] for r in snapshot_rows),
        paper_pulse_reconstructed=False,
        parameters_sha256=hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest())
    (STUDY/"data/summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary


if __name__=="__main__":
    print(json.dumps(generate(),indent=2))
