"""D4 reproducible momentum-flux figures and numerical diagnostics."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import numpy as np
from .stress import (periodic_grid, derivative, two_mode_velocity, two_mode_target,
    covariance, divergence, envelope, enveloped_velocity, envelope_covariance,
    envelope_mean_force, averaged_transport, force_from_mean_covariance, project_mean_force)

ROOT = Path(__file__).resolve().parents[3]
STUDY = ROOT/"code/d4_mean_stress"
WRITING = ROOT/"writing/d4_mean_stress"
TEAL, GOLD, INK, PURPLE = "#147d83", "#bc681b", "#1d2b3a", "#7558a4"


def load_parameters():
    return json.loads((STUDY/"parameters.json").read_text(encoding="utf-8"))


def write_csv(name, rows):
    with (STUDY/"data"/name).open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig,name):
    for ext in ("png","svg"):
        options={"metadata":{"Date":None}} if ext=="svg" else {}
        fig.savefig(WRITING/"figures"/f"{name}.{ext}",dpi=180,bbox_inches="tight",**options)
    plt.close(fig)


def diagnostics(p=None):
    p=p or load_parameters()
    points,x,y=periodic_grid(p["domain"]["grid_size"])
    constant=two_mode_velocity(x,y,**p["constant"])
    target=two_mode_target(p["constant"]["lambda_plus"],p["constant"]["lambda_minus"])
    q=covariance(constant,(0,1))
    w=enveloped_velocity(x,y,**p["envelope"])
    exact,leading=envelope_covariance(points,**p["envelope"])
    meanq=covariance(w,0)
    analytic_force=envelope_mean_force(points,p["envelope"]["epsilon"],p["envelope"]["slope"])
    force=force_from_mean_covariance(meanq)
    projected=project_mean_force(force)
    summary=dict(grid_size=len(points),
        constant_mean_velocity_max=float(np.max(abs(constant.mean(axis=(0,1))))),
        constant_stress_relative_error=float(np.linalg.norm(q-target)/np.linalg.norm(target)),
        constant_divergence_max=float(np.max(abs(divergence(constant)))),
        envelope_mean_velocity_max=float(np.max(abs(w.mean(axis=0)))) ,
        envelope_stress_max_error=float(np.max(abs(meanq-exact))),
        envelope_divergence_max=float(np.max(abs(divergence(w)))),
        envelope_force_max_error=float(np.max(abs(force-analytic_force))),
        transport_force_max_error=float(np.max(abs(-averaged_transport(w)-analytic_force))),
        envelope_covariance_min_eigenvalue=float(np.linalg.eigvalsh(meanq).min()),
        projected_force_max=float(np.max(abs(projected[:,0]))),
        net_mean_force_max=float(np.max(abs(force.mean(axis=0)))),
        analytic_finite_frequency_correction_max=p["envelope"]["epsilon"]**2/(2*p["envelope"]["k"]**2),
        paper_stress_reconstructed=False, time_integration=False,
        parameters_sha256=hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest())
    return points,x,y,constant,q,w,meanq,exact,leading,force,projected,summary


def diagnostic_table(summary):
    items=[("Constant target: relative covariance error", "constant_stress_relative_error"),
           ("Envelope: maximum mean-velocity component", "envelope_mean_velocity_max"),
           ("Envelope: maximum absolute divergence", "envelope_divergence_max"),
           ("Envelope: maximum covariance error", "envelope_stress_max_error"),
           ("Mean force: error from differentiated covariance", "envelope_force_max_error"),
           ("Mean force: error from direct transport", "transport_force_max_error")]
    return "| Check | Computed value |\n| --- | --- |\n"+"\n".join(
        f"| {label} | {summary[key]:.3g} |" for label,key in items)


def generate():
    p=load_parameters()
    (STUDY/"data").mkdir(parents=True,exist_ok=True)
    (WRITING/"figures").mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    points,x,y,constant,q,w,meanq,exact,leading,force,projected,summary=diagnostics(p)
    target=two_mode_target(p["constant"]["lambda_plus"],p["constant"]["lambda_minus"])
    write_csv("sine_samples.csv",[dict(phase=s,sine=np.sin(s),sine_squared=np.sin(s)**2) for s in points])
    write_csv("constant_covariance.csv",[dict(row=i,column=j,target=target[i,j],computed=q[i,j],
              absolute_error=abs(q[i,j]-target[i,j])) for i in range(2) for j in range(2)])
    a,ap,app=envelope(points,p["envelope"]["epsilon"])
    velocity_mean=w.mean(axis=0)
    direct_force=-averaged_transport(w)
    rows=[]
    for i,s in enumerate(points):
        rows.append(dict(y=s,envelope=a[i],mean_wx=velocity_mean[i,0],mean_wy=velocity_mean[i,1],
            Qxx=meanq[i,0,0],Qxy=meanq[i,0,1],Qyy=meanq[i,1,1],
            exact_Qxx=exact[i,0,0],exact_Qxy=exact[i,0,1],exact_Qyy=exact[i,1,1],
            pressure=-a[i]**2/2,force_x=force[i,0],force_y=force[i,1],
            direct_transport_force_x=direct_force[i,0],direct_transport_force_y=direct_force[i,1],
            projected_force_x=projected[i,0],projected_force_y=projected[i,1],
            projected_curl=-p["envelope"]["slope"]*(ap[i]**2+a[i]*app[i])))
    write_csv("spatial_stress_force.csv",rows)
    # Exact plotting-grid values for the velocity image; coordinates are separate columns.
    write_csv("envelope_velocity_grid.csv",[dict(x=x[i,j],y=y[i,j],wx=w[i,j,0],wy=w[i,j,1])
        for i in range(len(points)) for j in range(len(points))])
    resolution=[]
    maxmode=max(p["constant"]["n"],p["constant"]["m"])
    for size in p["resolution_scan"]:
        _,xx,yy=periodic_grid(size)
        wc=two_mode_velocity(xx,yy,**p["constant"])
        qc=covariance(wc,(0,1))
        resolution.append(dict(grid_size=size,
            stress_relative_error=float(np.linalg.norm(qc-target)/np.linalg.norm(target)),
            mean_velocity_max=float(np.max(abs(wc.mean(axis=(0,1))))),
            divergence_max=float(np.max(abs(divergence(wc)))),
            velocity_modes_below_nyquist=size>2*maxmode,
            all_quadratic_modes_below_nyquist=size>4*maxmode))
    write_csv("resolution_scan.csv",resolution)
    frequency=[]
    grid_size=p["frequency_scan_grid"]
    pts,xx,yy=periodic_grid(grid_size)
    for k in p["frequency_scan"]:
        e={**p["envelope"],"k":k}
        wv=enveloped_velocity(xx,yy,**e)
        analytic,leading_k=envelope_covariance(pts,**e)
        qk=covariance(wv,0)
        error=float(np.max(abs(qk-leading_k)))
        exact_error=p["envelope"]["epsilon"]**2/(2*k*k)
        frequency.append(dict(k=k,grid_size=grid_size,
            shortest_velocity_axis_wavelength=2*np.pi/max(k,abs(e["slope"])*k+1),
            all_quadratic_modes_below_nyquist=grid_size>4*max(k,abs(e["slope"])*k+1),
            finite_frequency_covariance_correction=error,exact_correction=exact_error,
            covariance_max_error=float(np.max(abs(qk-analytic))),
            divergence_max=float(np.max(abs(divergence(wv)))),
            projected_force_max=float(np.max(abs(project_mean_force(force_from_mean_covariance(qk))[:,0])))))
    write_csv("frequency_scan.csv",frequency)

    fig,axes=plt.subplots(1,3,figsize=(11.8,3.6),constrained_layout=True)
    phase=np.r_[points,2*np.pi]
    for ax,values,title,mean,color in (
        (axes[0],np.sin(phase),"A velocity component can average to zero",0,TEAL),
        (axes[1],np.sin(phase)**2,"Its square keeps a positive average",0.5,GOLD)):
        ax.plot(phase,values,color=color,lw=2)
        ax.axhline(mean,color=INK,ls="--",label=f"Mean = {mean:g}")
        ax.set(xlabel="Phase (radians)",ylabel="Dimensionless value",title=title)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.15)
        ax.set_xticks([0,np.pi,2*np.pi],["0","pi","2 pi"])
        ax.title.set_fontsize(10)
    heat=axes[2].imshow(q,cmap="YlGnBu",vmin=0,vmax=3)
    axes[2].set(title="Two waves realize the target covariance",
                xticks=[0,1],xticklabels=["x","y"],yticks=[0,1],yticklabels=["x","y"])
    axes[2].title.set_fontsize(10)
    for i in range(2):
        for j in range(2): axes[2].text(j,i,f"{q[i,j]:.2f}",ha="center",va="center",fontsize=18,
                                      color="white" if q[i,j]>1.5 else INK)
    fig.colorbar(heat,ax=axes[2],label="Mean velocity product")
    save_figure(fig,"01_zero_mean_nonzero_stress")

    fig,axes=plt.subplots(2,2,figsize=(10.8,8.1),constrained_layout=True)
    bound=float(np.max(abs(w[...,0])))
    im=axes[0,0].pcolormesh(points,points,w[...,0].T,cmap="RdBu_r",vmin=-bound,vmax=bound,
                           shading="nearest",rasterized=True)
    axes[0,0].set(title="Velocity oscillates; its envelope varies with y",xlabel="x",ylabel="y",aspect="equal")
    fig.colorbar(im,ax=axes[0,0],label="Horizontal fluctuation wx")
    axes[0,1].plot(points,meanq[:,0,1],color=TEAL,lw=2,label="Mean wx wy")
    axes[0,1].plot(points[::8],exact[::8,0,1],"o",mfc="none",color=INK,label="Exact covariance")
    axes[0,1].set(title="Averaging over x leaves a varying flux",xlabel="y",ylabel="Qxy(y)")
    axes[0,1].legend(fontsize=9)
    axes[1,0].plot(points,force[:,0],color=TEAL,label="Horizontal wave force")
    axes[1,0].plot(points,force[:,1],color=GOLD,label="Vertical wave force")
    axes[1,0].set(title="The flux gradient acts on the mean flow",xlabel="y",ylabel="Components of -div Q")
    axes[1,0].legend(fontsize=9)
    axes[1,1].plot(points,projected[:,0],color=TEAL,label="Horizontal force remains")
    axes[1,1].plot(points,projected[:,1],color=GOLD,ls="--",label="Vertical part cancels with pressure")
    axes[1,1].set(title="Pressure removes the gradient contribution",xlabel="y",ylabel="Force after pressure balance")
    axes[1,1].legend(fontsize=9)
    for ax in axes.flat:
        if ax!=axes[0,0]: ax.grid(alpha=0.15)
        ax.title.set_fontsize(10)
    save_figure(fig,"02_spatial_transport")

    fig,axes=plt.subplots(1,2,figsize=(11,4.1),constrained_layout=True)
    axes[0].semilogy([r["grid_size"] for r in resolution],
                     np.maximum([r["stress_relative_error"] for r in resolution],1e-17),"o-",color=TEAL)
    axes[0].axvspan(0,2*maxmode,color=GOLD,alpha=0.12,label="Carrier not resolved")
    axes[0].axvline(4*maxmode,color=PURPLE,ls="--",label="All quadratic modes require N > 36")
    axes[0].set(xlim=(12,132),title="One accurate average can hide missing modes",
                xlabel="Grid points in each direction, N",ylabel="Relative covariance error")
    axes[0].legend(fontsize=8)
    kvals=np.array([r["k"] for r in frequency])
    axes[1].loglog(kvals,[r["exact_correction"] for r in frequency],color=INK,label="Exact: epsilon^2 / (2 k^2)")
    axes[1].loglog(kvals,[r["finite_frequency_covariance_correction"] for r in frequency],
                   "o",color=TEAL,mfc="none",label="Computed covariance correction")
    axes[1].set(title="The envelope correction shrinks as k^(-2)",
                xlabel="Integer carrier frequency, k",ylabel="Maximum entry of Q - leading covariance")
    axes[1].legend(fontsize=8)
    axes[1].set_xticks(kvals, [str(k) for k in kvals])
    axes[1].xaxis.set_minor_locator(NullLocator())
    for ax in axes: ax.grid(alpha=0.15)
    save_figure(fig,"03_resolution_and_envelope")
    summary["coarsest_grid_stress_relative_error"]=resolution[0]["stress_relative_error"]
    summary["frequency_scan_correction_ratios"]=[frequency[i]["finite_frequency_covariance_correction"]/
        frequency[i+1]["finite_frequency_covariance_correction"] for i in range(len(frequency)-1)]
    (STUDY/"data/summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary


if __name__=="__main__":
    print(json.dumps(generate(),indent=2))
