"""D6 stage A: source-defined moment repair and its conditioning limits."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import numpy as np
from .moments import (SmoothBump, moment_matrices, solve_moments, combine_bumps,
                      original_moments, repaired_moments)

ROOT=Path(__file__).resolve().parents[3]
STUDY=ROOT/"code/d6_moment_correction"
WRITING=ROOT/"writing/d6_moment_correction"
TEAL,GOLD,INK,PURPLE="#147d83","#bc681b","#1d2b3a","#7558a4"


def load_parameters():
    return json.loads((STUDY/"parameters.json").read_text(encoding="utf-8"))


def setup(p=None):
    p=p or load_parameters()
    return ([SmoothBump(c,p["bump_halfwidth"]) for c in p["u_bump_centers"]],
            [SmoothBump(c,p["bump_halfwidth"]) for c in p["e_bump_centers"]])


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
    u,e=setup(p)
    (STUDY/"data").mkdir(parents=True,exist_ok=True)
    (WRITING/"figures").mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    ef=p["patch_amplitude"]/(1+p["eta"]**2)
    before=original_moments(p["d_u"],p["d_e"],ef)
    alpha,beta=solve_moments(p["lambda"],u,e,p["d_u"],p["d_e"],p["quadrature_order"])
    BU,BE=moment_matrices(p["lambda"],u,e,p["reference_order"])
    after=repaired_moments(BU,BE,alpha,beta,before,ef)
    scale=np.maximum(abs(before),1e-30)
    write_csv("moment_balance.csv",[dict(moment=i+1,before=before[i],after=after[i],
                    absolute_remainder=abs(after[i]),relative_remainder=abs(after[i])/scale[i]) for i in range(5)])
    write_csv("moment_matrices.csv",[dict(block=name,row=i,column=j,value=float(B[i,j]))
        for name,B in (("U",BU),("E",BE)) for i in range(len(B)) for j in range(B.shape[1])])
    write_csv("coefficients.csv",[dict(block=name,index=i+1,coefficient=value,center=b.center,halfwidth=b.halfwidth)
        for name,coef,bumps in (("U",alpha,u),("E",beta,e)) for i,(value,b) in enumerate(zip(coef,bumps))])
    r=np.linspace(min(b.center-b.halfwidth for b in u+e)-.2,max(b.center+b.halfwidth for b in u+e)+.2,2401)
    du,de=combine_bumps(r,alpha,u),combine_bumps(r,beta,e)
    write_csv("correction_profiles.csv",[dict(radius=rr,delta_U=uu,delta_E=ee) for rr,uu,ee in zip(r,du,de)])
    convergence=[]
    for order in p["quadrature_orders"]:
        a,b=solve_moments(p["lambda"],u,e,p["d_u"],p["d_e"],order)
        rem=repaired_moments(BU,BE,a,b,before,ef)
        convergence.append(dict(order=order,reference_order=p["reference_order"],
            maximum_absolute_moment_remainder=float(np.max(abs(rem))),
            maximum_relative_moment_remainder=float(np.max(abs(rem)/scale)),
            max_coefficient_difference=float(np.max(abs(np.r_[a-alpha,b-beta])))))
    write_csv("quadrature_convergence.csv",convergence)
    scan=[]
    for lam in p["lambda_scan"]:
        M,N=moment_matrices(lam,u,e,p["quadrature_order"])
        Mr,Nr=moment_matrices(lam,u,e,p["reference_order"])
        a,b=solve_moments(lam,u,e,p["d_u"],p["d_e"],p["quadrature_order"])
        rem=repaired_moments(Mr,Nr,a,b,before,ef)
        scan.append(dict(lambda_value=lam,condition_U=float(np.linalg.cond(M)),
            row_equilibrated_condition_U=float(np.linalg.cond(M/np.linalg.norm(M,axis=1)[:,None])),
            condition_E=float(np.linalg.cond(N)),max_axial_coefficient=float(np.max(abs(a))),
            max_angular_coefficient=float(np.max(abs(b))),maximum_relative_moment_remainder=float(np.max(abs(rem)/scale))))
    write_csv("lambda_conditioning.csv",scan)

    fig,axes=plt.subplots(1,2,figsize=(10.8,4),constrained_layout=True)
    for ax,values,color,title in ((axes[0],du,TEAL,"Two signed bumps repair the axial moments"),
                                 (axes[1],de,GOLD,"Three signed bumps repair the angular moments")):
        ax.plot(r,values,color=color)
        ax.axhline(0,color="#8b959d",lw=.8)
        ax.set(xlabel="Profile radius R = sqrt(2X)",ylabel="Added profile coefficient",title=title)
        ax.title.set_fontsize(10)
        ax.grid(alpha=.15)
    axes[0].set_xlim(1.6,3.4)
    axes[1].set_xlim(3.6,6.4)
    save_figure(fig,"01_local_profile_repair")

    fig,axes=plt.subplots(1,2,figsize=(10.8,4),constrained_layout=True)
    ids=np.arange(1,6)
    axes[0].semilogy(ids,np.ones(5),"o-",color=GOLD,label="Before correction")
    axes[0].semilogy(ids,np.maximum(abs(after)/scale,1e-17),"o",color=TEAL,label="After correction; higher-order quadrature")
    axes[0].set(xticks=ids,xlabel="Moment number",ylabel="Absolute moment / its initial magnitude",
                title="The five supplied discrepancies are cancelled")
    axes[0].legend(fontsize=8,loc="center right")
    orders=np.array([row["order"] for row in convergence])
    axes[1].semilogy(orders,np.maximum([row["maximum_relative_moment_remainder"] for row in convergence],1e-17),"o-",color=TEAL)
    axes[1].set(xlabel="Quadrature order used to compute correction",ylabel="Maximum relative moment remainder",
                title="Reintegration distinguishes repair from solve error")
    axes[1].set_xticks(orders)
    for ax in axes:
        ax.title.set_fontsize(10)
        ax.grid(alpha=.15)
    save_figure(fig,"02_moment_cancellation")

    fig,axes=plt.subplots(1,3,figsize=(12,3.8),constrained_layout=True)
    lambdas=np.array([row["lambda_value"] for row in scan])
    axes[0].loglog(lambdas,[row["condition_U"] for row in scan],"o-",color=TEAL,label="Axial block")
    axes[0].loglog(lambdas,[row["condition_E"] for row in scan],"s-",color=GOLD,label="Angular block")
    axes[0].set(title="Distinct powers become nearly equal",ylabel="Matrix condition number")
    axes[0].legend(fontsize=8)
    axes[1].loglog(lambdas,[row["max_axial_coefficient"] for row in scan],"o-",color=TEAL)
    axes[1].set(title="Repairing fixed inputs needs larger bumps",ylabel="Maximum axial coefficient")
    axes[2].loglog(lambdas,np.maximum([row["maximum_relative_moment_remainder"] for row in scan],1e-17),"o-",color=PURPLE)
    axes[2].set(title="Finite precision limits the cancellation",ylabel="Maximum relative moment remainder")
    for ax in axes:
        ax.set_xlabel("lambda (decreasing)")
        ax.invert_xaxis()
        ax.set_xticks([.1,.001,1e-5,1e-7])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.title.set_fontsize(9)
        ax.grid(alpha=.15)
    save_figure(fig,"03_conditioning_limit")
    summary=dict(maximum_absolute_moment_remainder=float(np.max(abs(after))),
        maximum_relative_moment_remainder=float(np.max(abs(after)/scale)),
        baseline_condition_U=float(np.linalg.cond(BU)),baseline_condition_E=float(np.linalg.cond(BE)),
        baseline_max_axial_coefficient=float(np.max(abs(alpha))),
        smallest_lambda=scan[-1]["lambda_value"],smallest_lambda_condition_U=scan[-1]["condition_U"],
        smallest_lambda_relative_moment_remainder=scan[-1]["maximum_relative_moment_remainder"],
        pde_residual_evaluated=False,full_order_one_correction=False,
        parameters_sha256=hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest())
    (STUDY/"data/summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary


if __name__=="__main__":
    print(json.dumps(generate(),indent=2))
