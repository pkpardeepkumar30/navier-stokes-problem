"""D7 arithmetic, coefficient storage, and input-sensitivity experiments."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import mpmath as mp
import numpy as np
from .moments import SmoothBump
from .conditioning import (axial_system,stable_axial_solve,transformed_axial_system,
    high_precision_solve,coefficient_relative_error,reference_moment_remainder)

ROOT=Path(__file__).resolve().parents[3]
STUDY=ROOT/"code/d7_conditioning"
WRITING=ROOT/"writing/d7_conditioning"
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
    bumps=[SmoothBump(c,p["halfwidth"]) for c in p["centers"]]
    d=np.asarray(p["d_u"])
    digits=p["reference_digits"]
    (STUDY/"data").mkdir(parents=True,exist_ok=True)
    (WRITING/"figures").mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    scan=[]
    for lam in p["lambda_scan"]:
        reference_matrix,reference=high_precision_solve(lam,p["centers"],p["halfwidth"],d,digits)
        raw_matrix=axial_system(lam,bumps,p["quadrature_order"])
        transformed,_=transformed_axial_system(lam,bumps,d,p["quadrature_order"])
        stable=stable_axial_solve(lam,bumps,d,p["quadrature_order"])
        try:
            raw=np.linalg.solve(raw_matrix,-d)
            raw_error=coefficient_relative_error(raw,reference,digits)
            raw_remainder=reference_moment_remainder(reference_matrix,raw,d,digits)
            raw_status="solved"
        except np.linalg.LinAlgError:
            raw_error=raw_remainder=None
            raw_status="singular_float64_matrix"
        scan.append(dict(lambda_value=lam,raw_status=raw_status,
            raw_coefficient_relative_error=raw_error,
            stable_coefficient_relative_error=coefficient_relative_error(stable,reference,digits),
            raw_moment_relative_remainder=raw_remainder,
            stable_moment_relative_remainder=reference_moment_remainder(reference_matrix,stable,d,digits),
            reference_moment_relative_remainder=reference_moment_remainder(reference_matrix,reference,d,digits),
            raw_condition=float(np.linalg.cond(raw_matrix)),
            transformed_condition=float(np.linalg.cond(transformed)),
            stable_coefficient_1=stable[0],stable_coefficient_2=stable[1]))
    write_csv("lambda_accuracy.csv",scan)
    precision=[]
    lam=p["precision_scan"]["lambda"]
    matrix,reference=high_precision_solve(lam,p["centers"],p["halfwidth"],d,digits)
    for chosen in p["precision_scan"]["digits"]:
        _,values=high_precision_solve(lam,p["centers"],p["halfwidth"],d,chosen)
        precision.append(dict(decimal_digits=chosen,lambda_value=lam,
            coefficient_relative_error=coefficient_relative_error(values,reference,digits),
            moment_relative_remainder=reference_moment_remainder(matrix,values,d,digits),
            remainder_after_float64_storage=reference_moment_remainder(matrix,[float(v) for v in values],d,digits)))
    write_csv("precision_scan.csv",precision)
    sensitivity=[]
    s=p["sensitivity"]
    base=np.array(s["base_discrepancies"],float)
    base_coeff=stable_axial_solve(s["lambda"],bumps,base,p["quadrature_order"])
    for epsilon in s["perturbations"]:
        perturbed=base+np.array([0.,epsilon])
        coeff=stable_axial_solve(s["lambda"],bumps,perturbed,p["quadrature_order"])
        input_change=float(np.linalg.norm(perturbed-base)/np.linalg.norm(base))
        output_change=float(np.linalg.norm(coeff-base_coeff)/np.linalg.norm(base_coeff))
        sensitivity.append(dict(lambda_value=s["lambda"],nominal_perturbation=epsilon,
            actual_perturbation=float(perturbed[1]-base[1]),relative_input_change=input_change,
            relative_coefficient_change=output_change,measured_amplification=output_change/input_change))
    write_csv("input_sensitivity.csv",sensitivity)

    lambdas=np.array([row["lambda_value"] for row in scan])
    fig,axes=plt.subplots(1,2,figsize=(10.8,4),constrained_layout=True)
    for key,label,color in (("raw_coefficient_relative_error","Original float64 system",GOLD),
                            ("stable_coefficient_relative_error","Equivalent expm1 formulation",TEAL)):
        values=[np.nan if row[key] is None else max(row[key],1e-60) for row in scan]
        axes[0].loglog(lambdas,values,"o-",color=color,label=label)
    failed=[row["lambda_value"] for row in scan if row["raw_status"]!="solved"]
    if failed: axes[0].axvline(failed[0],color=GOLD,ls=":",label="Original matrix is singular here")
    axes[0].set(title="Equivalent equations can keep different digits",ylabel="Coefficient error vs 70-digit reference")
    axes[0].legend(fontsize=8)
    axes[1].loglog(lambdas,[row["raw_condition"] for row in scan],"o-",color=GOLD,label="Original moment rows")
    axes[1].loglog(lambdas,[row["transformed_condition"] for row in scan],"o-",color=TEAL,label="Transformed rows")
    axes[1].set(title="The row transformation separates the two weights",ylabel="2-norm condition number")
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.set_xlabel("lambda (decreasing)")
        ax.invert_xaxis()
        ax.set_xticks([1e-2,1e-6,1e-10,1e-14,1e-18])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.title.set_fontsize(10)
        ax.grid(alpha=.15)
    save_figure(fig,"01_equivalent_formulations")

    fig,axes=plt.subplots(1,2,figsize=(10.8,4),constrained_layout=True)
    for key,label,color in (("raw_moment_relative_remainder","Original solve, float64 coefficients",GOLD),
                            ("stable_moment_relative_remainder","Stable solve, float64 coefficients",TEAL),
                            ("reference_moment_relative_remainder","70-digit solve and retained coefficients",PURPLE)):
        values=[np.nan if row[key] is None else max(row[key],1e-60) for row in scan]
        axes[0].loglog(lambdas,values,"o-",color=color,label=label)
    axes[0].invert_xaxis()
    axes[0].set_xticks([1e-2,1e-6,1e-10,1e-14,1e-18])
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].set(title="Accurate large coefficients still need enough storage",xlabel="lambda (decreasing)",
                ylabel="Original moment remainder / supplied moment")
    axes[0].legend(fontsize=7)
    for key,label,color in (("moment_relative_remainder","Keep coefficients at working precision",TEAL),
                            ("remainder_after_float64_storage","Round coefficients to float64",GOLD)):
        axes[1].semilogy([row["decimal_digits"] for row in precision],
                         np.maximum([row[key] for row in precision],1e-60),"o-",color=color,label=label)
    axes[1].set(title="More precise arithmetic helps if its digits are retained",xlabel="Working decimal digits at lambda=1e-12",
                ylabel="Moment remainder against 70-digit integrals")
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.title.set_fontsize(9)
        ax.grid(alpha=.15)
    save_figure(fig,"02_precision_and_storage")

    fig,ax=plt.subplots(figsize=(7.4,4),constrained_layout=True)
    changes=np.array([row["relative_input_change"] for row in sensitivity])
    ax.loglog(changes,[row["relative_coefficient_change"] for row in sensitivity],"o-",color=TEAL,
              label="Computed response at lambda=1e-10")
    ax.loglog(changes,changes,ls="--",color=INK,label="Equal relative input and output changes")
    ax.set(title="An accurate solve can remain sensitive to its inputs",xlabel="Relative change in supplied discrepancy vector",
           ylabel="Relative change in correction coefficients")
    ax.legend(fontsize=8)
    ax.grid(alpha=.15)
    save_figure(fig,"03_input_sensitivity")

    selected=next(row for row in scan if row["lambda_value"]==p["precision_scan"]["lambda"])
    summary=dict(reference_digits=digits,
        largest_stable_coefficient_relative_error=max(row["stable_coefficient_relative_error"] for row in scan),
        raw_singular_lambda_values=failed,
        selected_lambda=selected["lambda_value"],
        selected_raw_coefficient_error=selected["raw_coefficient_relative_error"],
        selected_stable_coefficient_error=selected["stable_coefficient_relative_error"],
        selected_stable_original_moment_remainder=selected["stable_moment_relative_remainder"],
        sensitivity_amplification_range=[min(row["measured_amplification"] for row in sensitivity),
                                       max(row["measured_amplification"] for row in sensitivity)],
        classification="educational numerical-analysis note; no novel research claim",
        dynamical_stability=False,
        parameters_sha256=hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest())
    (STUDY/"data/summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary


if __name__=="__main__":
    print(json.dumps(generate(),indent=2))
