"""Stable evaluation and independent precision reference for the axial moment block."""
from functools import lru_cache
import mpmath as mp
import numpy as np
from .moments import quadrature


def axial_system(lam,bumps,order=128):
    if not np.isfinite(lam) or lam<=0 or len(bumps)!=2:
        raise ValueError("positive lambda and two bumps required")
    return np.array([[b.moment(p,order) for b in bumps] for p in (1,1-2*lam)])


def transformed_axial_system(lam,bumps,discrepancies,order=128):
    if not np.isfinite(lam) or lam<=0 or len(bumps)!=2:
        raise ValueError("positive lambda and two bumps required")
    d=np.asarray(discrepancies,float)
    if d.shape!=(2,) or not np.all(np.isfinite(d)):
        raise ValueError("two finite discrepancies required")
    nodes,weights=quadrature(order)
    first=np.array([b.center for b in bumps])  # Exact by symmetry and unit integral.
    second=[]
    for b in bumps:
        r=b.center+b.halfwidth*nodes
        second.append(np.dot(weights,r*np.expm1(-2*lam*np.log(r))/(-2*lam)))
    matrix=np.stack((first,second))
    rhs=np.array([-d[0],(d[1]-d[0])/(2*lam)])
    return matrix,rhs


def stable_axial_solve(lam,bumps,discrepancies,order=128):
    matrix,rhs=transformed_axial_system(lam,bumps,discrepancies,order)
    return np.linalg.solve(matrix,rhs)


@lru_cache(maxsize=64)
def reference_system(lam_text,centers_text,width_text,digits=70):
    """Independent adaptive integral of the original powers, using decimal inputs."""
    with mp.workdps(digits):
        lam,width=mp.mpf(lam_text),mp.mpf(width_text)
        def kernel(z):
            return mp.exp(-1/(1-z*z)) if abs(z)<1 else mp.mpf(0)
        normal=mp.quad(kernel,[-1,0,1])
        centers=[mp.mpf(c) for c in centers_text]
        matrix=mp.matrix(2)
        for j,c in enumerate(centers):
            matrix[0,j]=c
            matrix[1,j]=mp.quad(lambda z:(c+width*z)**(1-2*lam)*kernel(z),[-1,0,1])/normal
        return matrix


def high_precision_solve(lam,centers,width,discrepancies,digits=70):
    with mp.workdps(digits):
        matrix=reference_system(str(lam),tuple(str(c) for c in centers),str(width),digits)
        rhs=mp.matrix([-mp.mpf(str(d)) for d in discrepancies])
        return matrix,mp.lu_solve(matrix,rhs)


def coefficient_relative_error(values,reference,digits=70):
    with mp.workdps(digits):
        # mp.mpf(float) retains the actual stored binary coefficient.
        vector=mp.matrix([mp.mpf(v) for v in values])
        return float(mp.norm(vector-reference)/mp.norm(reference))


def reference_moment_remainder(matrix,values,discrepancies,digits=70):
    with mp.workdps(digits):
        vector=mp.matrix([mp.mpf(v) for v in values])
        d=mp.matrix([mp.mpf(str(v)) for v in discrepancies])
        remainder=matrix*vector+d
        return float(max(abs(remainder[i]/d[i]) for i in range(2)))
