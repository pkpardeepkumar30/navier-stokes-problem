# D3 source audit and implementation scope

Status: source audit in progress.

The first evaluated paper field will be the explicit heat exterior in Theorem 4.6(v),
Equation (4.29), page 33 of the version recorded for D1/D2. The profile integral is

    H(Z) = 1/Gamma(1+h) integral_0^infinity exp(-v) v^h (1+Zv)^(-h) dv.

With A=1/2+h and s=r^2/2, the physical exterior swirl is

    u_theta = c_infinity s^(-A) H(2 tau/s),  tau=1-t,
    u_r = u_z = 0.

Here s denotes squared radius/2, NOT the remaining-time ratio used in D2.
The article and code will distinguish these explicitly.

The exterior is independent of z. It solves
partial_t u_theta = partial_rr u_theta + (1/r)partial_r u_theta - u_theta/r^2
at viscosity one. The pressure satisfies partial_r p=u_theta^2/r and p(infinity)=0.
Thus radial centrifugal acceleration and pressure gradient cancel, and the angular
time derivative cancels viscosity.

This formula belongs to the paper's full profile only beyond its matching boundary Xb.
Xb and c_infinity have not been calibrated from the inner/annular construction.
An evaluation with c_infinity=1 and positive h will be a normalized example of this
explicit exterior family on r>0, not an assembled full profile. It is singular at the axis
and independent of z, so it is not a finite-energy whole-space replacement for the core.

The educational contrast will be a prescribed smooth Gaussian swirl with its required
momentum residual. It must be labelled a toy; its manufactured force must not be identified
with the paper's force smooth through the singular time.

## Source locations

- Equation (3.1), page 6: residual at viscosity one.
- Equations (3.2), (4.2)-(4.7), pages 7-8 and 25-26: coordinates, derivative operators,
  Cartesian regularity, incompressibility, and leading centrifugal-pressure balance.
- Equation (4.25), page 33: pressure normalized at radial infinity.
- Equation (4.29), page 33: explicit exterior heat solution.
- Appendix B, pages 144-145: the inner profiles depend on the outer axis-pressure datum
  and coupled analytic construction. They are not specified by the exterior alone.

Source: [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
The source version/hash is the same as the D2 source manifest.

## Numerical plan

Evaluate H and derivatives by generalized Gauss-Laguerre quadrature with weight exp(-v)v^h.
Check selected points with high-precision mpmath integration and finite differences in
physical coordinates. Pressure uses its convergent radial integral with a change of
variables to a finite interval. Keep absolute residual and component-normalized residual
separate, and quantify quadrature convergence before interpreting cancellation.

Primary software references:
[SciPy generalized Laguerre quadrature](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_genlaguerre.html)
and [SciPy Jacobi quadrature](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_jacobi.html).
