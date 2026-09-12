# D3 source audit and implementation scope

Status: explicit exterior and Gaussian comparison implemented and validated.

The paper field evaluated here is the explicit heat exterior in Theorem 4.6(v),
Equation (4.29), page 33 of the version recorded for D1/D2. The profile integral is

    H(Z) = 1/Gamma(1+h) integral_0^infinity exp(-v) v^h (1+Zv)^(-h) dv.

With A=1/2+h and s=r^2/2, the physical exterior swirl is

    u_theta = c_infinity s^(-A) H(2 tau/s),  tau=1-t,
    u_r = u_z = 0.

Here s denotes squared radius/2, NOT the remaining-time ratio used in D2.
The article and code distinguish these explicitly.

The exterior is independent of z. It solves
partial_t u_theta = partial_rr u_theta + (1/r)partial_r u_theta - u_theta/r^2
at viscosity one. The pressure satisfies partial_r p=u_theta^2/r and p(infinity)=0.
Thus radial centrifugal acceleration and pressure gradient cancel, and the angular
time derivative cancels viscosity.

This formula belongs to the paper's full profile only beyond its matching boundary Xb.
Xb and c_infinity have not been calibrated from the inner/annular construction.
The evaluation with c_infinity=1 and positive h is a normalized example of this
explicit exterior family on r>0, not an assembled full profile. It is singular at the axis
and independent of z, so it is not a finite-energy whole-space replacement for the core.

The educational contrast is a prescribed smooth Gaussian swirl with its required
momentum residual. It is labelled a toy; its manufactured force is distinct from
the paper's force smooth through the singular time.

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

## Numerical method

Evaluate H and derivatives by generalized Gauss-Laguerre quadrature with weight exp(-v)v^h.
Check selected points with high-precision mpmath integration and finite differences in
physical coordinates. Pressure uses its convergent radial integral with a change of
variables to a finite interval. Keep absolute residual and component-normalized residual
separate, and quantify quadrature convergence before interpreting cancellation.

Primary software references:
[SciPy generalized Laguerre quadrature](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_genlaguerre.html)
and [SciPy Jacobi quadrature](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_jacobi.html).

## Implemented derivatives and pressure

Write sigma=r^2/2, zeta=2 tau/sigma, A=1/2+h. The code calls the profile integral H.

    d^n H/dzeta^n = (-1)^n (h)_n / Gamma(1+h)
        * integral exp(-v) v^(h+n) (1+zeta v)^(-h-n) dv.

Here (h)_n is the rising factorial. At n=0 it equals one. The time derivative is
-2 c sigma^(-A-1) H'. Radial derivatives are assembled by the chain rule with
d sigma/dr=r. The radial vector Laplacian can be written

    c sigma^(-A-1) [2h(1+h)H + (4A+2)zeta H' + 2 zeta^2 H''].

The coefficient 2h(1+h) avoids the cancellation in 2A^2-1/2. The raw sum
u_rr+u_r/r-u/r^2 is exported as a conditioning diagnostic and is less accurate
for this small h. These are algebraically equivalent evaluations.

The exterior pressure integral is

    p = -c^2/2 sigma^(-2A) integral_0^1 v^(2A-1) H(zeta v)^2 dv.

Gauss-Jacobi quadrature handles the weight v^(2h). A separate finite-difference
derivative of this integral is compared with u_theta^2/r. The displayed exact radial
cancellation is imposed analytically; it is not presented as an independent numerical discovery.

## Examples and diagnostics

- Exterior: h=1e-6, c=1, viscosity one. These are normalized exterior-family parameters,
  not a constructed choice satisfying the complete matching/parameter hierarchy.
- Time-scan marker: xi=r^2/(2 tau)=16, at the midplane. We do not know whether the
  chosen xi exceeds the full profile's uncomputed matching boundary Xb.
- Gaussian: h=0 limiting example from D2, a=sqrt(tau), b=2sqrt(tau), U=tau^(-1/2).
  Its exact pressure is -U^2 exp(1-R^2-Z^2)/2. Density does not enter this kinematic pressure.
- The Gaussian sample is R=1, Z=0.5. Its nonzero angular and axial force has pointwise
  magnitude proportional to tau^(-3/2). The plotted L2 norm is integrated on
  0<r<8a, |z|<8b and scales as tau^(-3/4); its quadrature is checked independently
  against a higher-order reference. This is a declared finite domain, not an integral
  of the paper's force or a claim about all exterior/correction contributions.
- For component j, epsilon_j=abs(R_j)/sum(abs(each component term)). All-zero inactive
  components are assigned zero by convention. Componentwise normalization prevents large
  radial terms from hiding angular errors.
- The exterior exact residual is zero. Its computed remainder is quadrature/roundoff error;
  absolute error can grow while relative error remains small along shrinking-radius markers.
- Error values at or below 1e-17 are displayed at 1e-17 on logarithmic plots; raw CSV
  values are retained without that display floor.

The source equation on PDF page 33 was rendered locally and visually checked. The full
inner/annular profile and its nonzero stress remain unevaluated. D3's completed numerical
example is the explicit exterior, supported by a separately labelled manufactured toy.
