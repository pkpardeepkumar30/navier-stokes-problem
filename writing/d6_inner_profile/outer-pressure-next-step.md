# Next implementation: the source's axis pressure

Source review: 2026-09-12. This specification was subsequently implemented in the
[D6 stage-C pressure study](../d6_axis_pressure/article.html) for declared finite
schedule parameters, with a local inner-solver handoff and a sufficient B.2 bound.
The full outer hierarchy and global inner/joining requirements remain open.
The text below records the component's intended scope and acceptance checks.
It narrows the next D6 step using Appendix A.2–A.4, especially Lemma A.5 and
Eqs.A.21–A.23 of the [paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).

The axis pressure can be evaluated before implementing every later correction.
Lemma A.5 defines it from the ideal azimuthal schedule with the angular correction
bumps omitted. Those bumps preserve its total pressure integral. The later heat
replacement and compensation also restore that integral.

With y=log(X/X_R), the exact defining formula is

    Pi0(eta)=-1/2 * integral_-infinity^infinity E_id,sched(y,eta)^2 dy.

It is independent of X_R. The schedule has the form c(y)*f(eta)^theta(y),
where f=(1+eta^2)^(-1) and 0<=theta<=1. Its initial branch is
P_star*f*exp(y/10) for y<=0. That branch alone contributes
-5*P_star^2*f^2/2; it is not the whole pressure trace.

The implementation should use logarithmic radius and normalized pressure Pi0/P_star^2.
It should integrate the semi-infinite initial and terminal powers analytically,
and each finite transition in its own coordinate. Huge radial lengths must not be
turned into a uniform physical-radius mesh. The source's smooth step, transition
lengths and stopping condition in A.12–A.13 must be retained explicitly.

The input ledger must distinguish evaluated constants from existence thresholds.
The conditions P_star>exp(T_d) and h<exp(-T_d) do not by themselves certify all
the "sufficiently large/small" choices in Appendix A. These scales can also exceed
the practical range of a direct amplitude representation even when normalized
shape integrals remain computable.

Acceptance checks for this next numerical component:

- Recover the exact initial-branch pressure and its eta derivatives.
- Verify independence from X_R, evenness, negative pressure, and the strict
  derivative sign for eta!=0, with quadrature and truncation errors reported.
- Compare direct quadrature derivatives with the analytic f^(2theta) derivatives.
- Account separately for omitted tail bounds and numerical roundoff; dropping
  an exponentially small term is an approximation that needs an explicit bound.
- Supply eta Taylor coefficients to the inner recurrence, then repeat its
  convergence and physical-coordinate checks.
- Evaluate H_star and Z_star and choose delta and sigma against condition B.2
  before describing the resulting axis data as suitable for that joining step.

This component alone will not construct the axial pulse, correct all outer
moments, verify the stress cone, certify a complex-domain amplitude bound, or
compute the first q-power correction. Those remain separately tracked in the
[D6 implementation ledger](../d6_moment_correction/implementation-ledger.md).
