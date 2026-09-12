# D10 source, normalization, and numerical audit

Reviewed 2026-09-12. The calculation uses a known stationary isentropic Euler vortex,
not the OpenAI joined profile, correction pulses, or a continuation of its singular flow.

## Primary sources and normalization

[Spiegel, Huynh and DeBonis, A Survey of the Isentropic Euler Vortex Problem using
High-Order Methods](https://ntrs.nasa.gov/api/citations/20150018403/downloads/20150018403.pdf):
Sections III-IV, Eqs.16-23, pp.4-6 define the conservative Euler equations, calorically
perfect closure, Gaussian velocity, and isentropic density and pressure.
Source p.5 was visually inspected; p.6 text confirms the nondimensionalization.
The source normalizes speed by exterior sound speed. Its mean-flow Mach M_infinity
is zero for our stationary example; it must not be confused with our peak-swirl M_*.
Its sigma=1, length R=a, and Gaussian strength beta=sqrt(e)*M_* reproduce our field.
We retain D1's gas constant 287.05 rather than the source's 287.15 J/(kg K).

The source studies propagation and numerical methods beyond our radial benchmark.
We do not reproduce its flux-reconstruction solver, grids, periodic boundaries,
or long-time stability experiments.

[NASA isentropic relations](https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html)
supports p/rho^gamma=constant, ideal-gas sound speed, and temperature-density relations.
We derive radial enthalpy balance directly. Streamline stagnation formulas from a
nozzle are not applied across the vortex.

[NIST air reference](https://nvlpubs.nist.gov/nistpubs/jres/110/1/j110-1kim.pdf),
Eq.11 and Table 3, supplies exterior mean free path, viscosity, temperature and pressure.
The temperature/pressure correction is used as a local approximation, including at
states different from the source's reference state. This is not a new measurement.

## Governing equations and thermodynamic choice

Use length a, velocity c_infinity, density rho_infinity, and pressure
rho_infinity*c_infinity^2 as nondimensional units for the conservative flux test.
The article's pressure ratio p/p_infinity is gamma times the code's flux pressure.
Total energy density in the flux units is p_hat/(gamma-1)+rho_hat*|u_hat|^2/2.

Both balances prescribe the same swirl, omit viscosity, and have zero external force.
Their exterior states agree. Interior densities intentionally differ: these are
separate equilibria, not trajectories from identical initial data.

Because radial/axial velocities vanish and all scalars are axisymmetric, steady
mass and energy transport vanish. Azimuthal momentum transport also vanishes.
Radial pressure supplies the centripetal acceleration. The constant-density Euler
reference satisfies the same mechanical balance with density fixed.

The compressible closure is a calorically perfect ideal gas at uniform entropy.
Let K=(gamma-1)*e*M_*^2/2. Then
Theta=1-K*exp(-R^2), rho_ratio=Theta^(1/(gamma-1)), and
p_ratio=Theta^(gamma/(gamma-1)). Positivity requires K<1 throughout the plane.
The routine rejects K>=1 and does not silently clip an inadmissible temperature.

R=0 is handled by the analytic smooth limit. The swirl is proportional to R there.
The constant exterior state and axial invariance mean this benchmark is not a
finite-total-energy whole-space fluid configuration of the kind in the source theorem.

## Observables and error definitions

- Density deficit: 1-rho(0)/rho_infinity; evaluated with log1p/expm1 for small M_*.
- Pressure discrepancy: [p_comp(0)-p_inc(0)]/p_infinity. This is not divided by a
  small pressure deficit; the main scan starts at M_*=0.01.
- Five-percent marker: solve the explicit density relation, independently checked
  with a bracketed root search. It is a declared tolerance for this profile.
- Peak local Mach: maximize u_theta/sqrt(gamma*R_gas*T(r)) over radius.
- Local thermodynamic Kn proxies: lambda*|d log(f)/dr| for f=T,rho,p.
  Since gamma>1 and entropy is uniform, pressure gives the largest of these three.
- Shear rate: |du_theta/dr-u_theta/r|, twice the magnitude of the r-theta strain entry.
  It is computed from the explicit swirl and checked against finite differences.
- Collision/shear ratio: [lambda/sqrt(8*R_gas*T/pi)] times that shear rate.

The local mean-free-path correction is
lambda/lambda_inf=(T/T_inf)*(p_inf/p)*
(1+110.4/T_inf)/(1+110.4/T).
Temperature in kelvin and pressure in a consistent unit are required.
The maxima use bounded scalar optimization over intervals containing the resolved
peaks; radial profiles through R=6 provide the exported values and tail behavior.

The fixed a=1 mm is independent of D1's collapsing radius-speed anchor.
Changing a rescales the local Kn and collision/shear estimates inversely, while
the dimensionless equilibrium profiles at fixed M_* remain unchanged.

## Independent computation and checks

The radial RK4 solver integrates the nonlinear equation
q'=M_*^2*R*exp(1-R^2)*q^(2-gamma), q=rho/rho_infinity, inward from an outer radius.
The analytic density is used only for the optional exact outer datum, not to
update interior values. An alternative boundary fixes q=1 at finite radius.

For the latter, the exact finite-boundary center is
[1-K*(1-exp(-R_out^2))]^(1/(gamma-1)).
Its difference from the infinite-domain center is evaluated with log1p/expm1.
Signed integration and boundary errors can partially cancel; the total-error dip
near R_out=5 is identified as cancellation, not improved convergence.

Eight scientific checks cover an independent adaptive enthalpy integral, EOS and
entropy, all four conservative equations in Cartesian coordinates, radial refinement,
outer-boundary error, small-Mach orders, the density marker and local derivatives,
and the rest limit/positivity rejection.
The four-equation check uses fourth-order centered flux differences. Its exported
refinement scan evaluates a 19 by 19 diagnostic grid; this grid is not a CFD evolution.

Five CSV datasets store radial profiles (including physical proxies), Mach observables,
radial-grid errors, outer-boundary errors, and conservative residual refinement.
No simulation time step, random forcing, molecular dynamics, or dynamical perturbation
growth is computed.
