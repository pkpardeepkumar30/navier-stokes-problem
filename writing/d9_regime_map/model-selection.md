# D10 model selection: compressible radial vortex balance

Decision on 2026-09-12: begin with one smooth, steady rotating gas vortex.
The first observable is the density deficit required to support the prescribed
rotation, compared with the pressure obtained when density is held constant.

## Why this is the next bounded experiment

D9's L=a air scenario meets the Mach screen before the gas Knudsen screen.
It makes compressibility a useful first branch, while finer unknown scales prevent
claiming this ordering for the actual paper flow.

Use a fixed 1 mm vortex for the equilibrium benchmark, independently of D1's
contracting radius-speed anchor. At this radius the reference air mean-free-path
ratio is 6.73e-5. Evaluate local density/temperature variation lengths too; do not
infer their validity solely from this radius ratio.

Prescribe the same smooth azimuthal velocity in both balances:

    u_theta(r) = U (r/a) exp[(1-(r/a)^2)/2],  u_r = u_z = 0.

U is its peak speed and a the radius of that peak. This is a new equilibrium
benchmark, not the joined OpenAI profile or a continuation of its forced collapse.

## Models and consistent comparison

Start with stationary inviscid Euler mass, momentum, and total-energy equations.
Axisymmetry and purely azimuthal motion satisfy steady mass and energy transport;
the radial momentum equation is dp/dr = rho*u_theta^2/r.

- Constant-density reference: rho=rho_infinity and p->p_infinity at radial infinity.
- Compressible reference: calorically perfect ideal gas, p=K*rho^gamma, constant
  entropy chosen from the same exterior p_infinity, T_infinity, and rho_infinity.

Compare equilibria with identical velocity profiles and exterior states. These have
different interior densities by design; they are not time evolutions of identical
initial data. No externally supplied force is used in either inviscid model.
Check the mass, both momentum components, and total-energy equation, including
the equation of state and positivity of density and temperature.

This stationary reduction answers an equilibrium question without adding a general
time-dependent compressible CFD solver. It also removes viscosity from both models;
its findings must not be attributed to viscous dynamics or arrest of blowup.
For a later viscous experiment, torque, heating, heat conduction, and consistent
forcing/energy inputs would need specification.

## Implementation and acceptance

1. Integrate radial enthalpy balance independently of the analytic Gaussian integral.
2. Refine the radial quadrature/ODE grid and the outer boundary; report both errors.
3. Sweep peak Mach number using the exterior sound speed, and also report local Mach.
4. Measure center density deficit, normalized pressure discrepancy, and their
   small-Mach scaling. Define a declared 5% density-deficit marker.
5. Audit local gas mean-free-path and thermodynamic-gradient lengths under explicitly
   stated property approximations; exclude states that fail positivity.
6. Produce an executable article, figures, CSV data, and a technical/research decision.

A zero or negative constant-density pressure after fixing the exterior absolute
pressure means that thermodynamic interpretation cannot be carried over to this gas.
It is not an inconsistency in incompressible pressure as a mathematical multiplier.

The expected contribution is educational unless a later focused literature review
and new quantitative result establish more. No manuscript novelty is assumed.
