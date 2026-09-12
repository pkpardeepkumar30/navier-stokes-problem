# D5: source, coordinate, and validation notes

## Sources and scope

The exact plane-wave family is described by
[Singh and Sridhar, arXiv:1101.5507v2 (2017)](https://arxiv.org/pdf/1101.5507v2),
Equations (1)-(7), pp.1-2. Page 2 was rendered locally and visually inspected.
The source uses U=S*x1*e2. We set x1=y, x2=x, and suppress the third coordinate.
The code's vorticity convention is omega=partial_x(wy)-partial_y(wx).
Its velocity is B*sin(phase), so its real pressure is P*cos(phase).

The [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
cites this reference as item 19. The correspondence is contextual; D5 does not
evaluate its cylindrical pulse equation or calibrate a pulse from its profile.
Relevant source locations are Section 2.2 and Equations (7.5), (7.22).

## Formulas derived for the implemented specialization

    ky = ky0 - S*kx*t, K2=kx^2+ky^2
    I = integral_0^t K2(s) ds
      = t*[kx^2+(ky0-S*kx*t/2)^2+(S*kx*t)^2/12]
    zeta = zeta0*exp(-nu*I), zeta0=initial_speed*sqrt(kx^2+ky0^2)
    B = zeta*(-ky,kx)/K2
    P = 2*S*kx*By/K2
    w = B*sin(kx*x+ky*y), p=P*cos(kx*x+ky*y).

The positive integral expression avoids cancellation in the expanded cubic.
The vorticity amplitude decreases monotonically when nu>0.
The velocity amplitude need not, because |B|=|zeta|/sqrt(K2).

Mean over a complete phase (equivalently a full x period at fixed y):

    E=|B|^2/4, Qxy=Bx*By/2
    production=-S*Qxy, dissipation=nu*zeta^2/2=2*nu*K2*E
    dE/dt = production-dissipation.

The real pressure coefficient is checked by a physical-coordinate momentum residual,
including background transport, perturbation transport, and (w dot grad)U.
No background pressure gradient is required for this affine base flow.

The solution is exact on the unbounded plane. Its background and planar perturbation
do not have finite total whole-space kinetic energy. Phase-averaged energy is an
energy density per unit mass. Uniform mean stress has zero divergence and does not
correct a localized mean residual. The supplied oscillation is initialized at t=0;
no seed-force or smooth cutoff construction is made.

## Units, wavelengths, and resolution

All baseline values are dimensionless. If x=L_ref*xhat, t=T_ref*that,
then velocity=(L_ref/T_ref)*uhat, viscosity=(L_ref^2/T_ref)*nuhat, and
S=(1/T_ref)*Shat. No physical fluid or anchor is assigned.
Normal wavelength is 2*pi/|k| in reference-length units.
Coordinate wavelengths are 2*pi/|kx| and 2*pi/|ky| (infinite when ky=0).

- Baseline S=1, nu=.005, kx=1, ky0=4, initial peak speed=1.
- Analytic history: t in [0,20], 2001 samples, spacing .01.
- Velocity peak: bounded scalar optimization checked against the displayed history;
  these baseline parameters have one interior speed maximum.
- Independent projected-momentum ODE: classical fixed-step RK4 with steps
  .4,.2,.1,.05,.025. No constraint projection or analytic amplitude reset.
- Error metric: maximum Euclidean vector-amplitude difference at matching times.
  The constraint defect is max|kx*Bx+ky*By|, without renormalization.
- Energy-budget quadrature: composite Simpson on each reported time grid.
- Physical momentum check: fourth-order centered first/second differences in
  laboratory x,y,t at multiple points; normalized by the sum of term norms.
- Snapshots: [0,2pi)^2 window, N=256, t=0,4,12, equal axes, fixed
  vorticity/initial-vorticity colors [-1,1]. No advected marker paths.
- For these snapshots, at least 32 samples occur per coordinate wavelength.
  The solution is evaluated directly; there is no numerical spatial evolution
  or assertion of ordinary top/bottom periodic boundary conditions.

Four exported CSV files contain the time history, time refinement, snapshot parameters,
and phase profiles. Snapshot parameters and the explicit formula regenerate the full
images without storing a large repeated image grid. The summary includes parameter hashes.
