# D2 source and model audit

This installment explains an energy concentration mechanism. It does not verify the proposed proof or reconstruct its full velocity profile.

## Source version and locations

The paper is the same 166-page version used for D1, with SHA-256
0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f.
Its PDF creation metadata is 2026-09-08 19:06:26 UTC. The source URL, hash, and locations are in
code/d2_finite_energy/source-manifest.json; source access is recorded there.

- Theorem 1.1, page 1: the stated full velocity has uniformly bounded L2 norm and unbounded L-infinity norm approaching the singular time.
- Section 2.1, pages 3-4: a scales as tau^(1/2), b as tau^(1/2-h), and dominant azimuthal/axial speed as tau^(-1/2-h). Radial speed is O(tau^(-1/2)).
- Page 4: core volume scales as tau^(3/2-h) and core energy as tau^(1/2-3h).
- Figure 1, page 4: the displayed difference in shrinking radial and axial scales is explicitly exaggerated.
- Lemma 4.8, page 35: h<exp(-Td), Td=exp(Md)+10 implies h<exp(-10). This is necessary, not sufficient. D1 contains the detailed parameter audit.

Page 4 was rendered locally with pdftoppm and visually checked against the extracted equations.
The source is [OpenAI, *Finite time blowup for Navier-Stokes*](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).

## Energy derivation

For a fixed similarity domain and velocity profile, the change of variables
r=aR, z=bZ gives dV=a^2 b R dR dtheta dZ. The dominant contribution therefore has
the power U^2 a^2 b = tau^(-1-2h) tau^(3/2-h) = tau^(1/2-3h).
The radial velocity contribution is at most order tau^(1/2-h).
All permitted introductory h values have positive exponents; the detailed hierarchy is stricter.
Shape constants and a full profile are not evaluated here.

Bounded kinetic energy means a bounded squared L2 norm for constant positive density.
It does not mean energy is conserved, nor does it bound the L-infinity norm. The forced,
viscous system has work and dissipation terms. This study does not compute that balance.

## Three distinct models

1. Scalar 1D Gaussian: f_epsilon=epsilon^-p exp(-x^2/(2 epsilon^2)).
   Its whole-line squared integral is sqrt(pi) epsilon^(1-2p). It is not a velocity solution.
2. Axisymmetric Gaussian swirl: u_theta=U R exp((1-R^2-Z^2)/2), ur=uz=0.
   In Cartesian coordinates u=(U/a)exp((1-r^2/a^2-z^2/b^2)/2)(-y,x,0).
   This representation is smooth on the axis and divergence free. The peak is U.
   The exact whole-space energy is rho e pi^(3/2) U^2 a^2 b/2.
   It has Gaussian tails, not compact support. Its changing profile is prescribed; the
   time-dependent azimuthal-only velocity does not advect the shrinking envelope into itself.
3. Material deformation: v=(-x/(2 tau),-y/(2 tau),z/tau), tau=1-t.
   With tau0=1 the exact map is diag(sqrt(tau),sqrt(tau),1/tau), determinant one.
   This is a separate affine field with infinite whole-space energy, used only to demonstrate
   volume-preserving parcel deformation. Its trajectories are not superimposed on the swirl.

The Gaussian swirl's energy is a whole-space integral for that toy. The paper comparison is a
leading core scaling, not the full paper energy. Exterior and correction contributions require
separate estimates. None of these examples starts from rest and realizes the proposed smooth force.

## Normalizations and parameter choices

- s=tau/tau0; baseline h=0 is the h->0+ limit, not an admissible construction parameter.
- a0=1, b0=2, U0=1, rho=1 are consistent arbitrary units. b is an axial Gaussian width, not full length.
- Reported ratios are independent of these fixed prefactors.
- The curve at h=exp(-10) is an excluded necessary upper endpoint. Other positive h values
  accepted by the reusable class support mathematical energy-regime examples; they are not automatically admissible.
- At s=10^-12 the ceiling endpoint decreases normalized a/b by about 0.1254 percent.
  For any fixed permitted positive h, a/b tends to zero as s tends to zero; a finite animation
  does not show that entire asymptotic limit.
- No physical viscosity, metres, joules, or collapse duration in seconds is calibrated.

## Numerical checks and animation

Energy is integrated with Gauss-Legendre quadrature on 0<=r<=8a and |z|<=8b.
The exact omitted fraction is erfc(8)+65 exp(-64) erf(8), approximately 1.04e-26.
It is evaluated directly to avoid cancellation in 1-fraction.
Quadrature error reaches floating-point precision and need not decrease monotonically after that.
Cartesian tensor-product quadrature independently checks the cylindrical Jacobian and normalization.
No PDE mesh or time integrator is used.

The animation has 72 frames from s=1 to 0.01, uniform in log(s), played at 12 frames per second
with pause, reset, speed selection, and a frame slider. It does not autoplay.
The spatial field is sampled on a 321 by 321 display grid; analytic formulas, not image pixels,
provide peak and energy labels. Each frame is stored as lossless WebP inside the HTML player;
three representative frames and a static article figure are PNG files.

The left panel has fixed spatial coordinates in initial-radial-width units and fixed logarithmic
colors for |u_theta|/U0, from 1e-4 to the maximum endpoint speed ratio. Values below the color
floor appear in the lowest color. The right panel uses current similarity coordinates and fixed
linear colors for |u_theta|/U(s), from 0 to 1. Axes stretch to fit the panels, as disclosed.
No arrows, streamlines, or particle paths are shown; successive fields are prescribed snapshots.

The animation is an offline standalone HTML file linked from the article. Keep its folder when
sharing the article. Scientific figures and equations remain available in the static reading copy.
