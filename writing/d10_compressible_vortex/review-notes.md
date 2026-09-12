# D10 review record

Completed 2026-09-12 for the stationary equilibrium comparison.

- Eight scientific checks passed, including independent enthalpy quadrature,
  all four conservative Euler equations, EOS/entropy, radial convergence,
  outer-boundary error, small-Mach orders, and local derivative checks.
- Six notebook code cells executed from cleared outputs.
- Five CSV datasets and four PNG/SVG figures regenerated.
- All four figures visually inspected. Radial-refinement ticks were replaced
  with the actual interval counts; the gas-screen legend was moved clear of its line.
  The outer-boundary figure identifies cancellation between signed errors.
- HTML and LaTeX source rendered. Four embedded figures, local links, anchors,
  LaTeX figure paths, the computed table and inline values, and UTF-8 decoding checked.
- The project-wide scientific suite passed all 75 checks.

The chosen vortex reaches a 5% center density deficit at M_*=0.193274, about
66.68 m/s for the exterior reference state. Its finite radial solver converges
at fourth order, reaching a maximum normalized density error about 4.60e-12.
These results concern a known stationary inviscid Euler family, not a time-dependent
compressible Navier-Stokes solution or a test of blowup arrest.

Browser visual review was unavailable; HTML received structural checks and figures
were inspected directly. LaTeX source was exported; no PDF was compiled.
