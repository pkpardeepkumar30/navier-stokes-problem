# D5 review record

Completed 2026-09-11.

- Eight scientific checks passed, including primitive momentum evaluated with independent
  finite differences and fourth-order convergence of the amplitude integrator.
- Six notebook code cells executed from cleared outputs.
- Three PNG/SVG figures visually inspected. Crowded logarithmic time-step labels
  were replaced by the actual step sizes, regenerated, and rechecked.
- HTML and LaTeX source rendered. Embedded figures, inline values, local links,
  internal anchors, and LaTeX figure paths checked. No replacement characters found.
- The source's coordinate and pressure conventions were checked on rendered PDF page 2.

The baseline reaches a speed ratio of 3.6326 and an energy ratio of 13.1959 near
t=3.995, followed by decay. The finest RK4 amplitude error is about 1.32e-8.
The source/parameter hashes and complete computed values accompany the code.

HTML received structural checks; browser visual review was unavailable.
The scientific figures were inspected directly. LaTeX was exported without compiling a PDF.

The figures describe an unbounded affine-flow example with phase-averaged perturbation
energy. They do not represent the complete paper pulse, a finite-energy whole-space flow,
particle trajectories, or a fluid-specific dimensional calibration.
