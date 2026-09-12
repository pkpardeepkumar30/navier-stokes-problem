# D6 stage B review record

Completed 2026-09-12 for the local leading-inner numerical study.

- Eight scientific checks passed, including physical-coordinate finite
  differences, radial-degree and precision refinement, independently centered
  eta expansions, the large-Lambda comparison, and axis checks.
- Six notebook code cells executed from cleared outputs without errors.
- Four CSV datasets and three PNG/SVG figures regenerated.
- All three figures visually inspected. The radial-degree figure uses integer
  ticks at the evaluated degrees; the large-Lambda title calls its axial term
  a comparison, avoiding confusion with the later time-scale correction.
- HTML and LaTeX source rendered. Embedded figures, local links, anchors, LaTeX
  paths, computed inline values, current article prose, and UTF-8 decoding checked.
- The cached source PDF hash matches the manifest. Source p.147 was visually
  inspected alongside the extracted equations.
- The project-wide suite passed all 83 scientific checks in 50.405 seconds.
  Its output is retained in the matching code folder.

The sampled minimum normalized swirl is about 0.222. Maximum sampled relative
leading angular, axial, and pressure defects are about 2.65e-22, 4.57e-38,
and 1.74e-16 at degree 22. These are local leading-equation diagnostics, not
the full Navier–Stokes residual or rigorous error bounds.

The test inputs explicitly fail B.2: at the identified Z_star zero,
chi is about 0.00927, below the required 0.99. The article states this and
separates axis-data admissibility from local numerical convergence.
It also exposes the loss of tiny pressure increments if total pressures are
subtracted. No fully joined profile or first q-power correction is claimed.

Browser visual review was unavailable; HTML received structural checks and the
figures were inspected directly. LaTeX source was exported; no PDF was compiled.
