# D6 stage C review record

Reviewed 2026-09-12 for the ideal-schedule axis-pressure component.

- Nine scientific checks passed: smooth-step identities, independent adaptive
  pressure quadrature, quadrature refinement, an independent terminal Q ODE,
  pressure signs/bounds and derivatives, eta coefficients checked by direct
  high-precision differentiation, the inner-solver handoff, an algebraic B.2
  bound and radius-scale check, and input validation.
- The project-wide suite passed all 92 scientific checks in 29.307 seconds.
  Its UTF-8 output is retained in the matching code folder.
- Six notebook code cells executed from cleared outputs without errors.
- Seven CSV datasets and three PNG/SVG figures regenerated.
- The figures were inspected for labels, legends and interpretation. The
  contribution plot retains logarithmic values below binary64 range. The final
  handoff plot uses the measured positive defects directly, with no artificial
  floor; its axial curve reaches working roundoff.
- HTML and LaTeX source rendered. Embedded figures, local links, anchors and
  LaTeX figure paths checked by the build. Computed inline values and UTF-8
  decoding checked separately.
- Source p.130 was visually inspected; the pressure-preservation argument in
  Proposition A.7 was checked on pp.139–140.
- Stage B was rebuilt after adding the optional pressure adapter. Its original
  eight checks and six notebook cells still pass, and its article links to stage C.

The normalized center pressure is about -3.314622730; the exact initial branch
contributes 75.4% of its magnitude. The infinite-tail contribution has base-10
logarithm about -629.73. The terminal logarithmic radius corresponds to a ratio
about 10^537.06, which is not formed directly in binary64.

The 128/256-node normalized pressure difference is about 8.9e-16 on the declared
grid. The very small inner equation defects measure consistency with this
numerical input. Its propagated quadrature changes are stored separately:
about 3.4e-19 in Phi and 4.5e-18 in U at the selected point.

The analytic pressure bound gives a sufficient B.2 estimate with chi at least
about 0.999975 in the required narrow band. This does not certify the full outer
hierarchy, moment/stress construction, global complex amplitude bound, inner
exit inequality, annular match, or first q-power correction.

Browser visual review was unavailable; HTML received structural checks and the
figures were inspected directly. LaTeX source was exported; no PDF was compiled.
