# D7 review record

Completed 2026-09-12 for the isolated axial moment operator.

- Seven scientific checks passed. They include independent high-precision adaptive
  integrals, 70-versus-90-digit refinement, equivalent-formulation checks, retained
  versus rounded coefficients, and controlled discrepancy perturbations.
- Six notebook code cells executed from cleared outputs. HTML and LaTeX source
  rendered; three embedded figures, local links, anchors, and LaTeX figure paths checked.
- All three PNG/SVG figure designs were inspected directly. Computed article values
  and UTF-8 decoding were checked after rendering; no unresolved inline code remains.
- The combined D1-D7 scientific regression passed all 59 checks.

At lambda=1e-12, the transformed solve's coefficient relative error is about 6.45e-17,
but its stored float64 coefficients leave an original relative moment remainder of
about 3.05e-4. The reference reevaluates the coefficients as actually represented.
The separate perturbation experiment concerns algebraic sensitivity, not dynamical
stability or admissibility of complete fluid profiles.

The research decision classifies this as an educational numerical-analysis note.
No novelty claim or full correction-level PDE result is established.

Browser visual review was unavailable. HTML received structural checks; figures were
inspected directly. LaTeX source was exported; no PDF was compiled.
