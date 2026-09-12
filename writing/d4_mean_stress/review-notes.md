# D4 review record

Completed 2026-09-11.

- Ten scientific checks passed; six notebook code cells executed from cleared outputs.
- HTML and LaTeX source rendered. Three embedded figures, local hyperlinks, internal
  anchors, and LaTeX figure references checked. Inline numbers and the diagnostic table
  inspected in the generated HTML; no replacement characters found.
- All three scientific figures visually inspected. The frequency plot's crowded
  logarithmic tick labels were replaced with the actual integer frequencies and rechecked.
- The velocity image uses physical Cartesian x,y coordinates with equal aspect ratio.
  It is a snapshot, not a flow animation. All baseline units are dimensionless.
- Source PDF page 82 inspected locally, with the paper's cone and averaging
  conventions recorded separately from the toy.

Measured baseline values are in `code/d4_mean_stress/data/summary.json`.
The maximum absolute mean-force error is about 4.3e-14; the coarse-grid stress error
is approximately 0.949. Small numerical divergence does not certify a resolved carrier.

Browser-based visual review was unavailable. The HTML received structural checks,
not a browser screenshot review. LaTeX source was exported; no PDF was compiled.

The complete paper background, pulse dynamics, admissibility verification, and
spatially localized stress construction remain outside this educational experiment.
