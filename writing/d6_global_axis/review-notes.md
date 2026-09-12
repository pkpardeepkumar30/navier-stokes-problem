# D6 stage D review record

Reviewed 2026-09-12 for the complex-amplitude and sampled inner-exit component.

- Ten scientific checks passed: exact rational envelope/upward rounding;
  the primitive against independently split high-precision quadrature;
  guarded logarithmic amplitude precision; agreement with the original
  moderate-amplitude solver; the centered Taylor-jet optimization against
  polynomial differentiation; shear diagnostics against radial finite
  differences; persistence of a failed fine-band point under refinement;
  passing critical points and pressure sensitivity; broad/transition exit
  samples; and input validation.
- All 102 project scientific checks passed in 32.302 seconds. This includes
  the earlier studies after the shared inner solver changes. The UTF-8 log is
  retained as `code/d6_global_axis/validation-project-tests.txt`.
- Eight CSV datasets, two summary/bound JSON files, and three PNG/SVG figures
  regenerated. Fine eta coordinates are retained as decimal strings; plot
  coordinates use the separately specified scaled offsets.
- All three figures were visually inspected. The amplitude plot's numeric
  offset was removed for clearer labels and its legend moved away from the
  poles. The fine-band figure retains both signs on symmetric logarithmic
  axes. The exit plot identifies the capped lower bound and its threshold;
  the degree plot uses measured positive defects without a display floor.
- Printed source p.149 was rendered and visually inspected for B.19 and its
  normalizations. The derivation and scope of the separate rational B.16
  envelope are recorded in `source-notes.md`.
- All six notebook code cells executed from cleared outputs without errors.
  HTML and LaTeX source rendered; a final editorial revision defined the
  parameters and identified the convex-integration connection as a
  conceptual comparison. Both formats were rerendered. A text check was
  adjusted to allow ordinary LaTeX line wrapping, then passed.
- Embedded figures, local links/anchors, LaTeX figure paths, computed inline
  values, UTF-8 decoding, and source/parameter hashes passed checks. The
  separate `validation-artifacts.json` also confirms positive Phi and source
  margins at all 306 profile entries, and positive p1 at every positive-Y
  entry. The 51 eta entries contain 50 distinct decimal coordinates because
  the two targeted grids intentionally share their center.
- The denominator-factor bound is exported as `H_factor_modulus_lower`.
  It bounds `|H_star +/- i sigma|`, not Euclidean distance to a pole. The
  metadata label was clarified without changing the bound or numerical
  results, and the project regression suite was rerun afterward.

The selected `Lambda/P_star^2=1e24` case has a fine-band width about
`1.9641e-26`. All 19 coordinates in that band round to one binary64 value.
The comparison at `1e16` has negative angular shear at a shifted point near
the H_star zero. Its failure persists under radial-degree and precision
refinement; the larger scale recovers the signs at the declared points.
The minimum reported baseline exit lower bound is 3, above 2.2. The selected
128/256-node pressure refinement changes the source margin by about
`7.54e-17`.

The analytic amplitude envelope covers the entire stated complex rectangle.
The inner and exit samples do not provide a uniform remainder or a
between-sample certificate. The very small internal equation defects do not
upgrade the finite pressure quadrature's absolute accuracy. These are
leading-inner equations, not the fully corrected Navier–Stokes PDE.

Browser visual review was not performed; scientific figures received direct
visual inspection. LaTeX source is the intended export; no PDF is compiled.
