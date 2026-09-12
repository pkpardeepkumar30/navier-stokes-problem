# D3 review record

Ten scientific checks passed, including high-precision profile integration, independent
physical-coordinate derivatives, a pressure-integral derivative, Cartesian divergence,
the Gaussian vector Laplacian, and the force's pointwise/L2 scaling.
All six notebook code cells executed from cleared outputs.

The HTML and LaTeX exports passed checks for three embedded figures, the computed table
and headline value, local links/internal anchors, and LaTeX figure paths.
All three figures were visually inspected; long numerical-check titles were wrapped
and the final figure inspected again.

The completed paper example is the explicit heat exterior in Eq. (4.29), studied as a
normalized family at positive radius. The inner/annular background and matching constants
remain unevaluated. The Gaussian comparison is separately labelled as a manufactured toy.
This scope is stated in the article, notebook, parameter file, and source notes.

The normalized exterior residual reaches ordinary floating-point precision; the exact
identity, rather than a numerical extrapolation, gives its zero residual. Both absolute
and relative numerical remainders are recorded. No claim about smooth extension of
the full paper force follows from these finite computations.

A complete browser visual review was unavailable. The LaTeX source was exported without
compiling a PDF. See [validation.json](../../code/d3_residual/validation.json).
