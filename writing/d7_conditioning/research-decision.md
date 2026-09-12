# First research decision

Current classification: educational numerical-analysis note with a reproducible
source-specific example. No original research-paper claim is supported yet.

## What is established

The two axial moment weights from Eq. (5.16) become nearly dependent as the
reserved-patch parameter lambda approaches zero. The experiments distinguish:

1. Avoidable coefficient error from forming the original rows in float64.
2. Loss of moment cancellation when large accurate coefficients are stored in float64.
3. Input sensitivity that remains when the equations are evaluated accurately.

An equivalent expm1-based row transformation addresses the first issue. Retained
arbitrary precision addresses arithmetic/storage limits. Neither makes unknown
input data exact or changes the underlying input-output sensitivity.

These conclusions are supported only for the declared operator, bump family,
parameter ranges, and discrepancy vectors.

## Literature and novelty assessment

The paper defines the finite moment operator. Accurate elementary-function evaluation,
row transformations, forward error, and conditioning are standard tools; see
[Higham](https://epubs.siam.org/doi/book/10.1137/1.9780898718027) and
[NumPy's expm1 documentation](https://numpy.org/doc/stable/reference/generated/numpy.expm1.html).
This is a focused contextual check, not an exhaustive literature or priority search.
The current implementation does not establish novelty.

## Evidence needed for a research manuscript

- One fully specified source-consistent leading profile and actual first-order discrepancies.
- An error budget separating profile truncation, quadrature, differentiation, linear
  solve, coefficient storage, and final residual evaluation.
- A demonstrated need for the proposed numerical treatment in that realization,
  with improvements to a meaningful observable beyond this isolated synthetic test.
- Analytic estimates or validated bounds explaining the measured limitation and
  the range where a remedy is reliable.
- A focused literature search around that concrete result before drafting a novelty claim.

The missing profile inputs are recorded in
[D6's implementation ledger](../d6_moment_correction/implementation-ledger.md).

## Implications for later deliverables

D8's independent PDE solver is optional. It would not resolve the algebraic questions
answered here; a benchmark should wait for a specific unresolved PDE question and
an evaluable field/force pair.

D9 can distinguish known dimensional scales from uncalibrated profile and pulse scales.
The toy Kelvin wavelength must not be silently converted into the actual paper pulse scale.
D10's richer-physics study and D11's manuscript remain conditional on their stated
regime-selection and novelty gates.
