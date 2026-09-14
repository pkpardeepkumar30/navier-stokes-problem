# Editorial review of the general-public series

Reviewed 2026-09-12. The fifteen articles are ready for reader feedback.
Start with [the paper map](article.html) and
[the convex-integration chapter](../convex-integration/article.html).

## What changed

All twelve existing public articles were rewritten around a question, a simple
example, and an explicit explanation of their role in the paper. Three new
articles supply the whole-paper map, convex integration, and the completion
of the construction. The main route has nine articles; physical interpretation
and numerical reliability have optional routes.

The intended mathematical vocabulary is arithmetic, school algebra, graphs,
area and slope. Advanced concepts are introduced by their purpose before their
technical names. The numerical reconstruction and proof notation remain in
the separate technical articles. Existing computations retain their original
scope; this revision does not claim new research or an independent proof audit.

## Coverage and explanatory examples

| Article folder | Main teaching device and connection |
| --- | --- |
| paper-map | Six connected blocks, a Sections 1–10/Appendices A–C table, and a reading route. |
| d2_finite_energy | Halving lengths and doubling speed explains how concentration can keep energy bounded. |
| d6_inner_profile | A rotating disk explains the smooth centre; a simple polynomial introduces a local profile. |
| d6_axis_pressure | Adding pressure changes across thin rings connects the centre to the surrounding flow. |
| d3_residual | A cart's missing push introduces the fluid force-balance residual. |
| convex-integration | Fine zigzags, weighted slopes and mean squares connect local shear preparation to useful averaged stress. |
| d4_mean_stress | Opposite fluctuations have zero mean but a nonzero product; spatial variation turns transport into a force contribution. |
| d5_shearing_wave | Moving belts and tilting stripes explain shear, amplification and viscous decay. |
| d6_moment_correction | Two signed adjustments preserve a total while changing its weighted total. |
| completing-the-flow | Scheduling, repairing new errors, adding corrections and confining the flow explain Sections 5–10. |
| d6_global_axis | A narrow dip between sample points explains why local or sampled checks cannot supply a global guarantee. |
| d7_conditioning | Decimal subtraction and almost-identical equations explain lost digits and sensitivity. |
| d1_core_size | A worked inverse scaling compares formal core radii with molecular dimensions. |
| d9_regime_map | Counting molecules in smaller boxes distinguishes core size, gradients and other relevant lengths. |
| d10_compressible_vortex | Three density comparisons show the effect of allowing a stationary gas vortex to change density. |

Each installment except the map has a section named “Where this fits in the
bigger picture.” Captions identify schematic or arithmetic examples, and
reused numerical results are described within their stated assumptions.

## Accuracy points checked against the sources

The map states the paper's claimed **forced** setting, initial rest, smooth
motion before the singular time, bounded energy and unbounded peak speed.
It does not equate this mathematical construction with an attainable
molecular-scale fluid motion. The reviewed PDF and its hash are recorded in
[source-notes.md](source-notes.md) and the reader manifest.

Convex integration is central to the route. The chapter explains the direct
fine-shear construction in Appendix C: acceptable local states average to
the desired coarse state, rapid variation changes slopes substantially while
changing profile values little, and further corrections restore weighted
totals. It states that this step uses a sufficiently large finite frequency
and yields a smooth profile. The later use of waves to supply averaged stress
is explained as a related strategy with additional viscous and smoothness
requirements. The elementary zigzag has corners and is explicitly an analogy.

The wave chapter distinguishes the existing straight-shear numerical example
from the paper's rotating pulse mechanism. The completion chapter distinguishes
background corrections from wave corrections, explains auxiliary coordinates
as an organizing device, and describes why small amplitudes alone do not
control derivatives. A geometric sequence is an analogy for decreasing errors,
not a claim about the paper's exact stage parameters.

Physical examples retain their reference conditions and calibration limits.
The compressible calculation is a stationary comparison at fixed radius;
it does not simulate the collapsing construction or prove that compression
arrests blow-up. Numerical component checks remain distinct from a complete
PDE solution.

## Artifact and build verification

- All fifteen articles rendered to embedded-resource HTML and LaTeX source.
- All fifteen local-link and anchor checks passed, including linked files and
  destination anchors. All fifteen LaTeX figure paths resolve.
- Every article contains its current PNG figure, with alternative text. The
  embedded image bytes match the corresponding local PNG.
- All fifteen PNG illustrations were visually inspected for readable labels,
  clipping, figure/caption agreement and elementary arithmetic. Three small
  figure refinements were regenerated, inspected and rendered again.
- UTF-8 text, unresolved execution markers, navigation, basic article structure
  and a guardrail against unexplained advanced notation were checked. These
  mechanical checks are **not a readability score or a comprehension test**.
- The complete D4 scientific build passed ten checks and executed six notebook
  cells in a fresh kernel. It rendered the technical article with its three
  scientific figures, then the public article separately. This validates the
  revised build routing. The whole historical scientific test suite was not
  rerun for this prose revision.

The current all-article report is
[validation.json](../../code/reader_series/validation.json).
Reproduction instructions are in
[the reader-series README](../../code/reader_series/README.md).

## Preserved material and review limits

The archive manifest records all thirty-six original article files. Hash
comparison confirms that all twelve original Quarto sources remain unchanged
as `technical-article.qmd`. All archived exports still match their original
hashes except the D4 technical HTML, which was regenerated during the build
verification. Its Quarto source and LaTeX output remain unchanged. The archive
manifest intentionally retains the original snapshot hashes.

The browser tool returned “No browser is available,” so complete browser
page layout, responsive behaviour and browser-rendered mathematics have not
been visually checked. The temporary local preview server was stopped. PNG
inspection and structural HTML checks do not replace that review. LaTeX source
was exported; no PDF was compiled or visually reviewed.

The next editorial step is feedback from a reader on whether they can explain
what each block does and why the next block needs it. Research novelty and
further solver development remain paused.
