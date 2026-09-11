# D1 review record

Prepared for review on 2026-09-11.

- All nine scientific checks passed. The notebook executed all eight code cells from cleared outputs in a fresh kernel, with no error outputs.
- Quarto rendered the HTML reading copy and LaTeX source. Automated checks found three embedded figures, evaluated inline Python expressions, and resolving local article links and LaTeX figure paths. See [the build record](../../code/d1_core_size/validation.json).
- All three generated scientific figures were visually inspected for legibility, annotation placement, units, and consistency with the stated assumptions. The final extreme-speed labels were checked after regeneration.
- The source audit records the paper version, the stricter exponent ceiling, radius/speed definitions, property-table locations, and the incomplete dimensional calibration. See [source notes](source-notes.md).
- Browser access was unavailable, so the complete HTML page has not had a browser-based visual review. Its build and document structure were checked. The LaTeX source has not been compiled to PDF or checked for printed pagination.

The editable narrative is [article.qmd](article.qmd). Review the opening explanation, the illustrative reference pair, and the distinction between formal extrapolation and physical fluid modelling before publication. Preparing these files has not published them.

## Revision following reader feedback

The article now introduces the motivation, question, method, and main physical message before the derivation. It adds a worked calculation at 0.01c, explains molecular averaging and the Mach/Knudsen comparisons in plain language, and moves detailed exponent and viscosity restrictions to a technical note. Figure captions explain how to read the comparisons.

Both HTML and LaTeX were regenerated after the rewrite; all four article code cells executed in each render. The revised exports passed checks for evaluated inline values, three tables, three embedded figures, resolving local links and internal anchors, and valid LaTeX figure paths. The numerical implementation and input parameters are unchanged. The browser/PDF review limitations above still apply.
