# Articles and manuscripts

Use a matching folder for every deliverable, such as `d1_core_size/`. Keep the editable narrative, source notes, generated article/LaTeX source, and publication figures together. Numerical code, notebooks and data belong in the matching folder under `code/`.

Each article follows this convention:

- `article.qmd` is the canonical editable article. Its Python cells call `nscomp` to compute numbers and tables and regenerate figures.
- `article.html` is the rendered reading copy.
- `article.tex` is the generated LaTeX source; edit the `.qmd` rather than maintaining two independent texts.
- D1's `article.css` supplies the shared HTML styling.
- `source-notes.md` records notation, provenance, dimensional calibration and limitations.
- `review-notes.md` records completed checks and remaining visual-review limitations.
- `figures/` holds generated PNG and editable-text SVG figures.

Run `code/build.py` with the project environment from the repository root to rebuild. A later manuscript can have its own LaTeX source in its deliverable folder when it needs a different structure; it should still reuse the computational outputs rather than copy calculations into the prose.

Rendered articles are local review artifacts. This repository does not automatically post them to a blog or submit manuscripts.
