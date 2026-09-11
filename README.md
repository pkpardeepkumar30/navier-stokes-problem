# Navier-Stokes computational companion

Reproducible explanations of a proposed forced incompressible Navier-Stokes blow-up construction, with a separate assessment of its physical interpretation.

## First deliverable

**How small does a collapsing vortex become?**

D1 is built, validated, and ready for review. See the [review record](writing/d1_core_size/review-notes.md) for completed checks and rendering limits.

- [Read the article](writing/d1_core_size/article.html)
- [Article source](writing/d1_core_size/article.qmd) and [generated LaTeX](writing/d1_core_size/article.tex)
- [Executed notebook](code/d1_core_size/study.ipynb)
- [Source and calibration notes](writing/d1_core_size/source-notes.md)
- [Parameters](code/d1_core_size/parameters.json), [CSV data](code/d1_core_size/data/), and [validation record](code/d1_core_size/validation.json)

The first study evaluates an anchored scaling surrogate, not a CFD solution. With an illustrative 1 mm / 1 m/s reference pair, the inverse-law limit gives approximately 0.334 nm at 0.01c and 0.0334 nm at 0.1c. These are formal extrapolations; no physical air/water realization has been established, and the displayed fluid-model screens are crossed much earlier.

## Two main folders

```text
code/
  src/nscomp/           reusable Python functions
  tests/               scientific validation checks
  d1_core_size/        notebook, inputs, data, provenance, validation
  build.py             reproducible build
writing/
  d1_core_size/        article, LaTeX source, source notes, figures
```

Future pieces get matching `d2_<topic>/`, `d3_<topic>/`, etc. folders. Repository-wide planning stays in [project-ideas.txt](project-ideas.txt) and the ordered [deliverables.txt](deliverables.txt).

See [code/README.md](code/README.md) for environment setup and the one-command build. Open the HTML file locally for the formatted reading copy; GitHub's file viewer displays HTML source rather than serving the article.
