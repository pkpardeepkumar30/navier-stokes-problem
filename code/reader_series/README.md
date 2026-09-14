# Building the general-public reader series

Use the existing project environment. From the repository root:

```powershell
& code/.venv/Scripts/python.exe code/render_readers.py
```

This regenerates fifteen simple PNG/SVG illustrations and renders all fifteen
public articles to HTML and LaTeX source. It then checks local links, embedded
figures, alternative text, LaTeX figure paths, source hashes and basic editorial
structure. It does not rerun the technical notebooks or imply a new proof audit.

For one article:

```powershell
& code/.venv/Scripts/python.exe code/render_readers.py --study convex-integration
```

Use `--check-only` to inspect existing exports, or `--skip-figures` to render
after a prose-only edit. Rebuild the figure assets when their inputs change.
The figure script is `figures.py`; `figure-data.json` records elementary data,
the reused scientific inputs and their hashes.

The editorial source stays under `writing/<topic>/article.qmd`. No heavy
numerical execution is needed to read or render those explanations.
`manifest.json` defines the main and optional reading routes.
`technical-archive.json` records the originals preserved before the rewrite.
The scientific modules, parameter files and notebooks retain their existing
locations and stated scopes.

The full scientific build `code/build.py --study d4` now regenerates the
technical article and its original checks, then renders the public article
separately. Its `validation.json` identifies the technical source explicitly;
public prose has a separate `reader-validation.json`. LaTeX source is
exported; no PDF compilation or browser-layout verification is implied.
