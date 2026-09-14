# A reader's guide to the OpenAI Navier–Stokes paper

**Start with [the paper map](writing/paper-map/article.html).**
It explains the goal, the main building blocks, how they connect, and where
each section and appendix belongs.

The public articles use arithmetic, graphs, area, and slope to explain the
ideas. Each starts with its question, introduces the concepts through a
simple example, and ends with its role in the whole construction. The
paper's forced setting and the difference between a mathematical model and
a real fluid are stated explicitly.

## Main reading route

| Order | Article | The idea it contributes |
| --- | --- | --- |
| 1 | [How can speed rise while energy falls?](writing/d2_finite_energy/article.html) | Shrinking mass can outweigh growing speed squared. |
| 2 | [What happens at the centre of a vortex?](writing/d6_inner_profile/article.html) | A regular inner shape starts the construction. |
| 3 | [Why centre pressure depends on the outside](writing/d6_axis_pressure/article.html) | The inner and outer flow must fit the same pressure balance. |
| 4 | [The missing push: what a residual means](writing/d3_residual/article.html) | The force calculation identifies the repair target. |
| 5 | [Convex integration: useful fine ripples](writing/convex-integration/article.html) | Fine local slopes can satisfy stronger conditions while preserving a coarse shape. |
| 6 | [Zero average, nonzero momentum transport](writing/d4_mean_stress/article.html) | Opposite fluctuations can produce a useful averaged stress. |
| 7 | [How a small wave grows, then fades](writing/d5_shearing_wave/article.html) | Shear supplies energy; viscosity controls the later decay. |
| 8 | [How local repairs preserve important totals](writing/d6_moment_correction/article.html) | Weighted sums connect local adjustments to global requirements. |
| 9 | [How the pieces become one complete flow](writing/completing-the-flow/article.html) | Scheduling, repeated repair, summation, and localization finish the argument. |

The convex-integration chapter explains both the rapid shear modification
in Appendix C and the related stress-realization idea in the later pulse
construction. It also explains why the viscous dynamics and derivative
estimates still matter.

## Optional side routes

For physical interpretation:

1. [How small does a faster vortex become?](writing/d1_core_size/article.html)
2. [The vortex radius is not the only small length](writing/d9_regime_map/article.html)
3. [What changes when the gas can change density?](writing/d10_compressible_vortex/article.html)

For understanding numerical checks:

- [Why checking a few points is not enough](writing/d6_global_axis/article.html)
- [What computer digits can hide](writing/d7_conditioning/article.html)

These side routes are optional on a first reading of the paper's mechanism.

## Supporting mathematics and computations

All twelve original computational articles are preserved as
`technical-article.qmd`, `technical-article.html` and `technical-article.tex`
in their existing writing folders. Each public article links to its
technical version, notebook and source notes where available.

The notebooks and scientific data retain their stated scopes. They do not
constitute a reconstruction or independent verification of the whole proof.
The public guide explains the remaining blocks without requiring us to
implement a full corrected PDE solution first.

We are currently prioritizing explanation for the general public.
New research, novelty claims and further solver development are paused.
The current priorities and historical computational record are in
[deliverables.txt](deliverables.txt); the topic catalogue is in
[project-ideas.txt](project-ideas.txt).

## Files and rebuilding

`writing/` contains the editable Quarto articles, rendered HTML, LaTeX source,
notes and figures. `code/` contains the notebooks, existing scientific
modules, reader-figure script, and validation records.

Build the public series with:

```powershell
& code/.venv/Scripts/python.exe code/render_readers.py
```

See [code/README.md](code/README.md) and
[reader-series reproduction notes](code/reader_series/README.md).
Open an HTML file locally for the formatted reading copy. LaTeX source is
exported; PDF compilation is a separate step.

The [series review notes](writing/paper-map/review-notes.md) record the
coverage, completed checks and remaining browser/PDF review limits.
