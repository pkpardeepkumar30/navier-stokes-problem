# D4: mean stress from zero-mean oscillations

From the repository root:

```powershell
& code/.venv/Scripts/python.exe -X utf8 code/build.py --study d4
```

This runs ten scientific checks, executes the notebook from cleared outputs,
generates six CSV datasets and three PNG/SVG figures, and renders HTML and LaTeX source.
Use `--compute-only` to skip the article render. Export only data and figures with
`python -m nscomp.d4` in the project environment.

The experiment uses explicit periodic fields, spatial averaging, and Fourier
derivatives. It does not integrate Navier-Stokes in time or reconstruct the paper's
stress. The full build records the Python/dependency versions, hashes, and checks.

Edit `parameters.json` for the published baseline. If changing carriers, envelope,
or grid size, update the article's stated choices and verify the velocity/product
bandwidths. Keep intentional under-resolution tests distinct from baseline results.
Notebook exploration can pass alternate values directly without changing published inputs.

Source conventions and analytic checks:
[writing notes](../../writing/d4_mean_stress/source-notes.md).
Reading copy: [article](../../writing/d4_mean_stress/article.html).
