# D6 stage A: finite moment correction

From the repository root:

```powershell
& code/.venv/Scripts/python.exe -X utf8 code/build.py --study d6
```

Runs seven scientific checks, executes the notebook, writes six CSV datasets,
generates three PNG/SVG figures, and renders the article to HTML and LaTeX source.
Use `--compute-only` to omit rendering or `python -m nscomp.d6` to regenerate data/figures.

This is the source-defined five-moment substep in Eq. (5.16), tested with prescribed
discrepancies and representative bump supports. It is not the full first correction
and does not evaluate the full Navier-Stokes residual.

The API accepts arbitrary supplied d_U,d_E with an explicit positive lambda and
disjoint positive-radius supports. Current notebook experiments do not calibrate
those inputs to a joined paper profile. Parameter changes require revisiting the
normalization, quadrature accuracy, condition numbers, and article's stated baseline.

[Implementation ledger](../../writing/d6_moment_correction/implementation-ledger.md)
records the remaining profile and momentum-residual work.
