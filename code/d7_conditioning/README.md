# D7: finite moment conditioning

From the repository root:

```powershell
& code/.venv/Scripts/python.exe -X utf8 code/build.py --study d7
```

The build runs seven scientific checks, executes the notebook from cleared outputs,
exports three CSV datasets and three PNG/SVG figures, and renders HTML/LaTeX.
Use `--compute-only` to skip rendering, or `python -m nscomp.d7` for data/figures.

The study uses a 70-digit adaptive-quadrature reference, validated at a selected
extreme point with 90 digits. It compares original and equivalent float64 systems,
coefficient storage, and separately specified input perturbations.

The stable float64 API improves coefficient accuracy; it does not guarantee accurate
original moments after rounding large coefficients. The high-precision API returns
mpmath values that must remain at adequate working precision in subsequent calculations.

[Research decision](../../writing/d7_conditioning/research-decision.md):
an educational numerical-analysis result, with no dynamical-stability or novelty claim.
