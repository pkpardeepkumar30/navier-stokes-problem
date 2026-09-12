# D10 reproduction

From the repository root:

```powershell
& code/.venv/Scripts/python.exe -X utf8 code/build.py --study d10
```

The build runs eight scientific checks, regenerates five CSV datasets and four
PNG/SVG figures, executes the notebook from cleared outputs, and renders HTML
and LaTeX source. Use --compute-only to skip Quarto, or -m nscomp.d10 for data/figures.

This is a stationary inviscid Euler equilibrium comparison, not a time-dependent
Navier-Stokes solver. The source notes document normalization, the nonlinear
radial RK4 solve, physical assumptions, and a separate outer-boundary error.
Inputs snapshot D9's reference air state; radius is fixed, not collapsed with speed.
