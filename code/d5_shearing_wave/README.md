# D5: one Kelvin shearing wave

From the repository root:

```powershell
& code/.venv/Scripts/python.exe -X utf8 code/build.py --study d5
```

The full build runs eight scientific checks, clears and executes the notebook,
generates four CSV datasets and three PNG/SVG figures, and exports HTML/LaTeX.
Add `--compute-only` to skip article rendering, or run `python -m nscomp.d5`
in the environment to regenerate just data and figures.

The main curves use an analytic solution on an unbounded affine shear. A separate
fixed-step RK4 calculation verifies its amplitude evolution. Neither is the
OpenAI paper's full pulse. No physical fluid, box boundaries, or finite whole-space
energy are assigned.

The published choices live in `parameters.json`. Changing them can change whether
there is an interior maximum, and requires checking peak selection, narrative values,
time interval, and snapshot wavelength sampling. The reusable KelvinWave class also
supports the zero-shear heat reference and zero-viscosity comparison.

[Read the article](../../writing/d5_shearing_wave/article.html) and
[numerical/source notes](../../writing/d5_shearing_wave/source-notes.md).
