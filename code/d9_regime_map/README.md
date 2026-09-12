# D9 reproduction

From the repository root:

```powershell
& code/.venv/Scripts/python.exe -X utf8 code/build.py --study d9
```

This runs eight scientific checks, regenerates five CSV datasets and four PNG/SVG
figures, clears and executes the notebook, and renders HTML and LaTeX source.
Use --compute-only to skip Quarto, or -m nscomp.d9 for just the data and figures.

parameters.json snapshots the audited D1 reference properties. No full physical
profile, actual pulse wavelength, liquid relaxation time, or CFD mesh is supplied.
See the matching writing folder's source notes and model-selection note.
