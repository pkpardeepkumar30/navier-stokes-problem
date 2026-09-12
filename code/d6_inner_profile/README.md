# D6 stage B: leading inner profile solver

Compute regular radial coefficients for the source's leading inner equations
using explicit analytic validation axis data.

Run from the repository root:

```powershell
& code/.venv/Scripts/python.exe -X utf8 code/build.py --study d6b
```

Use --compute-only to skip HTML and LaTeX export. To regenerate data and figures
alone, run the same Python interpreter with -m nscomp.d6b.
Dependencies and Quarto setup are in [../README.md](../README.md).

- [study.ipynb](study.ipynb): editable notebook, executed by the build.
- [parameters.json](parameters.json): axis inputs, degrees, precision and scans.
- [data/](data/): four CSV datasets and summary JSON.
- [source-manifest.json](source-manifest.json): source locators and PDF hash.
- [validation.json](validation.json) and [validation-tests.txt](validation-tests.txt):
  latest build and scientific checks.
- [Article](../../writing/d6_inner_profile/article.html) and
  [source notes](../../writing/d6_inner_profile/source-notes.md).

The reusable local engine is nscomp.inner; nscomp.d6b assembles this study.
The main baseline uses radial degree 22 and 80 decimal digits. Its four profile
slices are not a two-dimensional CFD mesh. The independent physical-coordinate
check uses separate moderate-amplitude parameters, recorded in the test and notes.

Radial degree is not the paper's q^(2nh) correction order. Axis pressure is
an analytic test datum, and condition B.2 explicitly fails. The local leading
equations are validated; a fully joined paper profile and first PDE correction
remain absent. The [D6 ledger](../../writing/d6_moment_correction/implementation-ledger.md)
tracks those dependencies.

The CSV fields F_decimal and pressure_increment_decimal are decimal strings:
do not cast them to binary64 when their small magnitudes matter. Recovering a
pressure increment by subtracting two total pressures can erase it even at the
working precision. Inspect pressure_storage.csv for the measured loss.
