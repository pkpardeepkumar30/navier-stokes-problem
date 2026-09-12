# D6 stage D: bounded axis amplitude and sampled exit

Run from the repository root:

```powershell
& code/.venv/Scripts/python.exe code/build.py --study d6d
```

This runs ten component checks, regenerates the data/figures, executes the
notebook from cleared outputs, and renders the HTML article and LaTeX source.
Use `--compute-only` to omit publishing, or `-m nscomp.d6d` for only data and
figures. No network is needed after installing the declared environment.

The project-wide regression command is:

```powershell
& code/.venv/Scripts/python.exe -m unittest discover -s code/tests -v
```

The shared implementation is in `src/nscomp/axis_bounds.py` and `inner.py`;
`d6d.py` assembles this study. The notebook includes an optional local
parameter experiment that does not overwrite the baseline.

## Inputs and outputs

`parameters.json` records the stage-C pressure schedule, two radial scales,
the complex-bound convention, arithmetic precision, grids, and exit
threshold. `log_amplitude_C` overrides the ordinary `amplitude_C` field in
`AxisData`. It is a decimal string, not a floating-point approximation to C.
Critical eta coordinates must remain high-precision values.

Eight CSV datasets in `data/`:

| File | Contents |
| --- | --- |
| `amplitude_real_axis.csv` | Primitive and scaled logarithm at declared eta points |
| `complex_poles.csv` | Six poles and normalized coordinates for the local plot |
| `profile_checks.csv` | 306 profile/source/equation checks on 51 eta entries |
| `exit_checks.csv` | Pointwise lower bounds at the 51 exit entries |
| `critical_scale_scan.csv` | Fine-band signs for both radial scale parameters |
| `degree_convergence.csv` | Radial-degree differences and leading defects |
| `precision_checks.csv` | Two degree/precision refinement comparisons |
| `pressure_sensitivity.csv` | Effect of refining the pressure quadrature |

The targeted grids share one center, giving 50 distinct decimal coordinates
among the 51 entries. `summary.json` supplies computed article values.
`complex-bound.json` keeps
the exact rational envelope and upward-rounded logarithmic normalization.
Figures are exported in PNG and SVG under the matching writing folder.

The analytic complex-neighborhood bound establishes a sufficient B.16
amplitude estimate. The finite inner and exit checks are sampled numerical
results, not uniform B.17–B.19 certificates. See the
[article](../../writing/d6_global_axis/article.html) and
[source audit](../../writing/d6_global_axis/source-notes.md) for derivations
and the remaining construction dependencies.
