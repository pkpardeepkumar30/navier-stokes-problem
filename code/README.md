# Computational studies

The main project is now a general-public reader series. Build its fifteen
articles with `code/.venv/Scripts/python.exe code/render_readers.py`; see
[reader-series notes](reader_series/README.md). The existing computations are
supporting material, and new research/solver development is paused.

The original executable articles are preserved as `technical-article.*`.
The scientific build below checks and regenerates those technical exports,
then renders the corresponding plain-language article separately.

Select a deliverable with --study d1 through --study d7, --study d6b, --study d6c,
--study d6d, --study d9, or --study d10 when running code/build.py. D6b is the
leading-inner numerical study; d6c computes the ideal-schedule axis pressure;
d6d bounds the complex axis amplitude and samples the inner exit conditions.
D6 remains the finite moment-repair study.
The default remains D1. D2's full build also regenerates its offline animation;
see [D2 reproduction notes](d2_finite_energy/README.md).

Run commands from the repository root. The implementation uses Python 3.13; the package declares Python 3.11 or newer. The exact environment is in `requirements-lock.txt` and each study's `validation.json`.

```powershell
py -3.13 -m venv code/.venv
& code/.venv/Scripts/python.exe -m pip install -r code/requirements-lock.txt
& code/.venv/Scripts/python.exe -m pip install -e ./code --no-deps
& code/.venv/Scripts/python.exe code/build.py
```

On Linux/macOS use `python3 -m venv code/.venv` and `code/.venv/bin/python` for the subsequent commands. The frozen environment was verified on Windows; platform-specific dependency availability may differ.

Full publishing also requires [Quarto](https://quarto.org/docs/download/), verified here with version 1.10.18. The build accepts a Quarto installation on PATH, or a portable installation extracted to `code/.tools/` so that `code/.tools/bin/quarto.cmd` exists. Quarto is an external tool, not installed by pip. No LaTeX engine is needed to produce the HTML article and `.tex` source. Compile `.tex` separately if a PDF is wanted.

The portable Windows archive used during development is `quarto-1.10.18-win.zip`, SHA-256 `4e824652ff0da3f646868277582ed59c0872d1456e35350b7d7cdc4243ee18c2`, from the official Quarto GitHub release. Downloaded tools, source documents and Python environments stay outside Git.

For calculations, scientific checks and the executed notebook without publishing:

```powershell
& code/.venv/Scripts/python.exe code/build.py --compute-only
```

For only the exported D1 data and figures:

```powershell
& code/.venv/Scripts/python.exe -m nscomp.d1
```

## Layout and responsibilities

| Path | Contents |
| --- | --- |
| `src/nscomp/` | Reusable calculations; `scaling.py` contains the power-law API and `d1.py` assembles this study |
| `d1_core_size/parameters.json` | Explicit input values, assumptions and source locators |
| `d1_core_size/study.ipynb` | Editable notebook, saved with outputs after a clean execution |
| `d1_core_size/data/` | Generated CSV tables, exact plotting grids, and summary JSON |
| `d1_core_size/source-manifest.json` | Audited source URLs, metadata and hashes |
| `d1_core_size/validation.json` | Result of the latest build, including whether publication formats were rendered |
| `tests/` | Scientific checks, including independent time parameterization and units |
| `build.py` | Runs scientific checks and the notebook, renders the technical article, then checks the public article separately |
| `render_readers.py` | Renders the public reader series without rerunning numerical studies |

Add future studies as `d2_<topic>/`, `d3_<topic>/`, etc. Put a routine in `src/nscomp/` only when it has a clear reusable responsibility; keep study-specific parameters and outputs with the study. Do not create top-level folders for each new post.

The notebook and article use the same numerical functions. Figures are generated in the matching `writing/dN_<topic>/figures/` folder, so each written piece carries its own publication assets. D1 and D9 use SI reference values; the intermediate mathematical studies use explicitly dimensionless examples. Consult each study's parameter file for units and scope.

## Changing the experiment

Use the notebook's alternative-calibration cell to explore without changing the published baseline. `CoreScaling(r0_m=..., u0_m_s=..., h=...)` can also be imported directly. The class accepts h=0 and h=0.01 only to evaluate limiting cases; it does not validate the full proof hierarchy.

`parameters.json` records the baseline used in this article. If changing its reference states, update source notes and any fixed narrative assumptions as well as rerunning the build. Exponent/prefactor sweeps are scenarios, not confidence intervals or constructed Navier-Stokes flows. A complete dimensional profile/viscosity match remains outside D1.

## Reproducibility boundaries

The calculation and ordinary build need no network after dependencies and Quarto are available. Source inputs are retained as small, cited numerical records rather than requiring a fresh web download. The embedded HTML contains its figures and styles; external hyperlinks still need connectivity, and mathematical rendering may use Quarto's configured MathJax resources.

The first deliverable has no CFD time step, spatial mesh, randomness, or molecular simulation. Its validation concerns the algebra and numerical evaluation of the explicit surrogate. Numerical roundoff is small in the displayed range; uncertainty in physical calibration is a separate issue.
