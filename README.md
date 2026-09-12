# Navier-Stokes computational companion

Reproducible explanations of a proposed forced incompressible Navier-Stokes blow-up construction, with a separate assessment of its physical interpretation.

## First deliverable

**How small does a collapsing vortex become?**

D1 is complete and accepted. See the [review record](writing/d1_core_size/review-notes.md) for completed checks and rendering limits.

- [Read the article](writing/d1_core_size/article.html)
- [Article source](writing/d1_core_size/article.qmd) and [generated LaTeX](writing/d1_core_size/article.tex)
- [Executed notebook](code/d1_core_size/study.ipynb)
- [Source and calibration notes](writing/d1_core_size/source-notes.md)
- [Parameters](code/d1_core_size/parameters.json), [CSV data](code/d1_core_size/data/), and [validation record](code/d1_core_size/validation.json)

The first study evaluates an anchored scaling surrogate, not a CFD solution. With an illustrative 1 mm / 1 m/s reference pair, the inverse-law limit gives approximately 0.334 nm at 0.01c and 0.0334 nm at 0.1c. These are formal extrapolations; no physical air/water realization has been established, and the displayed fluid-model screens are crossed much earlier.

## Following installments

| Deliverable | Reading copy | Computational notebook | Status |
| --- | --- | --- | --- |
| D2: finite energy and collapse geometry | [Article](writing/d2_finite_energy/article.html), [animation](writing/d2_finite_energy/animations/core_collapse.html) | [Notebook](code/d2_finite_energy/study.ipynb) | Complete |
| D3: the momentum residual | [Article](writing/d3_residual/article.html) | [Notebook](code/d3_residual/study.ipynb) | Complete for the explicit exterior and toy comparison |
| D4: zero-mean motion and mean stress | [Article](writing/d4_mean_stress/article.html) | [Notebook](code/d4_mean_stress/study.ipynb) | Complete |
| D5: shear growth and viscous decay | [Article](writing/d5_shearing_wave/article.html) | [Notebook](code/d5_shearing_wave/study.ipynb) | Complete for the Kelvin-wave example |
| D6 stage A: finite moment correction | [Article](writing/d6_moment_correction/article.html) | [Notebook](code/d6_moment_correction/study.ipynb) | Stage A complete; full PDE comparison pending |
| D6 stage B: leading inner profile | [Article](writing/d6_inner_profile/article.html) | [Notebook](code/d6_inner_profile/study.ipynb) | Local engine validated; source pressure supplied in stage C |
| D6 stage C: source axis pressure | [Article](writing/d6_axis_pressure/article.html) | [Notebook](code/d6_axis_pressure/study.ipynb) | Ideal-schedule pressure and local handoff; global construction pending |
| D7: conditioning and retained precision | [Article](writing/d7_conditioning/article.html) | [Notebook](code/d7_conditioning/study.ipynb) | Complete for the isolated moment operator |
| D9: physical regime map | [Article](writing/d9_regime_map/article.html) | [Notebook](code/d9_regime_map/study.ipynb) | Complete for declared scale scenarios |
| D10: compressible vortex equilibrium | [Article](writing/d10_compressible_vortex/article.html) | [Notebook](code/d10_compressible_vortex/study.ipynb) | Complete for the stationary Euler comparison |

Each article folder contains source and review notes documenting its checks and limits.

## Two main folders

```text
code/
  src/nscomp/           reusable Python functions
  tests/               scientific validation checks
  d1_core_size/        notebook, inputs, data, provenance, validation
  d2_finite_energy/   same layout for the energy study
  d3_residual/        same layout for the residual study
  d4_mean_stress/     same layout for the stress study
  d5_shearing_wave/   same layout for the wave study
  d6_moment_correction/  finite correction substep and dependency ledger
  d6_inner_profile/  local leading-inner solver and convergence studies
  d6_axis_pressure/  ideal-schedule integral and inner-solver handoff
  d7_conditioning/   precision and sensitivity experiments
  d9_regime_map/     gas/liquid screens and scale audit
  d10_compressible_vortex/  density, pressure, and radial convergence
  build.py             reproducible build
writing/
  d1_core_size/        article, LaTeX source, source notes, figures
  d2_finite_energy/   article, figures, and offline animation
  d3_residual/        article and figures
  d4_mean_stress/     article and figures
  d5_shearing_wave/   article and figures
  d6_moment_correction/  article, figures, and implementation ledger
  d6_inner_profile/  leading-inner article, figures, and scope audit
  d6_axis_pressure/  pressure article, figures, and source audit
  d7_conditioning/   article, figures, and research decision
  d9_regime_map/     article, figures, and model selection
  d10_compressible_vortex/  equilibrium comparison and research decision
```

Each new piece gets matching `dN_<topic>/` folders under code and writing. Repository-wide planning stays in [project-ideas.txt](project-ideas.txt) and the ordered [deliverables.txt](deliverables.txt).

See [code/README.md](code/README.md) for environment setup and the one-command build. Open the HTML file locally for the formatted reading copy; GitHub's file viewer displays HTML source rather than serving the article.
