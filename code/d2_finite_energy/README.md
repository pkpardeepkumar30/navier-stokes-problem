# D2: finite energy and concentration

From the repository root, with the environment described in ../README.md:

~~~powershell
& code/.venv/Scripts/python.exe code/build.py --study d2
~~~

This runs the scientific checks, regenerates the figures and animation, executes the
notebook from cleared outputs, and exports the article to HTML and LaTeX source.
The animation's event handlers are checked with Node when available. That check uses
a DOM stub and does not replace a browser visual review.

For only figures, data, and the animation:

~~~powershell
& code/.venv/Scripts/python.exe -m nscomp.d2 --animation
~~~

Without --animation, the generator updates the static plots and numeric data only.
The article uses that lightweight route. The full build generates the animation's
static frame first, so notebook and article references resolve.

## Files

- study.ipynb: executed, editable explanation and experiments.
- parameters.json: powers, toy normalizations, display sampling, and quadrature choices.
- data/: five CSV datasets plus a summary JSON. Values use arbitrary consistent units or named ratios.
- source-manifest.json: audited source version and page references.
- animation-player.html: player template; generated frames go into the HTML output under writing/.
- animation-manifest.json: frame sampling, coordinates, color ranges, and interpretation.
- check-player.cjs: event-handler checks without external JavaScript packages.
- validation.json and validation-tests.txt: results of the latest full build.

Reusable energy models are in ../src/nscomp/energy.py; study assembly is in d2.py.
The Gaussian swirl and affine deformation are separate models. Neither is the paper's
full base flow. Larger h values accepted by the model class are mathematical scenarios,
not validated construction parameters.

## Reading and sharing

Open ../../writing/d2_finite_energy/article.html for the article.
Open ../../writing/d2_finite_energy/animations/core_collapse.html for the offline player.
Keep the folder structure when sharing the article, so its animation and notes remain reachable.
The player embeds its frames and requires no server or network. It starts paused.
Static PNG/SVG figures and three representative PNG frames are available alongside it.
