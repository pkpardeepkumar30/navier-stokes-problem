# D3: inspect the momentum residual

From the repository root:

~~~powershell
& code/.venv/Scripts/python.exe code/build.py --study d3
~~~

This checks the derivatives and residuals, regenerates four CSV datasets and three
PNG/SVG figures, executes the notebook from cleared outputs, and renders HTML/LaTeX.
The Python environment and Quarto setup are documented in ../README.md.

The reusable implementation is ../src/nscomp/residual.py. Study generation is d3.py.
SciPy supplies generalized Laguerre and Jacobi quadrature; mpmath provides
high-precision reference integration at selected points.

The numerical paper example is the explicit heat exterior in Eq. (4.29). It is evaluated
as a normalized family at positive radius. The inner core, matching annulus, global
localization, and correction fields have not been reconstructed. The Gaussian comparison
is a prescribed toy requiring a manufactured, singular force.

Read the article at ../../writing/d3_residual/article.html and the model/source audit
at ../../writing/d3_residual/source-notes.md. All quantities use dimensionless units,
with viscosity one. Keep pointwise force, component-normalized residual, and the
Gaussian force L2 norm distinct when interpreting the exported data.
