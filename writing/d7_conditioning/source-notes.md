# D7 source and accuracy audit

## Source operator and numerical references

- [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
  Eq. (5.16), pp.51-52: axial moment system; same version/hash as D6.
- [NumPy expm1](https://numpy.org/doc/stable/reference/generated/numpy.expm1.html):
  accurate evaluation of exp(x)-1 near zero.
- [mpmath integration](https://mpmath.org/doc/current/calculus/integration.html):
  adaptive integration with arbitrary precision.
- [Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed.](https://epubs.siam.org/doi/book/10.1137/1.9780898718027):
  finite-precision computation, conditioning, and perturbation analysis.

The row transformation is derived within this study. It is not attributed to the
Navier-Stokes paper as one of its construction steps.

## Exact transformation

    original rows: integral R*b_j, integral R^(1-2*lambda)*b_j
    new second row: (old second row-old first row)/(-2*lambda)
    integrand: R*expm1(-2*lambda*log(R))/(-2*lambda)
    transformed rhs: (-d1, (d2-d1)/(2*lambda)).

The first moment equals the bump center exactly by symmetry and unit normalization.
All lambda values are positive. Geometry and support are unchanged from D6.
The expm1 formulation evaluates the row difference before it is lost to rounding.
It does not remove the amplification of input errors through the transformed RHS.

## Independent reference and errors

For each lambda, the original moment matrix is integrated with mpmath at 70 decimal
digits over the local bump coordinate, split into [-1,0] and [0,1].
The normalized continuous bump is integrated afresh; no float64 quadrature nodes
or weights enter this reference. Original decimal parameter strings are used.

The matrix and coefficients retain their working precision. Tests also compare
the smallest-lambda 70-digit solution with a 90-digit solve and compare transformed
float64 quadrature at orders 128 and 256.

Coefficient relative error = Euclidean norm(computed-reference)/norm(reference).
When the computed coefficient is a float64, mpmath receives its actual stored binary
value. Original moment remainder = max_i |(B_ref*alpha+d)_i/d_i|. Both supplied
baseline d_i are nonzero. Moment errors are evaluated in the 70-digit reference
arithmetic, so summation roundoff from a float64 recheck does not conceal stored
coefficient errors.

The retained-reference self-residual is also shown. This is a finite-arithmetic
equation check, not an independent theorem or a mathematical error bound.
Independent reference refinement is covered separately by the 90-digit test.

The original float64 matrix is singular at lambda=1e-18. Its failed coefficient
solve has blank CSV error fields and an explicit status. Nonfinite SVD condition
estimates are retained as inf in CSV and omitted by logarithmic plotting.
At lambda=1e-16 the SVD condition estimate is already nonfinite even though LU returns
coefficients. That solve's coefficient error exceeds one; it is not treated as reliable.

## Controlled experiments

- Accuracy scan: lambda=1e-2 down to 1e-18, fixed discrepancies (1,-.4),
  two disjoint bumps centered at 2,3, halfwidth .12.
- Precision scan: lambda=1e-12, original equations solved at 20,30,40,50 digits.
  Compare retaining coefficients with conversion back to float64.
- Input sensitivity: lambda=1e-10, base discrepancy (1,1), change only its second
  component. Actual input changes are measured after float64 representation.
  Report relative Euclidean coefficient change divided by relative input change.
- No dynamical evolution or perturbation of an assembled full paper profile.
- Log display floor 1e-60 applies only to plots; CSV values retain smaller numbers.

The current stable_axial_solve API returns float64. Its accurate coefficients can
still have large original-moment errors in sensitive cases. Use high_precision_solve
and retain its precision when the original moments require it.
