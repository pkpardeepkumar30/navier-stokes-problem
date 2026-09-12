# D6 stage A: source and numerical conventions

Primary source:
[OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
same version/hash as earlier studies. PDF p.52 was rendered and visually checked.
Relevant locators: Eqs. (5.1)-(5.7), pp.46-48; moments (5.10)-(5.11), p.49;
repair (5.14)-(5.16), pp.51-52.

The scoped computation implements the finite moment solve. Source-to-code mapping:

    R=sqrt(2X), ef=e_star*f(eta), f=1/(1+eta^2)
    BU rows: integrals of R*b_U and R^(1-2*lambda)*b_U
    BE rows: integrals of R^2*b_E, R^(-2-2*lambda)*b_E, R^(-2*lambda)*b_E
    alpha=solve(BU,-d_U), beta=solve(BE,-d_E)
    d_U=(m1,m4/ef), d_E=(m2,m3/(2*ef),-m5/ef).

The corrections to the five original moments are

    delta m = ((BU alpha)_1, (BE beta)_1, 2 ef (BE beta)_2,
               ef (BU alpha)_2, -ef (BE beta)_3).

These are positive-order profile corrections with signed coefficients. They are
different from the positive squared amplitudes in the pulse covariance construction.
The patch exponent lambda differs from the correction-order power lambda_n=2*n*h.

## Representative inputs

All coordinates and coefficients are dimensionless. Lambda=.1, ef=1 at eta=0,
two axial bumps at R=2,3 and three angular bumps at R=4,5,6, all with halfwidth .12.
The input discrepancies are prescribed test data. The support locations have not
been embedded in the source's actual uncomputed annulus.

The bump is exp(-1/(1-z^2)) on |z|<1, zero otherwise, normalized to unit integral
in R. Disjoint supports stay away from the axis. Profile values use a fixed
128-point normalization; moment quadrature normalizes at the specified order.
These normalizations converge to the same continuous bump. Tests compare against
adaptive quadrature and verify support and all five independently weighted changes.

## Numerical checks and limits

- Float64 Gauss-Legendre quadrature on each bump's local coordinate z.
- Baseline order 96; independent-order reevaluation at 384.
- Refinement orders 12,24,48,96,192. Coefficients are recomputed at each order.
- Direct linear solves are used instead of explicit inverses.
- Maximum relative moment remainder means max_i |m_after_i|/|m_before_i|;
  all five baseline denominators are nonzero. No physical sum of unlike moments is used.
- Matrix conditioning is the 2-norm condition number in the declared coordinates.
  The U block is also row-equilibrated before a second condition estimate.
- Lambda scan .1 down to 1e-8 holds bump geometry and test discrepancies fixed.
  It is algebraic sensitivity, not flow stability or an admissible full-profile family.
- Logarithmic plots floor zero remainders at 1e-17; raw CSV values remain unchanged.
- Seven tests include adaptive integration of reserved-patch pressure/momentum factors,
  not only matrix multiplication with the matrix used in the solve.

No physical-time scan, completed first coefficient, or full PDE residual reduction is
reported. See implementation-ledger.md for dependencies.
