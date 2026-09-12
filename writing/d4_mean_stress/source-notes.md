# D4 source and numerical audit

Status: implemented and validated; see review-notes.md and the matching validation record.

## Source map

[OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
same cached version and SHA-256 as D1-D3. Source page 82 was rendered locally and inspected.

- Section 2.2, pp.5-6: oscillatory momentum transport.
- Section 3.2, pp.9-10: cylindrical stress divergence.
- Equations (4.20), (4.26), pp.30-33: shear variables and cone margin.
- Equations (7.23)-(7.26), p.82: angular/auxiliary averages at fixed slow variables,
  before physical phase evaluation; positive squared pulse amplitudes.

The explicit cone test is

    n_theta + t_s n_z >= kappa,
    (v_s-2)(n_z-t_s n_theta)^2 <= (2-kappa)(n_theta+t_s n_z)^2,

with n=T0/|T0|, t_s=-b_s/a, v_s=a(1+t_s^2), v_s>2, and 0<kappa<2.
The a,b_s here are paper shear variables, unrelated to our envelope a(y).
We record this map without evaluating it. Our Cartesian covariance is not a reconstructed
paper flux vector or a test of its admissibility.

## Exact toy identities

All examples are dimensionless on the flat square torus [0,2pi)^2, without boundaries
or time integration. The velocity is smooth and periodic. These examples are derived
within this study; no profile is fitted to paper data.

Constant example: average over the entire torus, with wavevectors (n,-n), (m,m),
polarizations (1,1), (1,-1), and amplitudes sqrt(lambda_plus), sqrt(lambda_minus).
The normalized covariance is half their weighted outer-product sum. Its eigenvalues
are lambda_plus and lambda_minus. This family does not realize every covariance orientation.

Envelope example: average over x at each y, psi=a(y)cos(k(x+m*y))/k,
a=1+epsilon*cos(y). Use w=(psi_y,-psi_x), retaining a'/k in wx.
The covariance and force formulas are derived explicitly in the article. Pressure is
kinematic pressure (divided by density), p=-a^2/2. The restricted periodic Leray projection
removes nonconstant vertical force components and preserves the zero mode.
The force integral is zero in this example. Nonzero curl establishes that the remaining
horizontal force cannot be absorbed into pressure.

The averaged momentum equation applies to a hypothetical mean-plus-fluctuation field;
the prescribed snapshot is not asserted to be a steady unforced solution.

## Numerical conventions

- Python/NumPy float64. Grid coordinates 2*pi*j/N, j=0,...,N-1; no repeated endpoint.
- Derivatives use FFT multipliers i*2*pi*fftfreq(N,d=2*pi/N).
  [NumPy frequency convention](https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html).
  The even-grid Nyquist derivative is zero, the conventional real interpolant choice.
- N=128 baseline. For the constant example, max carrier coordinate frequency is 9,
  max product frequency 18. Strictly sub-Nyquist products require N>36.
  The envelope has coordinate bandwidth k in x and |m|k+1 in y; its products have
  at most twice these frequencies. N=256 resolves every displayed k up to 32.
- N=16 deliberately loses the frequency-eight sine in the constant example. N=20
  recovers the mean without resolving all products. Divergence and covariance alone
  cannot certify a field's full bandwidth.
- Stress relative error uses the Frobenius norm. Other reported errors are maximum
  absolute component errors. The unit test also normalizes divergence by the carrier
  frequency times the sampled maximum speed, with tolerance 1e-13.
- Roundoff grows with differentiation and accumulation; further refinement need not
  decrease errors already close to floating-point precision.
- Figure 3 floors relative errors at 1e-17 only for log display; exports are unfloored.
- The k^-2 covariance difference is an exact finite-frequency contribution, not
  discretization error. Its predicted factor-four decrease is checked against sampled
  products on the separate resolved grid.

## Exported data

Six CSV files: scalar sine samples, constant covariance, spatial stress/force,
full baseline envelope velocity grid, resolution scan, and frequency scan.
The summary JSON contains actual measured errors and parameter-file hash.
The plotted scalar curves append their exact periodic endpoint only for display.
