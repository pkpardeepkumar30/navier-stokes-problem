# D6 stage D: source and calculation audit

This component bounds the normalized axis amplitude on an explicit complex
neighborhood and samples the leading inner inequalities. The amplitude bound is
analytic; the inner Taylor calculations and exit checks are numerical evidence
at declared points. Neither constitutes a full joined-profile construction.

## Source and notation

[OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA-256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Appendix B.1, pp.144–145 defines the axis data and the B.2 condition.
Equations B.13–B.16, pp.147–148 give the inner representation and amplitude
choice. Equations B.17–B.19, pp.148–150 give the source, shear, and exit
inequalities. Printed p.149 was rendered and visually checked for the B.19
normalization and its two cases. The pressure input is the stage-C numerical
implementation of Lemma A.5 / A.21.

Use `A=1/2+h`, `D=1/2-h`, `d=1-eta^2`, `L=1-2h eta^2`,
`U_star=4 eta+j` and

$$
H_*(z)=(D+4)z-4z^3+j-jz^2,\quad
\zeta(z)=-\frac{L(z)H_*(z)}{H_*(z)^2+\sigma^2}.
$$

Our normalized angular profile is `F=g Phi`, with `Phi(0,eta)=1` and
`g=exp(Lambda S)/C`. The source uses the unnormalized axis amplitude before
division by `C`. The radial coordinate is `Y=Lambda X`. All quantities are
dimensionless.

## A sufficient B.16 envelope

The following deliberately loose bound is our derivation, not a numerical
maximum or a sharper constant quoted from the paper. For
`0<h<0.01`, `0<j<=0.05`, and `sigma>0`, choose

$$
b=\min(10^{-3},\sigma/68),\qquad R=1+2b,\qquad
\Omega=\{z:|\Re z|<1+b,\ |\Im z|<b\}.
$$

On the closure of this rectangle, `|z|<=R` and

$$
|H_*'(z)|\le 4.5+12R^2+2jR<17.
$$

For `z=x+iy` the vertical segment from `x` lies in the rectangle, so
`|H_*(z)-H_*(x)|<=17b<=sigma/4`. Since `H_*(x)` is real,

$$
|H_*(z)\pm i\sigma|\ge3\sigma/4.
$$

This is a lower bound on the two polynomial factors, rather than a Euclidean
distance to a root. It excludes all poles in the rectangle. Decompose

$$
\zeta=-\frac{L}{2}\left[
  \frac1{H_*-i\sigma}+\frac1{H_*+i\sigma}
\right].
$$

Thus `|zeta|<=B=4(1+2hR^2)/(3sigma)`. The rectangle is convex and contains
zero, so the analytic primitive along the straight path satisfies
`|S(z)|<=RB`. Selecting

$$
\log C\ge\Lambda RB+1
$$

gives `|exp(Lambda S(z))/C|<=exp(-1)` throughout the rectangle, which
supplies the B.16 amplitude domination. The pressure trace is analytic on this
rectangle as well: `1+z^2` has no zero there and admits the branch required
by its fractional powers. Also `L` has no zero there.

`complex_amplitude_bound` uses exact `Fraction` arithmetic on the string
representations of the declared parameter values. The derived binary64
`Lambda` is treated as its declared decimal representation; the solver uses
that same representation. The exported `log_C_upper` is rounded upward
to 40 decimal places. The rational intermediate bounds are retained in
`complex-bound.json`. No interval root solver is needed for this envelope.
This estimate does not quantify the paper's separate fixed-point threshold or
give an inner-series remainder bound.

## Evaluating the exponential with enough digits

The rational integrand has six simple poles for these inputs: the roots of
the cubics `H_*(z)-i sigma` and `H_*(z)+i sigma`. At a root `r` its residue is
`-L(r)/(2H_*'(r))`. Therefore, along the real interval,

$$
S(\eta)=\sum_r \operatorname{Res}(\zeta;r)
                  \log(1-\eta/r).
$$

Each logarithm is continued from zero along the real path; the path crosses
no pole or logarithmic branch cut for these conjugate factors. Conjugate
terms make the sum real. Numerical cubic roots are used for this evaluation
and the pole plot, not for the analytic bound above. A separately split
high-precision quadrature of the original integrand checks the primitive.

The implementation stores `log g=Lambda S-log C` and adds arithmetic guard
digits according to the magnitude of the logarithm before exponentiating.
Computing `S` to 80 digits and then multiplying by a roughly 40-digit
`Lambda` would otherwise lose about 40 absolute digits in the exponent.
The arbitrary-exponent representation avoids binary64 underflow. It does not
make the finite pressure quadrature arbitrarily accurate.

## Targeting the two narrow bands

Let `eta_H` be the central zero of `H_*`. The wider local scale is
`sigma/H_*'(eta_H)`. The pressure-dependent scale used for the finer scan is

$$
\Delta_\eta=\frac{Z_*(\eta_H)}{\Lambda H_*'(\eta_H)},\qquad
Z_*=-A(1-2\eta U_*)U_*-(1-\eta^2)\Pi_*'+4A\eta\Pi_*.
$$

This scale is a diagnostic choice suggested by the coupled source terms,
not a theorem locating the global minimum. The grid combines 17 broad
coordinates, 15 coordinates on the wider scale, and 19 on the finer scale.
The two targeted grids share their center: these are 51 sample entries
at 50 distinct decimal coordinates, with the overlap intentionally retained.
The finest baseline width is about `1.9641e-26` and all 19 coordinates
round to the same binary64 value. Coordinates are constructed at 95-digit
precision, passed to the solver without a float conversion, and exported
as 80-significant-digit decimal strings. Float coordinates in the amplitude
plot are presentation values only. The fine-band plot uses the specified
scaled offsets, not subtraction of rounded float coordinates.

The comparison at `Lambda/P_star^2=1e16` fails angular positivity near,
but not exactly at, `eta_H`. The baseline at `1e24` passes the sampled
signs. Degree/precision refinement at a failed comparison point and a
passing baseline point preserves the observations. The two scans each
use their own pressure-dependent width; their horizontal coordinates are
normalized displacements, not identical absolute eta coordinates.

## Source and exit diagnostics

From the normalized profiles,

$$
p_1=-\frac{2Y\Phi_Y}{\Phi},\qquad
n_s=-2\Lambda U_Y,\qquad
p_2=\sqrt{\frac{Y}{2\Lambda}}\frac{n_s}{g\Phi},
$$

and

$$
S_q=-2L\Lambda\frac{Y\Phi_{YY}+2\Phi_Y}{\Phi}.
$$

The B.17 margin is `S_q-(2.5+0.95 Lambda L chi)`, with
`chi=H_star^2/(H_star^2+sigma^2)`. Angular positivity is checked for
`Y>0`; `p1=0` at the axis `Y=0` by definition. The exit is tested at
`Y=4` with the declared `c_ex=0.2` threshold.

To avoid forming an enormous axial ratio, compute
`log(p2^2/p1)=2 log|p2|-log p1` when `p1>0`. The reported pointwise lower
bound is `p1` for `chi>0.99` and `p1+min(p2^2/p1,3)` otherwise. Both
are no greater than the actual exit expression. A failing angular
positivity test has no reported exit lower bound. The cap is a mathematical
lower-bound device, not a floor inserted into small equation defects.

The leading-equation defects retain stage B's normalization. They omit the
full PDE's axial viscosity and unassembled corrections. Radial-degree
convergence at `eta=0.2,Y=4.1`, two degree/precision checks, and a
128/256-node pressure sensitivity check are recorded separately. The
pressure input still originates in finite binary64 quadrature; the very
small internal defects do not certify comparable absolute source accuracy.

## Remaining scope

The analytic B.16 bound and sampled B.17–B.19 checks have different scopes.
There is no validated bound between eta samples, no rigorous radial
remainder bound, and no certification of the complete outer parameter
hierarchy, corrected moments, or stress cone. Annular joining, the first
q-power correction, and a corrected full-PDE residual remain pending.

The stress-canceling oscillations in the paper's later stages share an idea
with convex integration. This component computes the background axis/inner
profile; it does not implement that oscillatory cycle. The paper's
introduction, pp.2–5, describes both prior convex-integration results and
its dynamical amplification mechanism. This study is not a proof audit.
