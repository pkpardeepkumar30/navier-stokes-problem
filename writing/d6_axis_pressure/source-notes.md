# D6 stage C source and numerical audit

Reviewed 2026-09-12 against the [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
Appendix A.2–A.4, pp.129–134, and the axis condition B.2 on p.144.
Page 130 was also visually inspected. The manifest records the exact PDF hash.

The pressure datum follows A.21 for the ideal azimuthal schedule. The selected
finite constants meet the displayed scalar hierarchy but are not certified against
the complete sufficiently-large/small requirements. No outer cone, complete moment
restoration, or joined flow is claimed.

## Exact schedule represented by the code

Let y=log(X/X_R), f=(1+eta^2)^(-1), and ell=partial_y log(sqrt(2X)*E).
Thus partial_y log(E)=ell-1/2. The smooth step is exactly A.5:

    sigma(s)=exp(-1/s^2)/[exp(-1/s^2)+exp(-1/(1-s)^2)], 0<s<1,

extended by zero and one. Its logit is 1/(1-s)^2-1/s^2. The code evaluates
this expression with a stable logistic function; exact endpoints use their limits.
Symmetry gives integral_0^1 sigma=1/2.

The ideal azimuthal schedule includes, in order:

1. The initial branch E/P_star=f*exp(y/10), y<=0.
2. One unit with ell=(3/5)*(1-sigma), followed by ell=0 for T_d.
3. One unit with ell=-lambda*sigma, then ell=-lambda for 60 log(1/lambda).
4. The axial-pulse interval of length 13/lambda, retaining that azimuthal slope.
5. Profile interpolation over T_f, with
   log(E/P_star)=log(c_start)-(1/2+lambda)y-sigma(y/T_f)log(2)
   -[1-sigma(y/T_f)]log(1+eta^2).
6. A uniform interval of length 30 log(1/lambda), with ell=-lambda.
7. One unit smoothly varying ell from -lambda to -1; a hold of length
   4 log(1/h); one unit varying ell from -1 to -h.
8. A positive wait with ell=-h, selected by the terminal Q target.
9. The three-unit terminal factor f_o=1-rho_o*[1-sigma((y-1)/2)],
   rho_o=c_o*h, then the infinite exterior power with ell=-h.

All coordinates listed for transitions are local to their interval. The axial
pulse velocity and its corrections are not evaluated: they do not enter A.21.
The angular bumps in A.11 are omitted because their pressure integral is zero.
Their interval remains present. The heat replacement in A.39 changes the pressure
integral, but Proposition A.7 restores it; its final datum is therefore the same.

The terminal timing uses the source-prescribed value Q=(lambda-h)/(1-lambda)
at exterior transition start. This value follows after the angular moment
correction in the complete construction. For the ideal schedule, it sets the
unchanged interval lengths; it is not recomputed from uncorrected outer moments.
This distinction is required by Lemma A.5's prescription.

Integrate Q'+(1+ell)Q=-ell-h through the two smooth transitions and the steep hold.
On the final constant-slope wait, Q decays as exp(-(1-h)y).
In A.13, the terminal integrating factor cancels f_o from its denominator, giving

    Q_p=rho_o/(1-rho_o)
        * integral_0^1 exp((1-h)*(1+2s))*sigma'(s) ds.

The waiting length is log(Q_before_wait/Q_p)/(1-h). An independent DOP853
calculation checks the evolution and the near-zero terminal endpoint.
For the selected c_o, the terminal correction f_o'/f_o remains below h/4.
The selected T_f=100 also keeps the additional interpolation slope below 0.1
in magnitude: its largest possible value is 8 log(2)/100.
These checks do not establish all other source inequalities.

## Quadrature, logs, and tail accounting

Compute Pi0/P_star^2 to avoid applying the large amplitude during integration.
The representation of -2Pi0/P_star^2 is a positive sum

    sum_j w_j*(1+eta^2)^(-2theta_j), 0<=theta_j<=1.

Each weight is retained as log(w_j). The initial branch has integral 5*f^2.
Constant-slope finite intervals have exact exponential integrals, evaluated
with expm1. The final infinite tail has integral
c_tail^2/(1+2h). Finite smooth pieces use Gauss–Legendre quadrature, including
the integrated step defining their log amplitude. The terminal interval is
split at y=1, where its factor begins to change.

The baseline uses 128 nodes and compares with 256; orders 8 through 128 give
the refinement data. Separate adaptive integration checks the pressure-dominant
early schedule. Its omitted later tail has an explicit upper bound below 1e-113
in the doubled normalized integral; this bound is far below that check's tolerance.
The exported baseline nevertheless retains all scheduled contributions and the
analytically integrated infinite tail in logarithmic form.

The log10 infinite-tail contribution at eta=0 is about -629.73. Beyond an extra
logarithmic distance L its exact magnitude gains the factor exp(-(1+2h)L).
The tail start has log10(X/X_R) about 537.06. Neither converting that radius to
binary64 nor converting every tiny stage contribution to a direct binary64 number
is part of the computational method.

Ordinary total-pressure evaluation cannot retain terms smaller than its working
precision. Stage logs and decimal-string exports preserve those separate terms.
These decimal strings express values computed from binary64 quadrature logs;
their many printed digits are not independently accurate integration digits.
Quadrature differences are empirical error diagnostics, not rigorous enclosures.

## Derivatives and inner handoff

For g(eta)=(1+eta^2)^(-beta), at center e and with ordinary Taylor coefficients a_n,

    (1+e^2)(n+1)a_(n+1)
      = -2e(n+beta)a_n-(n-1+2beta)a_(n-1), a_(-1)=0.

This follows directly from (1+eta^2)g'=-2beta*eta*g.
It supplies eta coefficients to nscomp.inner through the optional axis_pressure
argument. Without that argument, D6 stage B's original analytic test input remains
the default. Independent high-precision differentiation checks these coefficients.

The handoff fixes eta=0.2, Y=1, Lambda=100P_star^2, j=0.02, sigma=1e-4,
C=1e4 and 70-digit arithmetic. It compares radial degrees 4,8,12,16.
It also repeats the final calculation with the 256-node pressure while holding
the inner parameters fixed, so the measured changes isolate pressure quadrature.
Very small leading equation defects refer to the supplied numerical input.
They do not imply equally small error relative to the exact source integral.

This handoff is local. The global real/complex amplitude bounds and inner exit
condition have not been checked; C=1e4 is not claimed to satisfy B.16.
The independent physical-coordinate test from stage B applies to its original
moderate-amplitude example, not automatically to this new schedule.

## An explicit sufficient bound for B.2

Let A=1/2+h, D=1/2-h, U_star=4eta+j, and H_star=Deta+(1-eta^2)U_star.
Write Z_star=G-d*Pi0'+4Aeta*Pi0, with
G=-A(1-2eta*U_star)U_star-4H_star.
On |eta|<=1, U_max=4+j gives

    |G| <= G_max=A*(1+2U_max)*U_max+4*(D+U_max).

The exact initial branch and nonnegative remaining integrals give
|Pi0|>=2.5P_star^2/(1+eta^2)^2. Also eta*Pi0'>=0.
Both pressure terms in Z_star have sign opposite eta, hence

    |-d*Pi0'+4Aeta*Pi0| >= 2.5*A*P_star^2*|eta|.

Therefore |Z_star|<=delta confines eta to

    |eta| <= epsilon=(delta+G_max)/(2.5*A*P_star^2).

Within that band,

    H_star >= j-(D+4)*epsilon-j*epsilon^2 = H_min.

If epsilon<1 and H_min>0, then
chi>=H_min^2/(H_min^2+sigma^2). With delta=1 and the declared inputs,
epsilon is about 3.17e-15 and this lower bound is about 0.999975.
It exceeds 0.99 with a substantial margin. This sufficient estimate uses the
integral's analytic signs and exact initial contribution; it does not depend on
a floating-point grid resolving the narrow Z_star band. The displayed numerical
constants are evaluations of the algebraic formulas, not a formal proof artifact.

Passing this single axis condition does not certify the other choices in
Appendices A and B, the cone, the exit inequality, the annular match, or the
Section 5 q-power corrections.
