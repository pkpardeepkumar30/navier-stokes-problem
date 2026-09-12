# D6 stage B: source equations and numerical scope

Reviewed 2026-09-12 against the [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
The source manifest records the exact PDF hash. Pages 25–29, 35–37 and 144–149
were read; page 147 was also inspected as a rendered image.

This study evaluates the leading inner equations with declared analytic test data.
The pressure trace is not supplied by an assembled outer construction, and one
axis-data selection condition explicitly fails. This is a numerical method
validation, not a joined realization or independent verification of the theorem.

## Variables and leading balances

Write A=1/2+h, D=1/2-h, d=1-eta^2, and L=1-2h*eta^2.
The physical coordinates obey

    1-t=q*(1-eta^2), z=q^D*eta, X=r^2/(2q).
    u_theta=q^(-A)*sqrt(2X)*F, u_z=q^(-A)*U,
    r*u_r=V0, p=q^(-2A)*Pi.

Source Eqs.4.3–4.7 impose axis regularity on F, U, Pi, and V0/X.
In Appendix B the source writes F=phi/C. We use Y=Lambda*X,
g=phi_axis/C, and Phi=F/g. Phi(0,eta)=1.

Define avg(U)=(1/Y)*integral_0^Y U(s,eta) ds with its smooth value at Y=0,

    W=1-2D*eta*avg(U)-d*partial_eta(avg(U)),
    H=D*eta+d*U,
    V0/X=(2eta*U+W-1)/L,
    xi=partial_eta(log(g)).

Multiplying out the logarithmic terms of Eqs.4.13/B.15 gives

    2L*Lambda*(Y*Phi_YY+2Phi_Y)
      = W*(Phi+Y*Phi_Y)+h*(1-2eta*U)*Phi+H*(Phi_eta+xi*Phi),

    2L*Lambda*(Y*U_YY+U_Y)
      = W*Y*U_Y+A*(1-2eta*U)*U+H*U_eta
        +d*Pi_eta-4A*eta*Pi-2eta*Y*Pi_Y,

    Lambda*Pi_Y=g^2*Phi^2.

All derivatives in these equations hold the other similarity coordinate fixed.
The pressure relation is the leading centrifugal balance. The two tangential
balances include radial viscosity with viscosity normalized to one. They exclude
axial viscosity, as in the source's R0 definition. The remaining radial momentum
balance also belongs to the subsequent correction analysis. These three equations
must not be called the full Navier–Stokes residual.

## Coefficient recurrence and derivative accounting

For ordinary radial coefficients f=sum f_n Y^n, use convolution for products.
The coefficient of Y^n in Y*f_YY+2*f_Y is (n+1)(n+2)*f_(n+1);
in Y*f_YY+f_Y it is (n+1)^2*f_(n+1).

At step n, form the two right-hand-side coefficients from indices 0 through n,
then divide by 2L*Lambda*(n+1)(n+2) or 2L*Lambda*(n+1)^2.
The pressure coefficient is

    Pi_(n+1)=g^2 * sum_(j=0)^n Phi_j*Phi_(n-j) / [Lambda*(n+1)].

The average has coefficients U_n/(n+1). These formulas give an explicit
triangular radial recurrence after the analytic axis values are supplied.
They implement the radial inverses described in Eq.B.5.

Each radial coefficient is represented by ordinary Taylor coefficients in eta
about one declared center. Multiplication, reciprocals and differentiation are
performed on these coefficients. The initial eta degree is radial degree+4.
Each radial step consumes at most one eta derivative, leaving the derivatives
needed for center evaluation and leading-balance checks. Artificial zero padding
at the top of a derivative array does not certify those highest coefficients.
Nonzero-offset evaluation is a local polynomial approximation, validated only
by a small overlap comparison and the declared local physical-coordinate test.
It is not a globally accurate interpolant over all eta.

The axis amplitude uses Eq.B.3:

    U_star=4eta+j, H_star=D*eta+d*U_star,
    zeta=-L*H_star/(H_star^2+sigma^2), xi=Lambda*zeta,
    g(eta)=exp(Lambda*integral_0^eta zeta(w)dw)/C.

The axis value uses arbitrary-precision quadrature. Its eta coefficients follow
g_eta=xi*g by convolution, avoiding numerical high-order differentiation.
The baseline pressure is Pi_axis=-100/(1+eta^2)^2. It is even and negative
with the stated derivative sign pattern; this is not a derivation from Lemma 4.8.

## Numerical checks and observables

The main profile scan uses degree 22, 80 decimal digits, Y in [0,4.1] with 61
samples, and eta=-0.6, 0, 0.05, 0.2. Residuals are reevaluated from the polynomials,
including nonlinear products beyond the degree enforced by the recurrence.
Relative defect means abs(sum of equation terms)/sum(abs(terms)); if all terms
vanish it is zero. No fixed floor replaces a small denominator.

The degree scan truncates a degree-24 series at degrees 4,6,...,22 and compares
values and defects at eta=0.2, Y=4.1. The degree-24 value is a numerical reference,
not an exact solution or rigorous remainder enclosure.

Eight scientific checks cover independently derived first coefficients, degree
refinement, 45/80-digit agreement, independently centered eta expansions, large
Lambda comparison orders, axis regularity and sampled positivity, invalid inputs,
and independent physical-coordinate derivatives. The last reconstructs the field
at Lambda=64, eta=0.05, Y=0.4, q=0.8 with 70 digits and step 2e-5. Fourth-order
centered differences check divergence, the two leading tangential balances, and
radial pressure balance. This moderate-amplitude local check does not assert
finite-difference accuracy at every baseline point.

Appendix B's f0(z)=sum (-z/2)^n/[n!*(n+1)!] equals 0F1(;2;-z/2).
With chi=H_star^2/(H_star^2+sigma^2) and

    Z_star=-A*(1-2eta*U_star)*U_star-H_star*partial_eta(U_star)
           -d*partial_eta(Pi_axis)+4A*eta*Pi_axis,

Eq.B.13 compares Phi with f0(Y*chi), and U with U_star-Y*Z_star/(2Lambda*L).
The scan at eta=0.2, Y=1 uses Lambda=128,...,2048 and degree 18.
Observed errors decrease approximately as Lambda^-1 and Lambda^-2 respectively.
Here C stays fixed. The source's complex-neighborhood bound and its dependence
on Lambda have not been verified, so this is a pointwise numerical comparison
with the formal expressions, not a certified uniform application of Eq.B.13.

## Pressure storage and admissibility are separate issues

Some g^2-dependent increments are below the precision available when added to
the nonzero axis pressure. The solver stores Pi-Pi_axis by summing nonconstant
radial coefficients directly, and evaluates radial derivatives from coefficients.
Recovering this increment by subtracting two total pressure values can lose it
entirely. The pressure-storage CSV demonstrates that loss at 80 digits.
Decimal strings retain these small quantities in exported data; total Pi stored
as a binary64 CSV number is not a substitute for the separate increment.

Condition B.2 requires chi>0.99 wherever |Z_star| is sufficiently small.
For these test data Z_star has a zero near eta=-0.000145161 and chi there is
about 0.00927069. Thus B.2 fails for any positive neighborhood tolerance delta.
This is stronger than merely failing to check it.

The complete construction also requires the prepared outer pressure trace,
the full parameter hierarchy, a complex-neighborhood amplitude bound (B.16),
the inner exit condition (B.19), annular continuation and moment compatibility.
These have not been supplied or certified. Real-axis samples cannot establish
the required complex-domain bound.

Radial polynomial degree is unrelated to the q^(2nh) correction index in Section 5.
This stage contains no first inner correction, actual extension discrepancies,
complete V1/Pi1 reconstruction, or corrected full-PDE residual comparison.
