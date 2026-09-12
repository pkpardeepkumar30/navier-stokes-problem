# Next D6 component: global axis amplitude and the inner exit

Source review: 2026-09-12. This specification preceded
[stage D](../d6_global_axis/article.html), which now implements the logarithmic
amplitude and an analytic B.16 bound, and supplies targeted numerical exit
checks. Uniform inner/exit control remains open; this specification is
retained to show the intended scope. The relevant source is Appendix B, especially
[B.16 on p.147 and B.17–B.19 on pp.148–149](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).

Stage C supplies a source-schedule pressure datum and a sufficient B.2 bound.
Its one-point handoff retains C=1e4 only as a local numerical choice.
The axis amplitude phi_star can grow exponentially with Lambda near the
negative zero of H_star, so this value cannot be assumed suitable globally.

The next representation should accept log(C) directly and evaluate

    log(g)=Lambda*integral_0^eta zeta(w)dw-log(C),

without first converting C or phi_star to ordinary floating-point amplitudes.
The normalized Phi recurrence remains useful. Tiny g-dependent pressure terms
must retain separate logarithmic or arbitrary-precision representations.

B.16 requires C to dominate |phi_star| on a complex neighborhood, not merely
at sampled real eta values. A concrete bound must specify that neighborhood,
keep it away from zeros of H_star^2+sigma^2 and pressure singularities,
and control the real part of the complex antiderivative.
An upper bound on |zeta| can control a short vertical extension from the real
axis, but the real-axis integral also needs a bounded numerical error.
Choosing a numerical safety factor without those estimates is not a proof
of the complex-domain condition.

After that representation and its checks are in place, evaluate the inner
solution through Y=4.1 across the full eta interval. Include points near the
H_star zero and the narrow transition of chi, using local refinements scaled
to sigma. A coarse uniform eta grid is insufficient.

At Y=4, evaluate the source-defined p1 and p2 and the gap
p1+p2^2/p1-2, together with p1>0 and the source's angular-source condition.
Very large ratios may require logarithmic comparison. Refine radial degree,
eta resolution and arithmetic independently, and propagate pressure quadrature
uncertainty. Sampled positive gaps are numerical evidence; a uniform assertion
requires a between-sample bound or a separate analytic estimate.

The selected finite outer schedule still has uncertified sufficiently-large/small
thresholds. Even a successful B.16/exit study would leave those thresholds,
the outer moments and stress cone, annular matching, and the first q-power
correction as separate requirements.
