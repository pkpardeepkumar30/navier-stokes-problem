# Continuing D6 after the amplitude and exit study

Stage D supplies an analytic B.16 amplitude envelope and finite numerical
checks of B.17–B.19. The larger declared radial scale passes those checks.
The next work should connect this numerical inner profile to the annular
construction while preserving the distinction between evidence and a
uniform certificate.

## Inner approximation and exit margins

Partition eta into the broad region, the sigma/H_star' transition, and the
much finer pressure-dependent band. Retain local multiprecision coordinates
in the fine regions. A uniform binary64 eta grid cannot resolve the latter.
Keep the analytic amplitude envelope; do not replace it with sampled maxima.

Check radial convergence and pressure sensitivity on these patches, rather
than inferring their behavior from the single convergence plot. If claiming
an inequality on every eta in a patch, also control the eta interpolation
and radial remainder errors using an analytic estimate or validated
arithmetic. A denser positive sample set alone remains numerical evidence.

This is not a requirement to formalize the entire paper before continuing
numerical work. Profile continuation can be explored with explicit finite
resolution/error diagnostics; any uniform mathematical claim must wait for
the corresponding bound.

## Annular continuation and actual construction data

Revisit the Appendix A radial schedule and Theorem 4.6 with the stage-D
inner exit data. Assign the matching constants and reserved support
intervals explicitly. The stage-C ideal azimuthal schedule supplies pressure,
but does not supply the actual axial pulses, their moment repair, or the
stress-cone checks.

Implement a bounded annular continuation experiment first, retaining both
velocity components and pressure/incompressibility conventions. Its useful
outputs are profile overlap/matching errors, the actual uncorrected moment
discrepancies, and stress diagnostics. Label sufficiently-large/small
source thresholds that remain unquantified rather than treating the finite
parameter selection as a full admissibility certificate.

## First correction and residual comparison

Once actual leading inner/annular profiles are available, use their
derivatives in the coupled first-correction system, Equations 5.2–5.7.
Compute the real extension discrepancies before invoking the existing
five-moment repair solver. Then reconstruct V1 and Pi1 and evaluate the full
momentum residual against precision and spatial-resolution error.

A local inner-correction experiment can be useful earlier if clearly
labelled. It cannot substitute for the joined field required by D6's
before/after full-PDE comparison.
