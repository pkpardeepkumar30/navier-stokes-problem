# D6 implementation ledger

The finite moment-repair substep, a local leading-inner solver, the
ideal-schedule axis-pressure integral, and a logarithmic axis-amplitude
representation with an analytic complex-neighborhood bound are implemented.
The leading inner source/shear/exit conditions have targeted numerical checks.
D6's full order-one PDE residual
comparison is still incomplete. This ledger identifies the missing inputs rather
than presenting the educational fields from D2-D5 as paper corrections.

| Component | Source | Current computational state |
| --- | --- | --- |
| Explicit heat exterior | Eq. (4.29), p.33 | Evaluated in D3 as a normalized family; matching constants/boundary not calibrated |
| Analytic leading-inner numerical engine | Eq. (4.13); Appendix B, especially B.5 and B.15 | Evaluated in D6 stage B with analytic validation axis data; radial degree, precision, overlap and physical-coordinate checks |
| Ideal-schedule axis pressure | Lemma A.5, Eq. (A.21) | Evaluated in D6 stage C with logarithmic contributions, analytic tails, quadrature refinement and eta coefficients supplied to the inner engine |
| Construction-specific axis choices | Lemma 4.8; Appendix B.1 | Stage C satisfies an explicit sufficient B.2 bound for its selected pressure schedule; full outer hierarchy remains uncertified |
| Global axis-amplitude bound | Eq. B.16, pp.147–148 | Stage D gives an exact-rational sufficient envelope on an explicit complex rectangle and an upward-rounded log(C) |
| Inner source, shear and exit | Eqs. B.17–B.19, pp.148–150 | Stage D passes declared samples for Lambda/P_star^2=1e24 and resolves a failed narrow band at 1e16; no uniform series remainder or between-sample bound |
| Joined annular profile and moments | Theorem 4.6; Appendix A | Not assembled; actual support radii and reserved intervals are unassigned |
| First inner correction | Eqs. (5.2)-(5.7), pp.46-48; Lemma 5.1 | Coupled system identified; not numerically solved |
| Extension cutoff and discrepancies | Eqs. (5.10)-(5.15), pp.49-51 | Source definitions recorded; no actual inner extension to integrate yet |
| Finite five-moment repair | Eq. (5.16), p.52 | Implemented with explicit representative bumps and declared test discrepancies |
| Complete V1 and Pi1 reconstruction | Eq. (5.15) | Not evaluated for a joined first coefficient |
| First corrected momentum residual | Eqs. (5.8)-(5.9) and physical powers (5.1) | Not evaluated; no residual-versus-time or correction-order claim |
| Oscillatory correction cycle | Sections 7-9 | Not implemented; D4-D5 explain selected mechanisms with separate examples |

## Why the available exterior is insufficient

D3's explicit exterior solves its heat balance exactly. It does not supply the
inner pressure datum or the joined annular shear. Choosing arbitrary Gaussian data
for those fields would define a different manufactured problem.

At order one the unknowns phi1, U1, and Pi1 are coupled. V1 follows from incompressibility.
The axial-viscosity source and Omega0 use the actual order-zero profiles. Pressure
therefore cannot be set independently after solving only a scalar heat correction.
The first coefficients have zero axis values under Lemma 5.1.

The available paper specifies a mathematical construction with choices and bounds.
This ledger does not claim the profiles are impossible to compute or that the
construction is nonconstructive. It records what this repository has actually evaluated.

## Formalization repository inventory

[Official repository](https://github.com/openai/NavierStokesAndEuler/tree/f9e8bc5b38b6e212696e8a30e3e91517af887bbd)
was inventoried at commit f9e8bc5b38b6e212696e8a30e3e91517af887bbd, on 2026-09-11.
The complete tree contained 2,669 files and no .py, .jl, .csv, .npz, or .ipynb files.
It includes Lean files for axis profiles, moment repair, and corrections.
This filename inventory is not a claim that no computable realization can be extracted.
No Lean build, theorem-dependency audit, or independent formal proof check was performed.

The cached tree record is code/.cache/sources/openai-repository-tree.json.
The source manifest records the commit independently of that ignored cache.

## Next concrete work

The [stage-B study](../d6_inner_profile/article.html) supplies a local radial
coefficient engine for the leading inner equations. Its analytic pressure test
datum has not been matched to an exterior. Its small local equation defects do
not establish the construction's admissibility or cancel the full PDE residual.

The [stage-C pressure study](../d6_axis_pressure/article.html) now evaluates the
ideal-schedule pressure and supplies its eta coefficients to a local inner
calculation. Its algebraic pressure bounds also give a sufficient B.2 check.
The [pressure specification](../d6_inner_profile/outer-pressure-next-step.md)
records why pressure-preserving later corrections were not required first.

The [stage-D study](../d6_global_axis/article.html) now supplies the global
logarithmic amplitude and a sufficient complex-neighborhood B.16 bound.
It checks the inner source/shear/exit conditions on a broad grid and two
targeted scales. It resolves a pressure-dependent parameter band whose
19 baseline coordinates all merge into one binary64 value. The larger
declared Lambda passes the sampled signs and exit threshold after degree,
precision and pressure-input checks. These samples do not establish
uniform inner control. The earlier stage-C C=1e4 handoff remains a local
experiment; stage D supplies its own globally bounded amplitude.

Next quantify the inner approximation and exit margins between the
declared samples before making any uniform claim; the
[continuation specification](../d6_global_axis/continuation-next-step.md)
separates that question from numerical profile assembly.
Outer parameter thresholds and moment/stress-cone checks remain distinct.
Only after adequate leading-profile assembly should the first inner correction
and its actual extension discrepancies support a full before/after PDE comparison.

The D6 stage-A lambda scan supplied the completed D7 conditioning study:
it separates near-degenerate moment weights, large required coefficients, and
floating-point loss. Its improved linear solve is available for later actual
discrepancy data.
