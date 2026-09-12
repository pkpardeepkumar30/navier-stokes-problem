# D6 implementation ledger

The finite moment-repair substep is implemented. D6's full order-one PDE residual
comparison is still incomplete. This ledger identifies the missing inputs rather
than presenting the educational fields from D2-D5 as paper corrections.

| Component | Source | Current computational state |
| --- | --- | --- |
| Explicit heat exterior | Eq. (4.29), p.33 | Evaluated in D3 as a normalized family; matching constants/boundary not calibrated |
| Axis data and analytic inner leading profile | Appendix B, pp.144-145; Eq. (4.13) | Not evaluated; needs the outer pressure trace and choices including h, C, Lambda, j0 |
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

First prepare a numerical representation of a single leading inner profile, using
source-consistent outer pressure data and an explicit parameter ledger. Validate
regularity, incompressibility, and leading balances before joining the annulus.
Only then generate d_U,d_E from an actual first inner correction.

The D6 stage-A lambda scan also supplies a well-defined conditioning problem for D7:
separate near-degenerate moment weights, large required coefficients, and floating-point
loss. Improvements to this small linear solve can proceed while profile assembly is developed.
