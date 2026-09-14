# Finishing the construction: source locations and analogy limits

[OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
accessed 2026-09-12; SHA-256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.

| Reader explanation | Source |
| --- | --- |
| Two correction tasks | Section 5 / Proposition 5.5 versus Sections 8–9 |
| Separate auxiliary workspaces | Section 3.3; Section 6, especially Lemma 6.1 |
| Recompute after every operation | Section 3.4; Proposition 9.6; Figure 6 |
| Mean motion and compatibility | Section 8.3–8.6 |
| Sum with shrinking supports | Lemma 5.4; Sections 9.4–9.5; Proposition 9.9 |
| Derivatives of the residual | Eq.3.4; Eq.9.20; Lemma 10.2 |
| Cutoff preserving incompressibility | Section 3.5; Proposition 10.1 |
| Smooth force extension | Lemma 10.3 |
| Energy and comparison | Lemmas 10.4–10.5 |
| Whole claimed result | Theorem 1.1; Section 10.4 |

The auxiliary torus is a periodic organizing space used before physical
phase evaluation. It is not an extra physical dimension. Disjoint support
removes products of distinct pulse labels with overlapping slow supports;
same-label interactions and differentiation terms remain.

The geometric series illustrates controlled summation, not the actual
correction amplitudes or an empirically measured full-PDE convergence
sequence. In the source, the residual exponent improves with the cycle,
and cutoffs are chosen using derivative-tail estimates. Local finiteness
for positive concentration scale differs from behaviour at the singularity.

The force has nontrivial transition-region terms. Flatness near the singular
point does not assert that the complete external force is identically zero.
The final comparison is stated in the source's bounded-energy smooth class.

No additional solver or full corrected flow was computed for this overview.
