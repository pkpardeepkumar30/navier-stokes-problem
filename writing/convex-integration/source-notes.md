# Convex integration: interpretation and source map

The public article distinguishes three levels of explanation.

1. The zigzag is an elementary scalar example: the admissible slopes are
   +1 and -1, their weighted averages fill [-1,1], and the oscillating height
   tends toward a coarse line. The corners are not smooth, and it is not a
   Navier–Stokes solution. It illustrates relaxation, mixing, and small
   values with non-small derivatives.
2. The broader convex-integration viewpoint uses differential constraints
   and fine variations to realize a relaxed object. The primary fluid
   references are [De Lellis–Székelyhidi](https://arxiv.org/abs/math/0702079)
   and [Buckmaster–Vicol](https://arxiv.org/abs/1709.10033).
   Their weak-solution setting is explained separately from the target of
   smoothness before blow-up in the present paper.
3. The [OpenAI source](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
   has two relevant, connected mechanisms. Appendix C is a particularly
   direct realization of the fine-shear/nearby-profile idea. Section 7 uses
   squared amplitudes and averaged products to realize a stress, with
   specifically designed viscous pulse dynamics.

## Appendix C

The introduction on pp.157–158 explicitly separates dependence on radial
derivatives (shear) from dependence on profile values (cumulative moments).
The joined background has a relaxed cone condition where the stronger
admissible condition may fail.

Lemma C.1, Eq.C.1 constructs a smooth periodic loop of admissible shear
vectors whose mean is the original shear. Eq.C.2 keeps a uniform positive
margin inside the admissible conditions. The actual loop is not a
two-valued scalar zigzag, and simple rounding of zigzag corners does not
prove its existence.

Proposition C.2, Eq.C.12 evaluates periodic antiderivatives at N log X and
divides their profile effect by N. Values change by O(1/N), while radial
shear changes by order one. A reserved correction patch restores the five
moments. The source chooses a sufficiently large finite N and retains
smooth profiles and endpoint collars.

Calling this a direct illustration of the convex-integration idea is our
explanatory characterization. The paper explicitly uses the term “convex
integration” in its discussion of earlier work; we do not represent it as
the paper's label for every operation in its own proof.

## Pulses and the later correction cycle

Sections 2.2 and 3.3, Lemma 7.4, Propositions 7.5–7.6 and Lemma 7.7 cover
amplification, viscous decay, positive leading squared weights, signed
amplitude increments, and exact divergence-free construction. Products
create momentum flux; spatial variation of that flux creates its net
effect. A constant covariance alone does not create a mean force.

The broad stress-realization viewpoint does not establish that arbitrary
waves obey the viscous equations. Nor does Appendix C's finite-frequency
profile modification replace the successive correction and summation
arguments in Sections 8–9.

Source accessed 2026-09-12. Paper hash:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
