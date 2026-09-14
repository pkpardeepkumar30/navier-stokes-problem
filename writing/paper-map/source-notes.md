# Source and editorial scope of the reader's map

Source reviewed 2026-09-12:
[Finite Time Blowup for Navier–Stokes, OpenAI](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
SHA-256: `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.

Theorem 1.1 supplies the stated forced, zero-initial-velocity setting, smoothness
for t<1, bounded energy and velocity growth. Sections 2–3 and Figures 5–6
provide the physical description and construction map. The contents on p.1
were checked against the section-by-section reading table.

The main reading order follows conceptual dependencies, not the D-number
history of the computations. Appendix C is deliberately introduced before the
pulse-amplitude construction because its shear preparation supplies that
construction's conditions. The complete proof uses Appendices A–C within
Section 4, rather than performing them only after Section 10.

All elementary stories, arithmetic examples and simplified diagrams in the
reader series are explanatory illustrations. Statements about what the
paper proves are descriptions of its claimed result, not an independent
certification or an assertion about its reception.

The archived scientific articles, numerical data and notebooks retain their
original narrower scopes. The physical interpretation route is supplementary
and does not enter the proof. Research novelty, full CFD reconstruction and
an independent proof audit are not current project goals.

The reader manifest and build are in `code/reader_series/` and
`code/render_readers.py`. The user requested this general-public revision;
the prior technical material is preserved as `technical-article.*`.
