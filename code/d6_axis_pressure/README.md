# D6 stage C: axis pressure from the ideal azimuthal schedule

Run from the repository root:

```powershell
& code/.venv/Scripts/python.exe -X utf8 code/build.py --study d6c
```

Use --compute-only to skip HTML/LaTeX. To regenerate data and figures alone,
run the same interpreter with -m nscomp.d6c.

The reusable schedule and pressure routines are in nscomp.pressure. Its
PressureTrace.eta_coefficients method supplies the optional axis_pressure input
to nscomp.inner. The baseline uses 128-node binary64 quadrature; the inner
recurrence uses 70-digit arithmetic without claiming that precision for the
pressure integral.

Seven CSV datasets record axis pressure and derivatives, interval contributions,
the radial schedule, quadrature refinement, the exact exterior tail, radial
handoff refinement, and propagated quadrature changes. Tiny contributions have
separate logarithmic and decimal-string representations; do not infer accurate
extra digits from the printed string length.

The [article](../../writing/d6_axis_pressure/article.html),
[source notes](../../writing/d6_axis_pressure/source-notes.md), parameters,
executed notebook, source manifest, and validation records accompany the study.
The selected constants satisfy the displayed scalar hierarchy and an explicit
sufficient B.2 bound. The complete outer hierarchy, moment/stress construction,
global inner amplitude/exit conditions, joining and first correction remain
unverified or unimplemented as listed in the D6 ledger.
