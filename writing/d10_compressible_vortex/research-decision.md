# D10 research decision

Current result: a reproducible educational comparison using a known exact
isentropic Euler vortex family. No original research-paper claim is established.

The density-deficit marker, the distinction between exterior and local Mach,
and the comparison with thermodynamic-gradient gas scales make the physical
message explicit. The radial solver and separate boundary-error experiment
provide an independently checked computational example.

The [NASA vortex study](https://ntrs.nasa.gov/api/citations/20150018403/downloads/20150018403.pdf)
already places this Gaussian isentropic family within established computational
fluid dynamics benchmarks. Our derivation, radial integration, and parameter
sweep do not support a claim of a new equilibrium, solver, or singularity result.
This is a focused source check, not an exhaustive priority search.

## Scope of completion

The selected stationary comparison and its article/notebook builds are complete.
It answers the D9 equilibrium question and provides a usable
richer-physics educational installment.

Time-dependent viscous D10 research remains conditional. A comparison involving
the actual proposed collapse needs the joined leading profile and source-consistent
corrections that remain unfinished in D6. It also needs consistent density,
temperature, equation of state, initial/boundary data, forcing, and energy inputs.
Incompressible pressure alone does not supply a thermodynamic initialization.

## Manuscript gate

D11's research manuscript is deferred because no supported original result has
been established. The educational articles and local reproducibility artifacts
can be reviewed on their own merits. Preparing an archive, assigning a DOI,
submitting a paper, or publishing the articles would be separate actions.

The next technical dependency is the actual leading inner/annular profile and
its first correction, as recorded in D6's implementation ledger. A future
computational question should be tied to those fields and a defined observable
before adding a general compressible or kinetic solver.
