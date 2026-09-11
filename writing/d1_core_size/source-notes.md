# D1 source, notation, and calibration notes

This is an audit of the inputs needed for the scaling study, not an audit of the full proposed proof. Source access date: 2026-09-11. The article's computations are generated from `code/d1_core_size/parameters.json` and the reusable `nscomp` functions.

## Source record

| Source | Location used | Purpose |
| --- | --- | --- |
| [OpenAI, *Finite time blowup for Navier–Stokes*](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf) | 166-page PDF; creation metadata 2026-09-08 19:06:26 UTC | Scaling and parameter assumptions |
| Same paper | Theorem 1.1, p. 1; Section 2.1, pp. 3–4 | Forced incompressible setting and characteristic radius/speed |
| Same paper | Section 3.1, Equation (3.2), pp. 7–8; Equation (4.3), p. 25 | Similarity variables and leading-field normalization |
| Same paper | Lemma 4.8, p. 35 | Additional smallness condition on h |
| [Official repository](https://github.com/openai/NavierStokesAndEuler) | README theorem statements | Scope cross-check only; no Lean verification performed |
| [Kim et al., NIST J. Res. 110 (2005)](https://nvlpubs.nist.gov/nistpubs/jres/110/1/j110-1kim.pdf) | Equation (11), Table 3, printed p. 35 / PDF p. 5 | Air mean-free-path convention, viscosity, temperature and pressure |
| [IAPWS SR6-08(2011)](https://iapws.org/technical-guidance/release/LiquidWater.download) | Table 8, p. 11, 298.15 K column | Water density, sound speed and dynamic viscosity at 0.1 MPa |
| [NIST Chemistry WebBook: water](https://webbook.nist.gov/cgi/cbook.cgi?Name=water) | Molecular weight 18.0153 | Molar mass 0.0180153 kg/mol for spacing estimate |
| [CODATA 2022 table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt) | c, N_A, k_B, a_0 | SI constants and an atomic length ruler |
| [NASA ideal-gas sound-speed derivation](https://www.grc.nasa.gov/www/k-12/BGP/snddrv.html) | a² = gamma R T | Air sound-speed estimate; R and gamma are rounded modelling assumptions |
| [NASA-hosted low-Mach discussion](https://ntrs.nasa.gov/api/citations/19960045736/downloads/19960045736.pdf) | Printed p. 36 | Ma near 0.3 as an engineering screen, not a universal law |
| [NASA-hosted dimensionless-number discussion](https://ntrs.nasa.gov/api/citations/20000004526/downloads/20000004526.pdf) | Section 3, Equation (12) | Kn = lambda/L and conventional Kn=0.01 screening |
| [Gravelle et al., *Large permeabilities of hourglass nanopores*](https://arxiv.org/abs/1501.01476) | Abstract and study scope | Averaged permeability can agree with molecular simulations; no universal cutoff inferred |

The cached source PDFs are excluded from Git. Their URLs, hashes and metadata are retained in `code/d1_core_size/source-manifest.json`. Reproduction of the numerical study uses the recorded inputs and requires no source download.

## Definitions and derivation audit

| Symbol | Definition in this study | Limitation |
| --- | --- | --- |
| tau | Nondimensional time remaining, positive | No time unit in seconds is calibrated |
| r_c | Characteristic radial width of the intense-flow core | Not a material boundary or a measured vorticity half-width |
| U | Dominant azimuthal/axial characteristic speed | Not angular frequency; not radial inflow; not a reconstructed full-flow maximum |
| r_0, U_0 | Reference pair at one nonzero-speed stage | Illustrative physical normalization |
| h | Small similarity parameter | Subject to the full proof's parameter hierarchy |
| beta | Assumed L_grad/r_c | A proxy, not a derivative computed from the velocity field |
| d | n^(-1/3) from number density | Mean spacing, not kinetic diameter or a universal cutoff |
| a_0 | Bohr radius | Atomic reference length, not the radius of all atoms |

For the leading flow, the source defines X=r²/(2q), tau=q(1-eta²) and u_theta=q^(-A)E(X,eta), with A=1/2+h. At the midplane eta=0, q=tau. One possible consistent radial marker is fixed X=X_c>0, giving r_c=sqrt(2X_c tau). A fixed nonzero profile value E(X_c,0) then yields the stated normalized radius-speed law. This illustrates how constants cancel in a ratio; this study does not construct X_c or E numerically. Maxima and other profile diagnostics need their own normalization.

Algebra for the surrogate:

    s = tau/tau_0
    r/r_0 = s^(1/2)
    U/U_0 = s^(-1/2-h)
    s = (U/U_0)^[-1/(1/2+h)]
    r/r_0 = (U/U_0)^[-1/(1+2h)]

The code checks this result against independent time parameterization and a symbolic logarithmic-slope calculation. It also checks the inverse function, the reference anchor, monotonicity, unit changes and invalid inputs.

The paper's comparison estimates alone do not determine exact dimensional constants or certify any finite-time interval of approximation. The equalities above define an explicit surrogate with a fixed anchor.

## The h ceiling and the limits of the sweep

Lemma 4.8 states T_d=exp(M_d)+10 and 0<h<min(1/100,lambda,exp(-T_d)). Therefore T_d>10 and h<exp(-10), without knowing M_d. This is a derived necessary ceiling. It is not a sufficient condition or an explicit construction parameter. Other smallness choices in the proof remain unresolved here.

For q=U/U_0 >= 1, the exponent-only multiplier relative to the inverse law is

    r_h / r_inverse = q^[2h/(1+2h)].

Using exp(-10) as an endpoint bounds this multiplier from above within the anchored power-law family. The endpoint is excluded by the strict inequality. At q<1 the envelope ordering reverses; the implementation handles both directions. Unknown calibration factors are not covered by this envelope.

Main tables use h=0 as the h->0+ limiting surrogate. The broad sweep h=0.001, 0.005, 0.009 and 0.01 is deliberately illustrative: these values exceed the necessary ceiling. The h=0.01 curve represents the introductory endpoint, not a valid full-profile instance. No numerical plot is claimed to realize the proof's parameter hierarchy.

## Calibration ledger

| Item | Choice or calculation | Status |
| --- | --- | --- |
| Reference radius | 1 mm | Illustrative; no laboratory preparation specified |
| Reference speed | 1 m/s | Same reference stage as radius; not initial rest |
| Radius versus diameter | All main comparisons use radius; CSV also exports diameter | Explicit |
| h | 0 for main estimates; necessary ceiling and broad sweep exported | Limit / sensitivity only |
| Fluid state | Air: 296.15 K, 101.3 kPa; water: 298.15 K, 0.1 MPa | Frozen reference properties |
| Air density | p/(RT), R=287.05 J/(kg K) | Ideal-gas approximation |
| Air sound speed | sqrt(gamma R T), gamma=1.4 | Calorically perfect diatomic-gas approximation |
| Viscosity | nu=mu/rho | Derived from reference properties |
| Matching r_0 U_0 / nu to the dimensionless profile | Not performed | No fully calibrated air or water realization |
| Physical time-to-collapse | Not assigned | CSV gives only relative time tau/tau_0 |
| L_grad | beta r, baseline beta=1; sensitivity 0.1 to 10 | Proxy; tensor velocity gradient uncomputed |
| Mach criterion | Ma=0.3 marker | Screen, not a hard universal boundary; even less definitive for liquids |
| Gas Kn criterion | Kn=0.01 marker | Uses the stated lambda convention and gradient proxy |
| Water spacing | (M/(rho N_A))^(1/3) | Mean-spacing estimate; no wall confinement assumed |
| Liquid relaxation time | Not supplied | No quantitative relaxation threshold claimed |
| Molecular / atomic bands | 1–3 water spacings; 1–2 Bohr radii | Reference rulers, not validated failure bands |
| Correction pulse scales | Not reconstructed | Core screening alone cannot validate them |
| Thermodynamic evolution, heating, cavitation | Not evaluated | Could invalidate assumptions before the displayed screens |
| CFD spatial / temporal grid | None | No numerical mesh threshold or DNS claim |

Restoring units gives U_*=L_*/t_* and nu_hat=nu t_*/L_*². With nu_hat=1, U_*L_*=nu, so r_0 U_0/nu must match hat(r_0)hat(U_0). The two reference fluids have different values of this Reynolds-like product, reported in `data/summary.json`; they cannot both be called the same fully matched dimensionless profile at the same stage without further work.

A normalization factor F is represented by effective r_0=F r_0_nominal. It is exactly degenerate with changing the reference radius in this model. The sensitivity CSV reports both interpretations but they are not combined as independent uncertainties. Reference speed and radius may also be coupled once an actual profile and viscosity are fixed; the separate sweeps are calibration scenarios, not confidence intervals.

## Spatial and temporal continuum interpretation

For a continuum averaging volume, one needs a length ell containing many molecules but much smaller than the flow variation length. The estimate n ell³ only counts molecules in a local cubic averaging volume. It is not the number of molecules in the whole anisotropic core. Radius approaching d indicates loss of radial scale separation even if the axial region still contains many molecules.

For air, the ideal-gas mean thermal speed is estimated as sqrt(8RT/pi), and lambda divided by it gives a collision-time proxy. Multiplying that proxy by U/L_grad gives a velocity-variation diagnostic, not the actual symmetric strain rate; rigid rotation illustrates the distinction. This diagnostic is linked to Kn and U/thermal-speed and is not an independent validated cutoff. The proxy is exported in the summary but not used to classify the entire flow.

For water, no measured relaxation time for the specific observable and thermodynamic path has been selected. A dilute-gas mean-free-path formula is not used. Molecular correlations, relaxation, pressure changes, and fluctuations require a later case-specific study.

The paper's oscillatory corrections can introduce smaller wavelengths and faster variation than the visible core. D1 does not quantify those scales. The plotted two-screen comparison is therefore a limited, conditional ordering, not a declaration that the whole solution is physically valid before the first marker.

## Validation and reproducibility

See `code/d1_core_size/validation.json` for the actual build outcome and `code/requirements-lock.txt` for the captured Python environment. The article executes the same functions as the notebook and regenerates its figures. The canonical editable article is `article.qmd`; `article.html` and `article.tex` are generated outputs.

Source values are kept at their reported precision in the parameter file for traceability. Display tables are rounded; those digits should not be read as physical predictive accuracy. Every CSV column carries units in its name where applicable.
