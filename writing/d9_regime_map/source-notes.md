# D9 source and scale audit

Access/review date: 2026-09-12. This extends D1's physical screening calculation.
It does not evaluate the complete paper profiles, audit the proof, or integrate a flow.

## Audited inputs

The parameter file snapshots D1's air, water, SI constants, and illustrative core
anchor, recording the original input hash. Ordinary regeneration uses this snapshot,
so a later change to D1's parameters cannot silently alter D9.

| Source | Location | Use |
| --- | --- | --- |
| [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf) | Section 2.1; Lemma 4.8 p.35; p.81 after Eq.7.22 | Core powers, necessary h ceiling, carrier power |
| [NIST air measurements](https://nvlpubs.nist.gov/nistpubs/jres/110/1/j110-1kim.pdf) | Eq.11 and Table 3, printed p.35 / PDF p.5 | 296.15 K, 101.3 kPa, mean free path 67.3 nm, viscosity |
| [IAPWS SR6-08(2011)](https://iapws.org/technical-guidance/release/LiquidWater.download) | Table 8, p.11 | Water properties at 298.15 K, 0.1 MPa |
| [NIST water](https://webbook.nist.gov/cgi/cbook.cgi?Name=water) | Molecular weight | 0.0180153 kg/mol |
| [CODATA 2022](https://physics.nist.gov/cuu/Constants/Table/allascii.txt) | SI constants and Bohr radius | Molecular count and atomic ruler |
| [Ashby, NASA-CR-201482](https://ntrs.nasa.gov/api/citations/19960045736/downloads/19960045736.pdf) | Printed p.36 | Air Mach 0.3 engineering convention |
| [Kuczmarski, Dimensionless Numbers Expressed in Terms of Common CVD Process Parameters](https://ntrs.nasa.gov/api/citations/20000004526/downloads/20000004526.pdf) | Section 3, Eq.12, PDF p.7 | Kn definition and 0.01 screening convention |
| [NASA sound-speed derivation](https://www.grc.nasa.gov/www/k-12/BGP/snddrv.html) | Ideal-gas formula | Sound-speed estimate |

The OpenAI p.81, NIST p.5, and IAPWS p.11 source pages were visually inspected.
NASA PDFs were fetched directly after the web reader returned 403; local text
confirmed the cited passages. Hashes are in source-manifest.json.

The air heat-capacity ratio 1.4 and specific gas constant 287.05 J/(kg K) are
rounded ideal-gas assumptions. The derived density and sound speed are not extra
measured values from the NIST table. Mean thermal speed follows by integrating the
Maxwell speed distribution; a numerical integration independently checks it.

## Definitions and distinct types of length

- a: core-radius surrogate, anchored to 1 mm at 1 m/s in the h->0+ limit.
- L=beta*a: assumed variation length. No actual derivative or phase field is supplied.
- lambda_mfp: air mean free path in the cited convention; no dilute-gas value for water.
- d=(M/(rho*N_A))^(1/3): water cubic mean spacing; not a molecular diameter.
- epsilon*L: hypothetical local averaging-cube side. Its mean count is (epsilon*L/d)^3.
- ell_car: an anchored carrier scale retaining only the source's Q power.
- 2*pi/|grad(phi)|: a local full wavelength for a specified physical phase; unevaluated.

The carrier comparison eliminates Q/Q0 between
U/U0=(Q/Q0)^(-1/2-h), a/a0=(Q/Q0)^(1/2), and
ell/ell0=(Q/Q0)^((1+h)/2). This produces
ell/a=(ell0/a0)*(U/U0)^(-h/(1+2h)).
The source gives a scale sqrt(Q)/k up to comparison factors. No hidden 2*pi
conversion or actual phase-gradient normalization is inserted.

In using this comparison, the leading background and pulse are associated with
comparable concentration scales. This is an asymptotic power comparison, not a
trajectory following one material particle or one pulse through the entire collapse.
The code freezes all prefactors. Unknown scale-dependent envelope, phase,
localization, and correction factors are excluded from the numerical bound.

h=0 is a limit; exp(-10) is an excluded necessary ceiling from T_d=exp(M_d)+10.
h=1e-6 satisfies that necessary condition alone, not the full construction hierarchy.
The 0.0781% endpoint change concerns only the normalized bare-power ratio and the
declared formal range. It is not a bound on the shortest actual scale.

## Time and constitutive interpretation

t_coll=lambda_mfp/sqrt(8RT/pi), t_adv=L/U, t_sound=L/c_s, and t_diff=L^2/nu.
All have units seconds, but none calibrates time remaining to the proposed singularity.
chi=t_coll/t_adv=Kn*U/mean_thermal_speed is a proxy. The CSV's chi=0.1 marker is
an illustrative time-scale comparison, not a validated onset of nonequilibrium.
Re=t_diff/t_adv and Mach=t_sound/t_adv are checked independently.

U/L uses the dominant characteristic speed. It is neither the actual rate-of-strain
norm nor automatically the variation rate of a smaller-amplitude correction.
The relevant amplitudes, derivatives, forcing frequency, and thermodynamic path
would be required for a field-based assessment.

The gas and liquid branches are distinct. Water's relaxation time, equation-of-state
evolution, absolute pressure, and phase changes remain unavailable. The 1 and 10
spacing lines and 1000-particle example are geometric rulers, not certified thresholds.
The comparison does not apply dilute-gas Knudsen criteria to water.

## Data and validation

Five CSV files contain scale trajectories, air crossings, gas time scales, the carrier
power comparison, and averaging-window counts. Units appear in column names or are
dimensionless. Crossings below U0 are marked as earlier than the reference stage.
There are no random inputs, flow time steps, or CFD meshes.

Eight checks test independent Q parameterization, bracketed crossing roots, Maxwell
integration, time ratios, the ordering reversal, molecule counts and unit invariance,
the inverse-law limit, and invalid physical inputs. The plotted curves remain formal
extrapolations after any screen is crossed.
