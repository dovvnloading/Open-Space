# Pulsar back-traces to Galactic remnant shells

An angular back-trace through 120 pulsars and 286 shell or composite supernova remnants recovers two published pulsar–remnant associations and two additional projected passages. Neither additional passage establishes a birthplace.

| Pulsar | Shell | Closest approach | Assessment |
|---|---|---:|---|
| J0002+6216 | CTB 1 | 48.0 kyr, 0.059 shell radii | Known runaway from CTB 1. [Discovery and tail imaging](https://arxiv.org/abs/1904.07993) |
| J1101−6101 | MSH 11-61A | 20.3 kyr, 0.238 radii | Known Lighthouse nebula association. [Chandra study](https://doi.org/10.1051/0004-6361/201322588) |
| J0108+6608 | G126.2+1.6 | 177.3 kyr, 0.310 radii | Disfavored. VLBI puts the pulsar at 2.14 ± 0.15 kpc. The proposed remnant distance from nearby B stars is 4.4 ± 0.3 kpc and its model age is 41 kyr. Older remnant estimates span a wider range. [VLBI](https://safe.nrao.edu/vlba/psrpi/release/J0108%2B6608.html) · [Remnant study](https://doi.org/10.1093/mnras/stae1307) |
| J1745−3040 | G356.3−1.5 | 168.4 kyr, 0.025 radii | Sensitive to motion errors. Only 34.2% of independent Gaussian component draws satisfy the same time and impact cuts. No matching independent remnant age or distance was established. [Shell listing](https://www.mrao.cam.ac.uk/surveys/snrs/snrs.G356.3-1.5.html) |

The exploratory population uses the ATNF v2.9.0 catalog: positive characteristic age ≤3 Myr, two nonzero-reported-error proper-motion components and nonzero total motion. Green's October 2024 catalog supplies shell centers and geometric-mean angular radii. Each center within 8° is projected onto the pulsar tangent plane. The script minimizes its angular offset along a constant proper-motion trajectory with fixed remnant center. The loose scan keeps times below 190 kyr, between 0.02 and six times characteristic age, and impact under 1.3 radii. The exterior shortlist requires a pulsar currently outside the effective radius, closest approach below 0.35 radii and total proper motion divided by the quadrature component error ≥5. The cuts were chosen after exploration.

The direction-shuffled screen holds positions, ages and motion magnitudes fixed. In 10,000 full-screen trials the mean exterior shortlist is 0.97 and 178 trials return four or more. The real shortlist includes two published associations. The control measures chance projections under randomized directions. A shell's apparent center can differ from its explosion site. The age proxy, solar and Galactic motions, correlated astrometric errors, projection, asymmetric shell expansion and uncertain distances further limit the physical interpretation.

Run `python screen.py /path/to/cache` with NumPy. The script downloads the two source catalogs, verifies how many objects pass the cuts, and rewrites the three output tables. Compare the input hashes in `screen.json` when using a later, live ATNF release. [ATNF Pulsar Catalogue](https://www.atnf.csiro.au/research/pulsar/psrcat/) · [Green catalog](https://www.mrao.cam.ac.uk/surveys/snrs/)
