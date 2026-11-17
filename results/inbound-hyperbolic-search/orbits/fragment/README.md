# Pre-split C/1996 J1 (Evans–Drinkwater)

The intact comet's astrometry is consistent with parabolic entry. At inbound 250 au its gravity-only barycentric reciprocal semimajor axis is −109.4 ± 224.2 × 10⁻⁶ au⁻¹, 0.49σ from zero. Removing the 35 measurements at the most prolific observatory leaves −11.9 ± 230.0 × 10⁻⁶ au⁻¹. The later post-break fragment A solution therefore cannot establish an unbound original orbit for the parent.

## Completed

Sekanina inferred fragmentation approximately 70 days before the December 1996 perihelion, with the two nuclei first resolved in May 1997. The input filter ended on 1996 October 1. The 131 retrieved positions themselves were observed from May 10.44268 to July 18.40006 UTC, substantially earlier than this cutoff. All 131 carry the unfragmented packed designation `J96J010`. None carry a fragment A or B label. Seventy-four match the originally published discovery circular by UTC stamp, sky coordinates, and observatory. The 131 positions represent 22 stations and 44 station-nights over 68.96 days. [Sekanina 1998](https://ntrs.nasa.gov/citations/20060034802) · [MPEC 1996-K02](https://www.minorplanetcenter.net/mpec/J96/J96K02.html)

| Force model | Retained of 131 | RMS per coordinate (arcsec) | Inbound 1/a ± formal σ (10⁻⁶ au⁻¹) |
|---|---:|---:|---:|
| Gravity only | 130 | 0.71 | -109.4 ± 224.2 |
| Signed radial inverse-square | 130 | 0.71 | 2226.0 ± 8886.0 |
| Outward Marsden radial + signed transverse | 130 | 0.70 | -180.4 ± 642.0 |

The initializer used the frozen fragment B osculating conic only as a starting state for six free coordinates. Each solution fitted the pre-split positions. All three nominal trajectories reached 250 au inbound. The radial inverse-square model changed the sign, and each uncertainty interval contains zero. The gravity result also remained consistent with zero when the first and last observing night were withheld in turn.

## Evidence and limits

The exact selected MPC80 lines, input hashes, archive and ephemeris receipts, MPC crossmatches, fitted states, covariances, holdout predictions and boundary events are in [results.json](evidence/results.json). The selected astrometry SHA-256 is `4150fd718b2d7985f5d8977e6e110f264718306fcec615b49314a99f208bb00d`. The [replay](replay.py) also checked five contemporary monthly MPC PDFs. The original designation covers the entire sample independently of those PDFs. [Workflow run](https://github.com/dovvnloading/Open-Space/actions/runs/44036049021)

The intact-comet measurements cover only the inbound observing window. A fragmentation-impulse inversion and joint dynamical fit of the nuclei were outside this test.
