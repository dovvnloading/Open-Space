# 3I/ATLAS independent fit and outgassing recoil: computed data

**Claim scope:** this is a physical sufficiency grid, **not** an identification of a jet direction or a proof of anomalous forces. All numbers are reproduced from the frozen input rules and public sources; the independent fit omits JPL's peak-delay DT and astrometric debiasing.

Retained full-arc sky positions: **8301** from 8345 supported; independent weighted RMS **0.637 arcsec per coordinate**.

The independently fitted symmetric inverse-square RTN A coefficients, in 10^-8 au/day^2: +1.8563 +/- 0.3147, +3.8648 +/- 0.2404, -0.7445 +/- 0.0228.

Fitted unit RTN direction: [0.42657941712493513, 0.8881213949531609, -0.17108591032519888]. The nongravitational acceleration at heliocentric distance r is fitted as A/r^2.

JPL #55: RTN coefficients [4.618075728659534, -0.9559184898971008, -0.8294590264785864] with DT=28.46530795960495 days; the vectors differ by 74.7 degrees. JPL's 117 positions (2025 Oct 3–2026 Oct 16) are **not** our MPC sample (2025 May–2026 Apr), so this angular difference cannot be interpreted as a physical conflict.

### Flux-to-recoil scenario (density 500 kg/m^3; 500 m/s; fully collimated efficiency eta <= 1)

| Gas measurement epoch | Radius (km) | Required eta | Magnitude possible at eta <= 1? |
| --- | ---: | ---: | --- |
| jwst_aug06_nuclear | 0.42 | 0.305 | Yes |
| jwst_aug06_nuclear | 1.30 | 9.04 | No |
| jwst_aug06 | 0.42 | 0.166 | Yes |
| jwst_aug06 | 1.30 | 4.92 | No |
| juice_nov02 | 0.42 | 0.0487 | Yes |
| juice_nov02 | 1.30 | 1.44 | No |
| jwst_dec23 | 0.42 | 0.138 | Yes |
| jwst_dec23 | 1.30 | 4.1 | No |

The JWST near-nucleus and terminal-coma Aug 6 inputs are **alternative** extraction scenarios, not additive gas sources; both use the same observation date.

The numerical grid also covers 0.22–2.8 km radii, 300–1000 kg/m^3 densities and 200–1000 m/s exhaust speeds. A 'No' in one row is **not** proof of excess acceleration; missing gas species, uncertain active fraction, a smaller nucleus or alternative gas speeds remain possible.

Scalar measured Q cannot establish the radial, transverse and normal momentum vector: observed gas jets and nucleus spin/orientation would be needed. Large-aperture OH can arise from distributed icy grains. A fit with fixed r^-2 cannot resolve perihelion asymmetry, and fitting uncertainty is not exclusively statistical.

Detailed values: [orbit.json](orbit.json), [jpl_comparison.json](jpl_comparison.json), [recoil_epochs.json](recoil_epochs.json), [recoil_grid.json](recoil_grid.json) and [retrieval receipts](retrieval_receipts.json).
