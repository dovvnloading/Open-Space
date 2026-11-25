# New TESS transit search — 25 November 2026

## Question and observing unit

Does a star without a previously reported planet/candidate show a repeatable,
planet-sized transit in the latest public two-minute TESS sector, with a fixed
ephemeris predicting events in other observing sectors?

The observing unit is an individual SPOC sector light curve and, during source
validation, its underlying calibrated target pixels. SAP and PDCSAP reductions
share photons. Other sectors supply distinct observations. A screening trigger
is neither a discovery nor a calibrated planet probability.

## Selection fixed before inspecting the new light curves

Public manifest checks on 2026-11-25 found Sector 108 available; the standard
Sector 109 and 110 two-minute manifest URLs returned HTTP 404. Select Sector 108
targets with TIC radius 0.1–1.4 solar radii, effective temperature 2700–6750 K,
TESS magnitude 7.5–13.5, and catalog contamination ratio at most one when known.
Missing radius, temperature or magnitude is ineligible; absent contamination
metadata stays explicitly unknown. Remove TIC duplicate/artifact entries and
hosts already in the retrieved confirmed-planet, TOI, CTOI and public candidate
tables. These catalog cuts do not replace a later positional and alias check.

Freeze the complete selection ledger, including exclusions, and order eligible
targets by SHA256 of `open-space-20261125:TIC`. Process sequential batches, retain
all outcomes, and report the final searched denominator. Additional population,
period-range or filter choices must be recorded as new exploratory stages.

## Initial discovery screen

- Keep finite, positive PDCSAP flux with quality zero. Preserve the native
  arrays and uncertainties. Bin into ten-minute intervals with at least three
  valid cadences. Split at gaps over 0.2 d. Segments need 100 binned cadences.
- Divide by a one-day running median; remove upward excursions beyond six
  robust scatter units. No lower-tail clipping is applied. Require at least
  1000 retained bins over ten days.
- Use box least squares over 0.5–13 d, trial durations 1, 1.5, 2, 3 and 4 hours,
  and frequency spacing 0.008 / baseline. Refine and retain five separate
  peaks. Run the same coarse search on inverted flux as a diagnostic.
- Preliminary triage requires white-noise depth/error >=8, depth 50–30000 ppm,
  duty cycle <=0.12, radius proxy <=22 Earth radii, three covered events with
  positive local depths, and fixed-ephemeris depth/error >=3 in each time half.
  These are screening rules and deliberately do not label a likely planet.

## Adversarial follow-up

Retrieve earlier/adjacent sectors for survivors and measure at the discovery
period, epoch and duration before fitting a combined ephemeris. Record failed
and uncovered predictions. Check local event baselines, alternative detrending,
SAP/PDCSAP on the same cadences, individual event dominance, odd/even events,
secondary eclipses, double/half periods, stellar modulation, background,
centroids and quality/boundary associations. Quantify transit-duration noise
and phase controls; their conditional ranks are not full-search false-alarm
probabilities. A full-search block-preserving null and real-background
injections are required before a strong detection claim.

For a promising survivor, use calibrated pixels and Gaia neighbors to determine
where light is lost, compare aperture choices, and inspect catalogued nearby
eclipsing binaries and period/epoch matches. Fit a physical transit model with
exposure integration and realistic stellar uncertainties. A statistical
false-positive model must state its assumptions and incomplete constraints.

## Priority and reporting

Search NASA confirmed planets, TOIs, ExoFOP's fresher TOIs/CTOIs, Kepler/K2,
SIMBAD/variable catalogs, public candidate repositories, aliases and literature.
Keep exact query receipts and retrieval dates, including failed searches. No
catalog match is only a bounded negative search, not proof of priority. Use an
internal candidate identifier; do not assign an official planet letter.

The final PR must preserve selection denominators, original product identifiers
and SHA256 digests, all rejected leads, controls, repeatable measurements,
limitations, and a discriminating next observation. A result is called a
candidate only when its actual evidence supports that status. The user's desired
outcome does not change the measured evidence or the acceptance standard.

## References

- [MAST TESS products](https://archive.stsci.edu/missions-and-data/tess/data-products)
- [NASA TOI table definitions and update cadence](https://exoplanetarchive.ipac.caltech.edu/docs/API_TOI_columns.html)
- [Giacalone et al., TRICERATOPS](https://arxiv.org/abs/2002.00691)
- [Kunimoto et al., LEO-Vetter](https://doi.org/10.3847/1538-3881/ae070a)
- [Pont, Zucker & Queloz, correlated transit noise](https://doi.org/10.1111/j.1365-2966.2006.11012.x)

Data credit: NASA, TESS, MAST/STScI and the Science Processing Operations Center;
TIC and the original authors of all crossmatched candidate catalogs.
