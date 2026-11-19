# Euclid DR1-Foundation × blind H I surveys: dark-galaxy census

**Claim (provisional): no starless H I galaxy is established by this investigation until a frozen, independently inspectable source census and pixel-based upper limits pass the referee checks.** A catalog nonmatch alone is not a dark galaxy. This study will report explicit eligibility denominators and failed/uncertain coverage rather than treating archive errors or undetected emission as negative evidence.

## Question and test

Search public ASKAP WALLABY, ALFALFA alpha.100 and FAST FASHI H I detections over the **12 November 2026 Euclid DR1-Foundation** imaging footprint. For each genuinely covered H I beam, separate optical/NIR counterparts, blends, tidal debris, Galactic HVCs and RFI from candidates lacking an obvious stellar counterpart. Estimate stellar-light and stellar-mass upper limits only where calibrated VIS/NISP pixels, PSF, masks, aperture/background statistics and distance assumptions are available. Compute H I masses and dynamical constraints only for sources with secure velocities, fluxes, reliable widths, distances and geometries.

## Decision gates

1. Freeze exact release identifiers, catalogs, timestamps and query bytes/hashes; enumerate all input rows and why each is or is not in the Euclid image footprint. No circular cuts around a merely assumed footprint.
2. Use beam-aware association with Euclid, Legacy Surveys DR10 and unWISE; seek DESI/SDSS spectroscopic or justified photometric redshifts, and distinguish noncoverage from nondetection.
3. For counterpart-free sources, inspect calibrated images with noise and surface-brightness sensitivity; quantify aperture-corrected upper limits before quoting M_HI/M_* lower bounds.
4. Audit beam confusion, HI S/N and RFI, W50 vs distance/inclination, HVC velocities, neighboring galaxy redshifts and tidal structures. Compare aliases against AGC 229385, Cloud-9 and published ASKAP/FAST clouds.
5. An independently reviewed detection-level result is needed to claim a dark galaxy. Neither a pipeline output nor a promising sky location alone qualifies for `[REMARKABLE]`.

Data access and computation run in a single path-filtered push job on the study branch. Generated results are committed in batches outside the trigger paths. Results, receipts, audit and limitations will be added below after the run.

## Data credits

[Euclid DR1-Foundation at IRSA](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_dr1.html), [ALFALFA alpha.100 at Cornell](https://egg.astro.cornell.edu/alfalfa/data/), [WALLABY released data](https://wallaby-survey.org/data/), [FASHI survey tables](https://nadc.china-vo.org/res/r103421/?lang=en). This project has made use of the DR1 data from the Euclid mission of the European Space Agency (ESA).
