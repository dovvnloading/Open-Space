# No confirmed starless H I galaxy in a Euclid DR1-Foundation overlap *inventory*

**Measured claim (19 November 2026):** A frozen scan of **191,367 public blind-survey H I catalogue rows** against the official [Euclid DR1-Foundation](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_dr1.html) Wide/Deep sky MOCs identified **1,375 catalogue-entry positions** whose H I centroids or sampled half-power beam rims fall on the released tile footprint. These are **coverage overlaps, not 1,375 starless-galaxy candidates**. The experiment has **zero confirmed dark galaxies and zero validated stellar-mass upper limits**; because the image/counterpart checks were not completed for every position, it establishes *neither a negative census nor an upper limit on the dark-galaxy frequency*.

| Frozen H I release | Catalogue rows | H I centers in Euclid MOC | Sampled beam-edge only | Total eligible positions |
| --- | ---: | ---: | ---: | ---: |
| [ALFALFA alpha.100](https://egg.astro.cornell.edu/alfalfa/data/) (2018, Cornell revised CSV) | 31,502 | 116 | 1 | **117** |
| [FASHI DR2](https://arxiv.org/html/2606.31539v1) (2026) | 156,411 | 1,222 | 36 | **1,258** |
| [WALLABY pilot DR2](https://wallaby-survey.org/data/data-pilot-survey-dr2/) (2024, 30-arcsec standard source catalogue, includes DR1) | 3,454 | 0 | 0 | **0** |
| **Sum of catalogue rows** | **191,367** | **1,338** | **37** | **1,375** |

The counts are **entries in released tables**, not deduplicated astrophysical sources: ALFALFA and FASHI can see the same gas, and a WALLABY catalogue row need not be a distinct galaxy. Both Euclid MOCs were retrieved from [IRSA's DR1 archive](https://irsa.ipac.caltech.edu/Missions/euclid.html). All denominators, per-row coordinates and selection status are in [footprint_hits.csv](output/footprint_hits.csv) and [summary.json](output/summary.json). The **1,338** positions passing a *broad mass/linewidth/velocity triage* (114 ALFALFA, 1,224 FASHI) remain ordinary unclassified H I detections; this is **not** a dark-galaxy count.

## Frozen data and reproducibility

The prospectively frozen [freeze.json](freeze.json) records the exact release source URLs or TAP endpoint and ADQL; [receipts.json](output/receipts.json) records query SHA-256, response SHA-256, byte counts, execution UTC, status, and native survey column names, including failed responses from earlier archive trials. Released Euclid Wide and Deep [MOCs](https://irsa.ipac.caltech.edu/euclid-archive-documentation/data-products/) were not approximated by a hand-drawn sky rectangle. The parser keeps every input catalogue row in the selection denominator, propagates missing measurement values and records any missing releases as fatal for a complete census.

The overlap rule is: a row is inside if the **H I centroid** falls in either published Euclid Wide/Deep MOC; otherwise, **16 directional samples at half the H I beam FWHM** test a possible edge contact. Radio beam half-power radii were frozen at 105 arcsec (ALFALFA), 87 arcsec (FASHI) and 15 arcsec (WALLABY PDR2). An edge sample is neither an exact disc–polygon intersection nor verified pixel coverage. The MOCs describe tiles rather than per-band valid depth. DR1 calibration observations outside the specified Wide/Deep MOCs are not included. The *all-rows* H I and MOC matching is complete for these three fixed release tables and this geometric rule, **not for every public H I observation or all pixels in DR1**.

The single branch-specific, path-filtered [GitHub Actions job](../../.github/workflows/euclid-hi-dark-census.yml) obtains the public data and performs all remote-dependent computation, then commits derived outputs in one batch. The runner does not re-trigger when writing the `output/` files. Replay on a network-enabled Python environment (Python 3.12, requirements installed by the workflow):

```sh
python -m unittest discover -s results/euclid-hi-dark-census -p tests.py -v
python -u results/euclid-hi-dark-census/run.py
python -u results/euclid-hi-dark-census/followup.py
python -u results/euclid-hi-dark-census/imaging.py
```

The main data products are [all footprint hits](output/footprint_hits.csv), [a 32-row prioritization table](output/priority_for_followup.csv), [cross-catalogue per-position dispositions](output/followup.json), [follow-up TAP request receipts](output/followup_receipts.json), [native Euclid image access and exploratory pixel checks](output/native_pixel_checks.json), and the [machine referee gate summary](output/referee.json). These are small reanalysis products rather than copies of the entire upstream releases.

## Why this is not a dark-galaxy discovery or a negative census

The cross-catalogue program attempts a **32-source, non-random diagnostic subset** (16 from each of the two overlapping surveys), seeking Euclid MER sources, DESI DR1 spectra, Legacy Surveys DR10 catalogue objects and unWISE detections within each H I half-power beam. It records HTTP errors, blank results and saturated `TOP` responses separately. **A missing catalogue object is not an image nondetection, and even a source at the same sky location is not a host without a concordant redshift.** SDSS spectroscopy and robust photometric redshifts were not surveyed exhaustively. This is not a 1,375-object optical-counterpart census.

The native-image audit selects only the first two prespecified unclassified positions in that subset. It attempts Euclid VIS/Y/J/H cutouts and background-aperture noise checks. Such empirical aperture statistics lack verified segmentation/bad-pixel masks, point-spread-function and extended-light injection recovery, sky-systematics assessment, full H I beam apertures, and per-galaxy stellar population/dust/distance modelling. Any displayed local surface-brightness diagnostic is **exploratory**, not a calibrated upper bound on stellar flux or stellar mass. Accordingly, the machine outputs intentionally report `Mstar_upper_limit=null` and `MHI_over_Mstar_lower_limit=null`; **no M_HI/M_* lower bounds are established**.

Some catalogue line widths and masses are galaxy-like, but **W50 is not automatically a rotation speed**. Neither H I radii nor inclinations, source cubelets, line-profile asymmetry, beam confusion and RFI, high-velocity Galactic gas, or projected and velocity-matched tidal companions have been validated. A dynamical halo mass or a starless classification based on W50 alone would be unjustified. Nearby H I neighbours within 20 arcmin and ±500 km/s are recorded as a **descriptive diagnostic**, not a tidal-debris veto. The broad screen of H I mass ≥10⁷ solar masses (where known), W50 ≥30 km/s and optical cz ≥600 km/s only prioritizes, and cannot reject all HVCs or recombination-line contaminants.

## Prior literature and novelty boundary

These cases are **comparison examples, not newly reidentified objects in our selected sample**. Full source-by-source identity and literature matching remains a required follow-up, so no novelty priority is claimed.

- [Janowiecki et al. 2015](https://doi.org/10.1088/0004-637X/801/2/96): ALFALFA AGC 229385 looked dark in shallow surveys, but deeper data showed a g-band faint stellar counterpart near 26.4 mag arcsec⁻²; AGC 229383/229384 remained unmatched at reported limits.
- [Anand et al. 2025, Cloud-9](https://doi.org/10.3847/2041-8213/ae1584): independently detected compact nonrotating gas with Hubble nondetection of stars to about 10⁴ solar masses; [2026 follow-up](https://doi.org/10.1093/mnras/stag1293) addresses its M94 environment. Cloud-9 is already known, not a novelty target.
- [O'Beirne et al. 2024](https://doi.org/10.1093/mnras/stae215): WALLABY H1032−2819 has a DR10 low-surface-brightness counterpart and a plausible group-interaction origin.
- [FAST and Dark 2026](https://doi.org/10.1093/mnras/stag732): an existing list of **70 FASHI DR1 dark candidates** with optical, position, RFI and environment vetting. [AGC 322753](https://doi.org/10.1093/mnras/stag1874) is a distinct 2026 ALFALFA/FASHI candidate with both radio detections and a published `M_HI/Mstar > 40` inference.
- The repository's separate [FASHI Galactic recombination contamination study](../fashi-recombination-lines/README.md) identifies another reason that apparently extragalactic single-dish H I rows can be misleading.

The independent critical assessment is in [REFEREE.md](REFEREE.md). **Disposition: measured and reproducible footprint eligibility only; no `[REMARKABLE]` tag or findings-log entry.** The decisive next experiment is all-position forced VIS/NISP/Legacy/unWISE imaging with injections and per-band masks, then redshift/environment/HI-cube checks for genuinely counterpart-free sources.

**Euclid credit:** This work has made use of the Data Release (DR1) data from the Euclid mission of the European Space Agency (ESA). This study also uses archival publicly released Arecibo ALFALFA, FAST FASHI and ASKAP WALLABY data.
