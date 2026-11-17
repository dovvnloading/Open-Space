# SPHEREx DR1 and stellar waste-heat candidates

**No star in this study has a validated infrared excess requiring an artificial source.** In a 12.565 deg² high-latitude field, zero of 404 nearby stars met the W3/W4 excess criteria. A separate search on the central DR1 tile found zero qualifying 4.42–4.90 µm excesses among 55 stars with usable matched comparisons. The samples overlap. A featureless thermal continuum alone cannot distinguish circumstellar dust from waste heat.

## Published targets

JWST/MIRI imaging and spectroscopy identify the sources of the long-wavelength excess in Hephaistos D and E as galaxies projected within about 1 arcsec of the stars. The galaxy at D has redshift 0.922 and an obscured-AGN-like spectrum. E has redshift 0.410 and a dusty starburst spectrum. Both fields look substantially stellar near 5 µm. High-resolution radio imaging identifies an AGN off the star at G. Its presence makes WISE confusion a strong explanation, although radio imaging does not decompose the infrared flux. B and C have strong evidence of contamination. A has a suggestive radio association and F has marginal evidence. The weaker additional targets H and J have marginal evidence of contamination, while I lacks a convincing interloper in the published archival inspection. None of the latter cases establishes an intrinsic stellar excess. [Original list](https://doi.org/10.1093/mnras/stae1186) · [Optical follow-up and H–J](https://arxiv.org/abs/2607.25701) · [Archival astrometry and imaging](https://arxiv.org/abs/2607.03619) · [Radio interferometry at G](https://academic.oup.com/mnrasl/article/538/1/L56/7954734) · [Resolved JWST spectroscopy](https://arxiv.org/abs/2607.09460)

We extracted DR1 aperture spectra for the ten Hephaistos positions, five Contardo–Hogg objects, their Gaia-matched comparisons and two pilot field stars. Individual channels have unequal coverage and reliability. For candidate E the strict aperture is 0.43 mJy above a local continuum at observed 3.293 µm, in a single channel. Its background galaxy's 3.29 µm PAH band is expected near **4.64 µm**. There the blended aperture contrast is +0.30 mJy with local blank scatter 0.15 mJy. Including outlier-flagged contributing pixels gives +0.25 mJy with 0.30 mJy scatter. This does not identify PAH. Stellar CO, stellar water absorption, isolated bad channels, neighboring light and diffuse background all affect the other feature tests. [Channel spectra](spectra/channel_spectra.csv) · [Feature and blank-aperture measurements](spectra/feature_tests.csv) · [Continuum measurements](spectra/continuum.csv)

The reported central temperatures for the original seven are 100–275 K. A 275 K blackbody emits 0.678% of its bolometric luminosity between 0.75 and 5 µm. For A's published 138 K model, the fraction is 0.000142%. These are model leverage calculations. They explain why even a well-measured photosphere at 3–5 µm is compatible with a large W3/W4 excess from another source. [Temperature inputs and integral](wise/thermal_leverage.py) · [Results including temperature uncertainties](wise/thermal_leverage.csv)

## WISE photometry and final-year variability

For A–J the median absolute difference between CatWISE and corrected unWISE magnitudes was 0.015 mag in W1 and 0.014 mag in W2. These are separate reductions of overlapping exposures. CatWISE resolves a second catalog source 4.33 arcsec from B's propagated stellar position. Even a close W1/W2 match cannot localize the emission in W3/W4. [Catalog photometry and neighborhood counts](wise/catalog_photometry.csv) · [Raw inputs](wise/archive_inputs.zip) · [40 query receipts](wise/queries.json)

The final NEOWISE release supplied 173 single-exposure rows inside the ten 12 arcsec neighborhoods. Only A had two visits with at least three accepted detections per band. Its changes over 161 days were **+0.007 ± 0.025 mag in W1** and **+0.006 ± 0.072 mag in W2**, including visit-difference calibration floors of 0.02 and 0.06 mag. The other nine had one qualifying visit each. W1/W2 monitoring does not measure variability of the longer-wavelength component. [Per-visit medians](wise/visits.csv) · [Variability tests](wise/variability.csv) · [NEOWISE calibration](https://irsa.ipac.caltech.edu/data/WISE/docs/release/NEOWISE/expsup/sec2_2.html)

## Blind controls and measured sensitivity

The preselected field is a 2° circular cap centered on ICRS (180°, +30°). There were 596 Gaia DR3 stars with parallax ≥3.333 mas, parallax S/N ≥10, G <17 and BP−RP between 0.5 and 3.5. The main-sequence proxy and Gaia quality requirements removed 115, 19 lacked an AllWISE match, five had nonunique matches, and 53 failed W1 or morphology quality. Of the remaining 404, all had usable signed W3 photometry and peer models and 324 also had usable W4. There were 205 W3 detections at ≥5σ and zero W4 detections at ≥5σ. For each band the trigger required a ≥5σ measurement, a ≥5σ excess above the matched photosphere, and ≥0.25 mag of excess. [Frozen criteria](blind/selection.md) · [Every Gaia disposition](blind/screen.csv)

| Injected blackbody | Recovered at 1% of stellar luminosity | Recovered at 10% |
| --- | ---: | ---: |
| 100 K | 84/334 (25.1%) | 232/334 (69.5%) |
| 300 K | 202/334 (60.5%) | 334/334 (100%) |
| 600 K | 157/334 (47.0%) | 330/334 (98.8%) |

These injections add one physically normalized blackbody to each real control's measured fluxes. The star's W1 flux fixes the common normalization across bands. A preliminary calculation normalized W3 and W4 separately to their empirical peer medians. The common-W1 calculation above supersedes it. Catalog inclusion, image blending, noise realizations and quality flags are held fixed. The experiment measures recovery after those cuts. The field and the injection sample cannot support a Galaxy-wide occurrence constraint. [All temperatures and luminosities](blind/injections.csv) · [Control identities](blind/injection_controls.csv)

The fixed central-tile cohort included 58 field stars with W1 <12. Fifty-six yielded sufficient photometry and local sky. Fifty-five had at least four comparison stars. Zero passed the 25% and 5σ red-continuum excess thresholds, either with the strict mask or when including pipeline outliers. The observable compares 4.42–4.90 µm against 1.51–1.72 µm, so it is sensitive to warmer excesses and to stellar molecular structure. [Blind spectral protocol and limitations](blind/README.md)

## Independent published list

We included the 53 FGK stars of Contardo & Hogg. CatWISE corroborates a W1 and W2 excess of more than 0.10 mag against the authors' respective photospheric predictions in all 53, with median excesses of 0.208 and 0.346 mag. CatWISE and unWISE share input exposures. Only two have inverse-parallax distances within 300 pc. Both have published natural counterparts, the optically dipping debris-disk system TIC 43488669 and extreme debris-disk system TYC 5940-1510-1. [FGK sample](https://doi.org/10.3847/1538-3881/ad6b90) · [Dipper](https://doi.org/10.3847/1538-4365/abbc17) · [Extreme debris disk](https://doi.org/10.3847/1538-4357/abdc26)

One red aperture provides a check against a real astrophysical excess. Gaia DR3 **3183786280737192704**, the dipper TIC 43488669, has a strict DR1 red-to-H flux-density ratio of **0.289 ± 0.005**, compared with 0.189 ± 0.005 and 0.162 ± 0.006 for two clean comparison stars. These uncertainties describe local sky scatter. The published disk interpretation and observed dips explain why it fails the requested quiet-star selection. [Spectra, control quality, and source identities](spectra/README.md) · [Literature review and astrometric checks](literature/README.md)

## Pending validation and scope

The source spectrum of a subarcsecond interloper is inseparable from the M dwarf in a 19 arcsec DR1 aperture. Candidate-level dust or companion exclusion requires source-resolved mid-infrared photometry or spectroscopy, independent epochs and a photospheric model with quantified systematics. A low RUWE does not exclude all companions. A sparse optical time series cannot exclude permanent dimming or dips between visits. No population rate or artificial-source claim follows from this study.

The DR1 spectral all-sky cube combines surveys 1 and 2. [IRSA DR1 product and acknowledgments](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/overview_dr1.html) · [Cube DOI](https://doi.org/10.26131/IRSA679)

`wise/fetch.py` requires an explicit `--output` campaign directory; cached response/receipt pairs must match the requested URL, query, byte count, and SHA-256. Partial, corrupt, or changed-query pairs and preserved errors stop the request; use a separate new campaign for changed inputs. A zero-row successful response remains a valid receipt.

For a new D: campaign, fetch the frozen response pairs and write derived tables to a separate new directory:

```powershell
python results/spherex-waste-heat/wise/fetch.py --output D:/Open-Space/campaigns/spherex-wise-20261117/archive
python results/spherex-waste-heat/wise/measure.py --archive D:/Open-Space/campaigns/spherex-wise-20261117/archive --output D:/Open-Space/campaigns/spherex-wise-20261117/derived
```

The NEOWISE single-exposure request is restricted to MJD 60291–60524; it is not a full 2010–2024 motion baseline.

For the deposited historical `wise/archive_inputs.zip`, extract it to a new directory and explicitly permit receipts without a recorded HTTP status:

```powershell
Expand-Archive results/spherex-waste-heat/wise/archive_inputs.zip D:/Open-Space/campaigns/spherex-wise-legacy/archive
python results/spherex-waste-heat/wise/measure.py --archive D:/Open-Space/campaigns/spherex-wise-legacy/archive --output D:/Open-Space/campaigns/spherex-wise-legacy/derived --allow-legacy-status
```

`replay-provenance.json` lists the receipts with unknown HTTP status; it does not modify them. A zero-row historical response describes the preserved bytes and does not by itself establish a successful empty query or completeness. New retrievals and fetch-cache reuse require recorded HTTP 200.
