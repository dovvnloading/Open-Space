# Independent ground-based photometry of OGLE-SMC-ECL-06090

The DELVE-MC DR3 catalog separates the Gaia-identified binary from a star 1.50 arcsec away. Across two nights in different years, the target exhibits a primary dip in two bands and a secondary dip in two bands. DELVE measures calibrated point-spread-function magnitudes. The source identifiers are 0210m732.7445 for the binary and 0210m732.48525 for the neighbor. [DELVE DR3 data release](https://datalab.noirlab.edu/data/delve)

| Observing night | Phase | Band | Excess over off-eclipse magnitude | Comparison nights |
| --- | ---: | --- | ---: | ---: |
| 2019 July 10 UTC | 0.0107 | g | 0.481 mag | 10 |
| 2019 July 10 UTC | 0.0125 | r | 0.425 mag | 10 |
| 2018 December 24 UTC | 0.5091 | i | 0.374 mag | 3 |
| 2018 December 24 UTC | 0.5109 | r | 0.393 mag | 10 |

On the 2019 night, earlier points at phases 0.9678 and 0.9695 showed 0.214 and 0.191 mag dimming in g and r. On the 2018 night, later measurements at phases 0.5522 and 0.5539 showed 0.165 and 0.144 mag dimming in i and r. Band-specific baselines are medians of nightly medians at least 0.11 cycles from phases 0 and 0.5, using finite photometry with quoted error ≤0.10 mag. The script retains the full measurements and baseline definitions.

The resolved neighbor is fainter than the target by 3.565, 2.696, and 2.279 mag in g, r, and i. Even disappearance of its entire catalog-average flux from an aperture containing both stars would yield only 0.040, 0.087 and 0.126 mag, respectively. This restricts the contribution of this particular neighbor under that flux assumption. An unresolved companion, blending biases, changes in brightness, and instrumental effects remain possible.

## Timing and additional archives

- A two-eclipse fit to 45 DELVE g/r/i observations on 19 nights from 2017 through 2020 gave a separation of 0.5215 cycles for a fixed Gaussian width of 0.022 cycles. Varying the width from 0.017 to 0.035 changed the estimate from 0.5287 to 0.5046. The primary core was covered on just one night. This is a template-dependent check of the phase of the dimming. DELVE's final i exposure in 2022 fell at phase 0.071, outside the core.
- The independent SkyMapper instrument measured the secondary eclipse in four bands in a single visit on 2016 October 25 UTC, spanning phase 0.514–0.523. Approximate dimming was 0.272 g, 0.373 r, 0.359 i and 0.491–0.625 v mag. Only one off-eclipse comparison night passed the baseline cut in each band. Nineteen exposures were retrieved overall. [SkyMapper DR4](https://skymapper.anu.edu.au/data-release/dr4/)
- The NOIRLab Source Catalog has 40 DECam measurements. Thirty-five share exposure numbers with DELVE. SMASH has ten measurements, eight shared with DELVE and nine with the NOIRLab Source Catalog. Comparisons between these catalogs must account for shared photons. [NSC DR2](https://datalab.noirlab.edu/data/nsc) · [SMASH DR2](https://datalab.noirlab.edu/data/smash)
- OGLE-III supplied 48 V measurements, and OGLE-IV supplied 47. Only one OGLE-IV V point falls after HJD 2459700, at phase 0.196, giving no late minimum time. [OGLE archive](https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/)
- ESO's VISTA Magellanic Survey matched source 558360887601, 0.055 arcsec from the OGLE coordinate. The release contains 5 Y, 6 J and 19 Ks epochs. All have nonzero published quality-bit integers. They were retained without using their fluxes for a timing claim. A nearby VISTA source is 1.03 arcsec away and is 0.54 mag fainter in catalog Ks, further complicating interpretation. The 8″ EROS-2 cone and 3″ 2MASS point-source cone produced no catalog rows. The published OGLE-III identification has no OGLE-II or MACHO match. [VISTA release](https://doi.org/10.18727/archive/64)
- The separate NOIRLab file search found 27 unique v1 calibrated images on 2025 October 24 UTC near the target sky position in N501 and N540, with approximate phases 0.531–0.663 calculated from their filenames. The subsequent target-position detector-footprint query in ../decam2025 matches 22 exposures (11 per filter). The archive sets a release date of 2027 April 16. Image fluxes could not yet be measured. The archive receipts preserve their MD5 identifiers. [NOIRLab Astro Data Archive](https://astroarchive.noirlab.edu/)

## Reproduction

Run `python check.py` with NumPy, SciPy and Astropy installed. This recreates results.json from the files in inputs. Magnitudes use their original passbands. Phases use P = 1.1696719 d and OGLE epoch HJD 2460000.0369. Ground-based MJD exposure start times were shifted to midpoint with their published integration times and heliocentrically corrected from CTIO or Siding Spring using the target coordinates. OGLE times are already heliocentric. Photometric zero points and source associations can differ across surveys.

The input hashes and per-epoch measures are in results.json. Read queries are recorded in sources.json. Credit belongs to OGLE, DELVE-MC, NSC, SMASH, SkyMapper, ESO/VMC, EROS-2, 2MASS and their contributing observatories and archives.
