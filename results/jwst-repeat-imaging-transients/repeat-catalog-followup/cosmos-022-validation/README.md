# COSMOS-022 detector validation and archival host spectrum

At (150.122033°, +2.323386208°), the compact central signal increases in independent, unresampled JWST CAL exposures. The host is SMUVS 90474 / COSMOS2020 Classic 953994, with a pre-event VLT/MUSE spectroscopic redshift **0.5067**. The variability mechanism remains undetermined: **an unclassified nuclear brightening**, consistent with weak AGN variability or a nuclear transient. A new classification-quality spectrum during the bright state and a later NIRCam epoch in matching filters are the decisive missing observations.

## Native pixels

The fixed 0.16″ radius aperture was placed on each original detector science plane using its own WCS. Each exposure has its own local sky annulus. Fluxes are aperture measurements in nJy with good pixels only, without point-source aperture corrections or empirical gain rescaling.

| Filter | Early exposures | Early range | Early median | Later exposures | Later range | Later median |
|---|---:|---:|---:|---:|---:|---:|
| F162M | 9 | 572–666 | 611 | 8 | 710–900 | 843 |
| F182M | 9 | 595–677 | 624 | 8 | 698–880 | 824 |

Three distinct detector visits contribute to each epoch per filter. Odd- and even-dither medians increase in both bands. All 34 cores have zero SATURATED, DROPOUT, PERSISTENCE and OUTLIER flags in the CAL DQ planes. The input sample does contain DO_NOT_USE, HOT and JUMP_DET flags. After rejecting every exposure with even one DO_NOT_USE core pixel, six early and four late remain per filter. Their ranges are 572–666 versus 858–900 nJy in F162M and 603–677 versus 830–880 nJy in F182M. Removing every flagged pixel likewise retains nonoverlapping early/late flux ranges for exposures with ≥90% clean weighted aperture coverage. A persistence detection algorithm need not flag every real artifact. Independent locations on multiple detectors and the surviving disjoint subsets are the stronger checks. [CAL fluxes and DQ counts](native_cal.csv) · [Split dithers](native_splits.csv) · [34 original exposures, byte ranges and checksums](native_receipts.json)

## Pre-event host spectrum

Rosani et al. report a VLT/MUSE redshift of 0.5067 for SMUVS 90474, 0.124″ from this position, at original quality flag 3. The COSMOS Spec-Z Compilation DR1.1 independently cross-identifies that entry with COSMOS2020 Classic 953994, using corrected coordinates 0.068″ away. The compilation's 95% confidence class refers to redshift identification. The 3D-HST galaxy COSMOS 14157 lies 0.136″ away and has photometric redshift 0.5189. [Rosani et al.](https://doi.org/10.1051/0004-6361/201935782) · [Catalog row](smuvs.tsv) · [Compilation match](specz.json)

An independent extraction from the 2015 ESO MUSE cube uses a source-centered 0.7″ aperture and 1.2–2.0″ local sky annulus. Narrow [O II], Hβ, and both [O III] lines recover the published redshift. [O III] λ5007/Hβ is 2.91 ± 0.30, and [O II]/Hβ is 1.73 ± 0.22, with local-template fitting errors. [Ne V] λ3426 and He II λ4686 have absolute S/N below 1. Historical broad Hβ templates with FWHM 2,000 and 4,000 km/s are consistent with zero. Their respective conservative 3σ flux limits are 1.56 and 2.82 × 10⁻¹⁷ erg s⁻¹ cm⁻². The spectrum confirms a narrow emission-line host. Its wavelength range omits Hα and [N II] at this redshift, and these data cannot exclude a weak or obscured AGN. This spectrum predates the JWST variability. [1D spectrum](muse_spectrum.csv.gz) · [Narrow lines](muse_lines.csv) · [Broad-line test](muse_broad_hbeta.csv) · [ESO input](muse_receipt.json)

## Prior activity

At this position the nearest Chandra-COSMOS Legacy entry is 41.25″ away, nearest XMM-COSMOS entry 82.49″, and nearest VLA-COSMOS 3 GHz source 11.48″. Neither De Cicco VST optical AGN sample queried (2015 and 2022) has an entry within 100″. The nearest Milliquas entry is 82.49″ away. These individual catalogs do not provide a counterpart at the nucleus or establish novelty for the 2026 event. [Queries, sources and checksums](catalog_receipts.json)

## Reproduction

The detector reader takes the parent campaign directory containing `individual_apertures.csv` and `pixel_check.py`, followed by the chosen output directory. The spectroscopy and catalog scripts write beside themselves.

```sh
python native_confirm.py /path/to/repeat-catalog-followup /path/to/output
python specz.py
python catalog_check.py
python muse_extract.py
python emission_lines.py
python broad_hbeta.py
```
