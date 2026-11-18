# DR11 field snapshots

Both rectangles are in the southern Legacy Surveys processing region. Their boundaries were fixed in `footprint.json` before downloading any object rows or ranking positions. Boundaries are half-open. The known systems in the footprint motivated the field selection.

`acquire.py` intersects the official DR11 brick-boundary index with each rectangle and reads every overlapping southern Tractor FITS file. Each source must be marked BRICK_PRIMARY, fall inside the rectangle, have positive finite measured g and r fluxes with finite positive magnitude errors, and have 17 ≤ g₀ < 24. All measured morphologies passing those criteria are included. No row cap is used. `brick_index_receipt.json` and `brick_receipts.csv` record the source URLs, byte lengths, unmodified SHA-256 values, row counts, retrieval times and the HTTP Last-Modified values. `snapshot.json` records the SHA-256 values and validation counts for the reduced catalogue files.

| Field | Definition |
| --- | --- |
| id, release, brickid, objid | LS_ID_DR11 and the three native integer identifiers. Read IDs as integers or strings. |
| ra, dec | ICRS coordinates in decimal degrees. |
| g, r | Extinction-corrected AB magnitudes, 22.5 − 2.5 log₁₀(FLUX / MW_TRANSMISSION). FLUX is Tractor model flux in nanomaggies, a PSF fit for point-like objects. Extended-object comparisons use the respective fitted galaxy models. |
| gerr, rerr | Linearized one-sigma magnitude uncertainty, (2.5 / ln 10) / (FLUX × √FLUX_IVAR). |
| star, type | star = 1 for the Tractor PSF morphology and 0 for every other model. type retains the native label. Compact galaxies can pass the PSF classification. |
| clean | 1 when both bands have NOBS ≥ 1, flux S/N ≥ 5, FRACFLUX < 5, FRACMASKED < 0.4 and FRACIN > 0.3, and no MASKBITS from {1,2,3,5,6,10,11,12}. |
| maskbits, fitbits | Native bit fields. The selection retains cluster and resolved-dwarf flags 13 and 17 to allow recovery of known clusters and galaxies. |
| gdepth, rdepth | Dereddened AB magnitude of the local idealized five-sigma point-source depth calculated from PSFDEPTH. These per-object estimates do not measure empty sky, blending or catalogue classification completeness. |
| nobs_g, nobs_r | Native counts of exposures. |
| ref_cat, ref_id | Native reference catalogue and ID. G3 denotes Gaia DR3. A blank reference label indicates there was no reference-catalogue association. |
| pmra, pmdec, pmra_ivar, pmdec_ivar | Native reference proper motions (mas/yr, the RA component includes cos declination) and inverse variances ((mas/yr)⁻²). Require appropriate G3 flag and positive inverse variances before interpreting. |

The gzip files use concatenated gzip members in a fixed brick order. To reproduce, install Python 3.12 with NumPy and fitsio, then run `python acquire.py --work /path/to/fresh/scratch` from the analysis directory on Linux. Temporary storage contains only per-brick compressed projections and metadata. Compare the resulting per-brick and per-catalogue SHA-256 columns. The retrieval timestamps will differ.

Sources: [DR11 description](https://www.legacysurvey.org/dr11/description/), [Tractor catalogue definition](https://www.legacysurvey.org/dr11/catalogs/), [bit masks](https://www.legacysurvey.org/dr11/bitmasks/), [release files](https://www.legacysurvey.org/dr11/files/), [NERSC release directory](https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr11/).
