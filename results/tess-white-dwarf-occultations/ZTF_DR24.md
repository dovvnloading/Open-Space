# ZTF DR24 epochs for the four timed-out cones

The IRSA cloud catalogue supplied 11,477 photometric measurements for these four target fields. After the same zero-catflag, positive-error and 2-arcsec proper-motion association cuts used for the other targets, 9,354 measurements enter the independent optical checks.

| TIC | DR24 HATS epochs | With measured single-epoch coordinates | With reference-object coordinates | Accepted after all cuts |
|---:|---:|---:|---:|---:|
| 357686299 | 6,511 | 4,777 | 1,734 | 5,961 |
| 408000326 | 1,248 | 828 | 420 | 1,094 |
| 446165080 | 2,530 | 2,530 | 0 | 2,064 |
| 55997020 | 1,188 | 528 | 660 | 235 |

For each target we queried IRSA's DR24 objects table at the Gaia 2016 position with a radius equal to 11 years of measured proper motion plus 2.5 arcsec. We selected the returned object IDs in the official AWS HATS partitions. Each partition's HEALPix directory came from the catalogue partition index. The Parquet measurements supply the magnitudes, uncertainties, quality flags and heliocentric modified Julian dates (HMJD). IRSA defines these HMJDs at **mid-exposure**, reported to five decimal places in days. [IRSA DR24 HATS documentation and schema](https://caltech-ipac.github.io/irsa-tutorials/ztf-lightcurves/) · [IRSA DR24 release notes](https://irsa.ipac.caltech.edu/data/ZTF/docs/releases/dr24/ztf_release_notes_dr24.pdf)

For 8,663 epochs, direct by-ID light-curve responses supplied the exposure-start UTC MJD, actual exposure duration, and single-epoch sky position. Joins required consistent five-decimal HMJD, magnitudes, errors and flags. For the other 2,814 epochs we matched the HATS HMJD to the unique physical exposure in IRSA's science-image metadata. We computed the expected heliocentric midpoint from the recorded UTC shutter start, half the exposure time, Palomar's location and the HATS object's direction with a 1-second matching tolerance. One May 2018 shutter on TIC 357686299 was absent from its CCD quadrant's image index. All 42 other returned quadrant products from the same field and band gave one identical exposure start and duration, within 0.18 s of the HATS HMJD after conversion.

For those 2,814 epochs the HATS coordinates are the object's reference-image centroid. We tested that measured position against the Gaia prediction at each epoch. Slow proper motion, crowding and blending can make a reference-position match less specific than a single-epoch astrometric measurement. The fallback CSV preserves the coordinate source on every row. The unchanged downstream time conversion receives a physical UTC shutter start and exposure duration in both cases.

The receipts preserve the TAP query, data-product URL, individual HTTP byte ranges and hashes, and every terminal calibration or photometry response used. To reconstruct the gzip CSVs and their receipts in a new input cache, run from the result directory with NumPy, Astropy and PyArrow installed.

```bash
python ztf_native_fallback.py vetting/catalogue_checks.json.gz ztf-downloads catalog-cache
```

[Reconstruction and byte checks](ztf_native_fallback.py) · [Target receipts and folds](vetting/catalogue_checks.json.gz) · [Matched optical epochs](vetting/optical_epochs.json.gz)
