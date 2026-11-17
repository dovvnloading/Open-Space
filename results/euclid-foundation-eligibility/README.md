# Euclid DR1-Foundation: bounded visit eligibility check

The IRSA image service returns **21 VIS science mosaics** intersecting a 10-arcsecond circle at ICRS RA 273.1670 degrees, Dec 64.3687 degrees in the Euclid Deep Field North. They have 21 distinct archive observation IDs across two overlapping mosaic tiles. None has both time bounds populated in the returned SIA metadata. This does not demonstrate 21 independent visits: overlapping tiles and coadds can share exposures. The number of eligible three-visit sets is **unknown**, not zero.

This is an archive eligibility result, not an astrophysical discovery or a variability measurement. We did not run a source search. The precondition for freezing a disjoint third-visit holdout has not been established, and no sensitivity or false-alarm result is claimed. A download of pixel cutouts would not by itself establish independence.

## Selection, observations and access

The coordinate and radius were fixed from the [IRSA Euclid mosaic tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-mer-images/) before querying. The single bounded SIA query selects collection `euclid_dr1_wide_DpdMerBksMosaic`, instrument `VIS`, and all returned products within the circle. It returned status `OK` and 168 rows. Filtering the returned rows by science subtype and VIS band gives the inventory above. The complete response and exact URL are in `inputs/`; `inventory.csv` retains every science product ID, tile, reported times and source URL. Null time fields are blank; filename dates are processing dates, not observation dates.

The first three science observation IDs in lexical order were chosen without reference to their science contents for a bounded access test. HTTP byte-range reads returned status 206 for each corresponding public FITS product. Their captured primary headers contain none of `DATE-OBS`, `MJD-OBS`, `MJD-BEG` or `MJD-END`; this observation applies only to those three files, not the other mosaics or other header conventions. The [tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-mer-images/) describes wide mosaics made from individual deep-field passes and an official wide survey mosaic, so these labels alone are not sufficient exposure provenance. Time ranges or exposure ID associations from an authoritative release metadata product are needed before assessing the holdout.

The TAP schema response enumerates the release tables and column definitions for observation, plane and artifact metadata. It confirms the existence of time-bound columns, not that they are populated for this field. Exact ADQL, URLs and individual HTTP range receipts are in `inputs/requests.json`; SIA and TAP responses and the three 5760-byte primary headers are deposited verbatim. The rows are a footprint-intersection query, so a later pixel analysis must also establish coverage at each selected source and inspect masks.

Next, obtain official mosaic-to-exposure associations and observation times, prove exposure sets do not overlap, and then freeze a development/control split and a later independent holdout. Retrieve calibrated cutouts, masks and noise products only after that test is possible. If a service or product remains unavailable, preserve its failed response and count it as missing coverage.

## Offline reproduction

From this checkout, run:

```sh
python results/euclid-foundation-eligibility/replay.py
```

The standard-library script reads the frozen VOTable and headers and regenerates `inventory.csv` and `summary.json`. Verify the deposited bytes with `SHA256SUMS`. The raw SIA response is a service snapshot, not all Euclid coverage, and the outputs can change if the live archive is queried again. The archive's [DR1-Foundation supplement](https://euclid.esac.esa.int/dr/dr1/expsup/) describes the release and its multiple deep-field passes.
