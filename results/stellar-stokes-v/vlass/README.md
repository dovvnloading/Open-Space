# VLASS single-epoch circular polarization

**No VLASS Stokes V detection reached 5σ.** Thirteen catalogue positions fall within the twelve available images. Eleven primary beam-separated, motion-corrected positions give zero ≥5σ measurements, as do 88 offset positions. The largest absolute target signal is 2.79σ. Local 5σ sensitivities are 0.892–1.045 mJy/beam. All are single-epoch measurements on 2020-08-05. None of the 4,862 confirmed-planet host positions lies inside these V images (4,480 hosts have measured motion). The VLASS release therefore contributes no host upper limits or orbital-phase test.

At the time-resolved astrometric position of the T dwarf CFBDS J005910.90−011401.3, V = +0.127 ± 0.209 mJy/beam (5σ sensitivity 1.045 mJy/beam). CNS5 277 has the same 103.2 mas parallax and (+884.7, +44.0) mas/yr motion, but its AllWISE-derived coordinates carry epoch 2000.0 in CNS5. The dwarf census gives epoch 2009.65. The propagated positions differ by 9.2 arcsec. Both rows are archived. The census position is the primary test. Two distinct Gaia sources also land within one synthesized beam. Both forced measurements are archived and their shared beam is counted once. [Position and pixel receipts](extraction.json) · [Frozen dwarf census](../samples/dwarfs.csv.gz)

## Coverage

The physical CADC storage inventory exposes 12 averaged VLASS2.1 Stokes V FITS images in the searched naming conventions. All 12 are included. The NRAO demonstration release combines common-beam channel images using a min/max rejection. These are 2–4 GHz continuum measurements. The V product is a two-dimensional array. Its identity comes from the survey's release documentation and the product name. The FITS brightness unit is Jy/beam. [NRAO single-epoch cube guide](https://science.nrao.edu/vlass/vlass-data/single-epoch-cube-users-guide) · [Science validation memo](https://library.nrao.edu/public/memos/vla/vlass/VLASS_018.pdf)

The linked demonstration webcache was offline on 2026-11-17. The CADC ObsCore polarization column had 156,189 I, one IQU and nine unspecified entries. NRAO TAP had 159,126 VLASS images marked I. A filename search in CADC's CAOM artifact index returned zero V-mean files. The independent storage inventory gave working FITS downloads. `products.csv` contains the full-product identifiers, sizes and publisher checksums. `inventory_receipt.json` records queries, timestamps and response hashes. The accessible footprint pertains to the archived filenames queried. Unprocessed visibilities and other surveys' cubes are outside this count. [NRAO outage and data products](https://science.nrao.edu/vlass/vlass-data/basic-data-products) · [CADC storage inventory](https://www.opencadc.org/storage-inventory/)

## Measurement

The observing unit is one catalogue position in one V image. `extract.py` moves each position with known proper motion to the FITS `DATE-OBS`, treating pmRA as μα cos δ. It projects every nearby position through the image WCS and checks the exact image boundaries. Rows with missing motion have a reference-coordinate footprint check. They are excluded from the motion-qualified count. Perspective acceleration, parallax, and orbital astrometric motion are omitted. The accompanying wideband Stokes I measurement uses the corresponding single-epoch continuum field and date, with its native beam.

`measure.py` samples the nearest native pixel at the predicted position (0.6 arcsec pixels), subtracts the 12–40 arcsec annulus median, and estimates local noise as 1.4826 times the median absolute deviation. A usable annulus contains at least 50 synthesized-beam areas. The search threshold is absolute forced V ≥5 local noise units. The tabulated 5σ sensitivity applies to steady, unresolved emission averaged over the released image. It does not measure instantaneous burst amplitude. Eight samples 25 arcsec away at 45° intervals use the same background, noise and threshold. A source repeated in different catalogue lists is kept in the ledger and counted once per field if positions are within one major-axis FWHM. Short offsets sample local image fluctuations and residual artifacts. They do not establish the probability that any putative stellar counterpart is astrophysical.

`extraction.json` identifies each original field and records observation times, PM flags, original header hashes, row ranges, row hashes and cutout pixel hashes. `cutouts/` contains native pixels in V and I. `measurements.csv` and `offsets.csv` preserve every measured position, including aliases. The RMS, flux and threshold columns are in mJy/beam. The result and input denominators are in `summary.json`.

## Replay

Run from this directory with Python, NumPy, pandas and Astropy. `measure.py` verifies the archived pixel hashes before regenerating numerical results. To repeat extraction, first regenerate the four astrometric target files using `../samples/freeze.py`. Use the regenerated sample directory for `--samples`.

```sh
python measure.py
python inventory.py
python extract.py --samples /path/to/frozen/samples
python measure.py
```
