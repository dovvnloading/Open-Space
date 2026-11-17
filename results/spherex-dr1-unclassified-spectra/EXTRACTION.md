# Sample, measurement and spatial checks

**The northern ecliptic pole cap supplied two detector bands. The supplemental cap supplied all six across the sample.** The 54 and 59 extracted targets give 1,561 and 5,269 strict finite aperture measurements, respectively. The spectral coverage of a particular source can be smaller. The starting caps, AllWISE catalog magnitude window and tile edge cut precede residual ranking. The second tile was selected after establishing the first tile's limited wavelength coverage.

| Selection/measurement | Northern cap | Supplemental cap |
|---|---:|---:|
| ICRS center (deg) | 270.0, +66.56 | 209.53125, +30.0 |
| NESTED order-5 tile | 4010 | 2582 |
| AllWISE in 0.50° cap, 7 ≤ W1 < 12 | 148 | 59 |
| Outside designated tile | 84 | 0 |
| Less than 16 pixels from tile edge | 10 | 0 |
| Extracted | 54 | 59 |
| Strict channels per aperture, median (range) | 31 (17–34) | 98 (39–100) |
| Strict rows, median center NHITS | 57 | 9 |
| Unique Gaia DR3 match within 2 arcsec | 50 | 54 |
| Gaia GSP-Phot effective temperature available | 32 | 42 |

The supplemental center is the WCS position of pixel (255.5, 255.5), zero-based, in the order-5 tile containing (210°, +30°). Sources enter by catalog coordinates and brightness. Membership does not require successful identification, a good model fit, or a particular color. Each row of `sources.csv` retains its selection outcome. Catalog-depth and surface-brightness limits, epoch differences and WISE blending affect the underlying population.

## Input and extraction

The IRSA TAP AllWISE query freezes the full parent catalog in `inputs/allwise.csv`. The 0.51° ancillary AllWISE (W1 < 16) and Gaia (G < 19) queries cover the entire 19 arcsec source aperture around each target. Neighbors as far as 60 arcsec can fall outside those ancillary caps near the sample boundary. `inputs/bright.csv` extends to 0.75° for W1 < 9. `inputs/phot.csv` gives AllWISE-associated 2MASS photometry for the original selection. A unique Gaia source within 2 arcsec at the approximate AllWISE epoch 2010.5 supplies the 2025.9 aperture position, using linear catalog proper motion. Gaia sources lacking proper motion are held fixed. Other apertures use the WISE position. 2MASS and Gaia identifications can themselves be confused.

`acquire.py` uses the repository's nested HEALPix and tile URL functions and requests the full IMAGE, FLAGS and NHITS FITS object in 8 MiB byte ranges. Each response must return the requested range, size, ETag and HTTP 206. The assembled MD5 must equal the single-part ETag. The complete SHA-256 and range hashes are in `receipts/cube.json`. The small raw query responses and query text, provenance, byte counts, SHA-256 and UTC times are under `inputs/` and `receipts/`. The supplement keeps its own files under `wide/`. Both cubes were acquired on 2026-11-17 from the DR1 release. The immutable byte checks for both cube objects and all ten input catalog responses are in `input_hashes.csv`. [DR1 all-sky cube, DOI 10.26131/IRSA679](https://doi.org/10.26131/IRSA679)

`extract.py` calls the existing `read_cubes.unpack` and `photometry.extract` implementations in `results/spherex-waste-heat/spectra/`. Those helpers validate pixel/sky transformations, use the FITS channel grid and aggregate 6×6 fractional sky pixels over a 19 arcsec circular aperture. The 36–60 arcsec annulus supplies its local median background. The aperture retains at least 97% unmasked area and at least twelve background pixels. Flux densities are aperture mJy, without a wavelength-dependent aperture correction. The error proxy is 1.4826 × MAD across the available offset sky apertures (24 positions, at least eight valid, using the more permissive mask). Diffuse backgrounds, confusion and correlated sky structure affect it. It omits point-source photometric and time-dependent systematics.

Strict pixels have positive NHITS, finite intensity, and no bits other than SOURCE (21). `flux_allow_outlier` also permits OUTLIER (19). The cube rejects several bad input-pixel classes during coaddition. It retains the bitwise OR of OVERFLOW, SOURCE and OUTLIER from contributing pixels. Rejected contributions can distort a point source even where the output flags show only SOURCE. NHITS counts input pixels and cannot establish the number or independence of observing visits. The table stores central-pixel NHITS and flags separately from aperture coverage. Channels 51 and 52 were intentionally blanked by the mission. The two northern detector ranges are 1.113–1.631 and 3.833–4.405 µm. [SPHEREx Explanatory Supplement, §§2.5, 3.2](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp.pdf)

## Coadd tests

`scene.py` repeats each ranked source's winning channel or three-channel feature using 13, 19 and 25 arcsec radii and two annulus choices under both masks. At width one, the straight-line wavelength continuum uses offsets 2–4 channels on each side. At width three it uses offsets 4–6 channels on each side. Every window is constrained to one detector. At least two valid fluxes on each side and all feature bins are required. Because this test fits an independent local straight line, its contrast can differ from the template residual. The peak and sideband intensity centroids are computed on the sky in a 25 arcsec core after annular sky subtraction, clipping negative intensity to zero for centroid weights. These are positions of the detected light, with no PSF fit and no formal error. Sparse or spatially extended backgrounds can displace them.

`coadd_checks.csv` stores these contrasts and centroids, center-to-annulus NHITS, catalog neighbors, the nearest other W1 < 9 source in the bright-star query, and a broad flux-scale test. The latter divides the measured median 1.20–1.30 µm flux density by the AllWISE-matched 2MASS J flux density (1594 Jy zero point). Different passbands, epochs, saturation, missing core pixels, source extension and blends can change this ratio. Across measurable apertures its median is 0.959 in the north and 0.982 in the supplement. A ratio of 0.09 signals a large flux inconsistency; alone it cannot identify the cause. [2MASS zero points](https://irsa.ipac.caltech.edu/data/2MASS/docs/releases/allsky/doc/sec6_4a.html)

## Replay

With the repository checked out, use Python with NumPy, SciPy and Astropy. From this result directory run:

```bash
python verify_inputs.py
python acquire.py cube
SPHEREX_FIELD=wide python acquire.py cube
python verify_inputs.py --require-cubes
python extract.py
SPHEREX_FIELD=wide python extract.py
# After producing the residual rankings:
python scene.py nep_two_band/ranks.csv
SPHEREX_FIELD=wide python scene.py wide/full/ranks.csv
```

Omit the `cube` argument to requery the catalogs as well. Frozen hashes guard against silently using changed archive bytes. Each cube needs about 321 MB of temporary storage. The full FITS and blank-aperture cache can be regenerated. Classification, independent visits, and injection/recovery have their own outputs.
