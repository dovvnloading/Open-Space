# A second observing night for the original starburst galaxy

One of the 12 galaxies in the first campaign has a second SDSS spectrum under a search for the same photometric object. It is the bundled starburst galaxy, fiber 0139. The two coadds contain disjoint exposures. A bright [O III] line lands on the same sampled wavelength in both. A forced Ca II K measurement has the expected absorption sign in both, but the earlier spectrum is too noisy for a firm calcium result. The catalog search returned only one visit for each of the original three galaxies with candidate calcium features. [Original investigation](FIRST_INVESTIGATION.md) · [Machine-readable measurements](../results/repeat-observation.json)

## Select the repeat using metadata

A SkyServer query joined the 12 original spectra to all DR17 spectra sharing the same positive `bestObjID`. It returned the 12 originals and one additional row. The query includes secondary spectra. A qualifying comparison required a different coadd MJD, a target separation of at most 2 arcseconds, a BOSS `v5_13_2` galaxy solution with zero galaxy-redshift warning, and a good plate in both FITS files. The three published calcium candidates had priority, followed by the remaining objects in ascending original fiber order. An imaging-ID search can miss repeats assigned different photometric counterparts. The result is specific to this archive and matching rule. [Frozen SQL](../data/repeat-catalog.sql) · [Complete query response](../data/repeat-catalog.csv) · [SDSS repeat guidance](https://www.sdss4.org/dr17/spectro/caveats/) · [How SDSS matches spectra to photometry](https://www.sdss4.org/dr17/algorithms/match/)

Both spectra correspond to imaging object `1237665180057927873` and FITS `THING_ID 364543359`. The recorded positions coincide at RA 214.316940°, Dec +29.275123°. The coadd MJDs differ by 228 days. We parsed every `EXPIDnnn` header card, required the count to equal `NEXP`, and compared unique science integration numbers across cameras. The two sets have no overlap. Blue and red camera records for a single exposure contribute one integration to that count. Both visits used the SDSS 2.5-m telescope and the same pipeline release. Independent photons still share calibration and reduction methods. [SDSS spectral products](https://www.sdss4.org/dr17/spectro/pipeline/) · [SDSS raw exposure numbers](https://classic.sdss.org/dr7/dm/flatFiles/sdR.php)

| Spectrum | Coadd role | Unique science integrations | Camera records | Galaxy redshift |
| --- | --- | ---: | --- | ---: |
| 6497-56329-0139 | Original campaign anchor | 7 | b1, r1 | 0.08566851 |
| 6012-56101-0585 | Earlier repeat | 6 | b2, r2 | 0.08571332 |

## Measure a prediction made before opening the repeat pixels

Both spectra use the original galaxy redshift for the predictions. Ca II K at 3934.775 Å in vacuum should appear at 4271.861 Å, while [O III] at 5008.239 Å should appear at 5437.29 Å. The calcium window, wider continuum sidebands, both mask policies and two displaced control windows were fixed before requesting the repeat pixels. After the metadata identified the starburst galaxy, we added the launch pilot's [O III] positive control before opening those pixels. [Laboratory calcium wavelength](https://physics.nist.gov/PhysRefData/Handbook/Tables/calciumtable4.htm) · [SDSS vacuum convention](https://www.sdss4.org/dr17/spectro/spectro_basics/) · [Launch oxygen prediction](DATA.md)

The fixed rest-frame core extends ±3 Å around each reference. Two sidebands lie 8–15 Å on either side. Their median wavelength and median flux define two continuum anchors for linear interpolation through the core. Define fractional depth as one minus the ratio of inverse-variance-weighted core flux to inverse-variance-weighted continuum. Positive depth denotes absorption. A measurement requires at least three usable core pixels, two on each side, and positive modeled continuum. The default uses positive inverse variance and zero AND mask. The strict variant requires zero OR mask. No estimate uses a fitted shift or a recentered line window.

| Measurement | Original, MJD 56329 | Earlier visit, MJD 56101 |
| --- | ---: | ---: |
| Ca K fractional depth, default | +22.8% | +12.4% |
| Ca K fractional depth, strict | +25.3% | +12.4% |
| Displaced depth at K − 25 Å, default | −5.8% | −5.4% |
| Displaced depth at K + 25 Å, default | +8.3% | +0.7% |
| [O III] core flux above continuum, default¹ | +12.14 | +14.31 |
| Brightest eligible pixel within ±4 observed Å of [O III] | 5437.509 Å | 5437.509 Å |

¹ Flux density in 10⁻¹⁷ erg cm⁻² s⁻¹ Å⁻¹. The two oxygen peaks fall in the same bin of a shared logarithmic wavelength grid. The strict OR mask removes every oxygen core and sideband pixel in both coadds. An OR flag indicates at least one contributing exposure flagged a pixel. The displaced windows are neighboring spectral regions that may contain real lines. They provide a local comparison with the calcium window.

## What the second night supports

The calcium deficits are larger than the displaced-window deficits under the default mask on each visit. To quantify sensitivity to quoted pixel noise, we then made 4,096 perturbations per spectrum and mask, with independent zero-mean Gaussian errors scaled by `1/sqrt(ivar)`, recomputing each median continuum and index. The 2.5th–97.5th percentile range was +12.7% to +30.3% for the original and −11.0% to +28.7% for the earlier visit. These descriptive ranges use the diagonal pixel errors, omit correlations from resampling and treat the window choice, wavelength and redshift as fixed. The repeat allows zero depth even under that simplified model.

The bright oxygen emission provides a positive check across disjoint integrations and different BOSS cameras. The calcium window shows consistent signs, with inconclusive evidence in the noisier visit. A higher-sensitivity observation or a full stellar-template fit would test calcium in this starburst. This observation supplies no replication for the three galaxies in the earlier calcium investigation.

## Reproduce

From the checkout, with Astropy installed, run the following command. The script verifies the frozen metadata digest and both original FITS SHA-256 hashes, computes the indices, and prints complete JSON. An existing cached FITS file is checked before reuse.

```sh
python -m pip install -e '.[archive]'
python -m examples.repeat_observation --cache runs/repeat-fits > runs-repeat.json
```

Use `--refresh-catalog` to rerun the original SkyServer query and compare the returned bytes with the published snapshot. The JSON includes the raw camera exposure cards, source hashes, flags and coverage counts. [Reproduction code](../examples/repeat_observation.py) · [Original FITS](https://data.sdss.org/sas/dr17/eboss/spectro/redux/v5_13_2/spectra/lite/6497/spec-6497-56329-0139.fits) · [Earlier FITS](https://data.sdss.org/sas/dr17/eboss/spectro/redux/v5_13_2/spectra/lite/6012/spec-6012-56101-0585.fits)
