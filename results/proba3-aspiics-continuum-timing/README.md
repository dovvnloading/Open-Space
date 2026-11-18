# Continuum timing in Proba-3/ASPIICS Fe XIV images

**Noteworthy calibration result.** The 63 flagged Fe XIV images in the September 17 orbit used wideband continuum donors from a different acquisition cycle. All 63 donors were closer in representative time than the same-cycle donors in the orbit's two unflagged images. The quality bit and the actual time difference carry distinct information. Across ten frozen image comparisons, the tested interpolation changed the sign of at most **0.0367%** of originally positive retained line pixels in any one annulus after the alignment stress test. This measures sensitivity to one alternative continuum, with no new coronal feature established.

## What the released files link

We matched the complete three-file wideband parent set in each Fe XIV FITS header to a single released total-brightness product. The values below use the donor's DATE-AVG minus the line's DATE-AVG. Negative lags mean the donor is earlier. DATE-AVG is the representative exposure epoch chosen by the pipeline. It need not equal the mean exposure time at each pixel. The data release notes define QUALITY bit 2 by sequence iteration and warn that non-contemporaneous continuum subtraction can leave artifacts. [ASPIICS v03 notes](https://www.sidc.be/proba-3/version-03-release-notes)

| Orbit | Date in 2026 | Fe XIV QUALITY | Header matches | Different-cycle donors | Lag range (s) |
| --- | --- | ---: | ---: | ---: | --- |
| 795 | September 17 | 2 | 63 | 63 | −35.137 to −34.136 |
| 795 | September 17 | 0 | 2 | 0 | −50.768 to −49.768 |
| 772 | August 29 | 0 | 65 | 0 | −50.772 to −26.968 |

All **130** Fe XIV headers in the two frozen queries matched exactly one total-brightness parent set. The September flag primarily distinguishes the two acquisition schedules in that orbit. It cannot be translated to a fixed lower bound on elapsed time across the observed sequences. The two clear September line products are too few for the primary five-epoch control, so the pixel comparison uses the August orbit. Within the held archive snapshots, the measured delays depend on how exposures were scheduled. [Header audit table](lineage.csv) · [Exact header blocks](fe-headers.zip) · [Byte-range receipts](header-receipts.json.gz)

## What the continuum substitution changed

The [frozen protocol](PROTOCOL.md) selected five bracketing image triples per orbit. A line product, the total-brightness product matching its actual donor, and the first later total-brightness product form each triple. Nominal helioprojective grids and complete donor parents match. We estimated the continuum at the line epoch by interpolating between the two bracketing images and applied the resulting signed change using the conversion of 3.7420 × 10¹⁸ photons s⁻¹ cm⁻² sr⁻¹ per MSB in the line headers. The calculation was restricted to originally positive, finite line pixels on finite positive continuum support in the 1.20–1.40 and 1.40–1.70 solar-radius annuli.

| Image group | Tested triples | Nominal sign flips across annuli | Flips under all four half-pixel shifts | Originally positive retained pixels across annuli | Highest annular robust fraction |
| --- | ---: | ---: | ---: | ---: | ---: |
| September, QUALITY=2 | 5 | 363 | 240 | 2,211,573 | 0.0367% |
| August, QUALITY=0 | 5 | 0 | 0 | 2,588,658 | 0% |

The largest annular fraction came from **95** of **259,172** available positive pixels in the outer annulus at 08:50:15.656 UTC on September 17. The largest 30° sector fraction was 0.230%. Neither group reached the prespecified 5% annular materiality threshold. The focused inner annulus had at least 99.5% geometric positive-line support in the flagged sample, while the flagged outer annulus had only 74.6–77.3%. The August outer annulus had at least 99.995% support. A negative residual or a missing pixel censored to NaN adds no evidence about its sign in this comparison.

The largest 95th percentile of the absolute perturbation relative to the donor K-corona intensity was 3.03% across the September annuli and 1.16% in August. A separate, post-analysis intensity diagnostic found that the 95th percentile of the absolute shift relative to the original positive line intensity reached 4.90% in the September sample and 1.93% in the August sample. Small sign-change fractions can coexist with several percent intensity shifts. The line-intensity diagnostic was specified after the sign results and has no detection threshold. [Measurements and sector coverage](measurements.json) · [Exploratory intensity ratios](intensity-diagnostic.json)

The published total-brightness products have already had a modeled F-corona removed. The test estimates the wideband shift from their difference on a common delivered grid. It assumes the modeled F-corona terms cancel closely over 60–91 s and uses linear time interpolation. The largest solar-radius variation within a triple was under 0.001 arcsec. Neither the small radius variation nor the half-pixel stress test gives a formal bound on F-corona model differences. Ghosts, detector bands, residual registration, different solar states and the different acquisition schedules remain. No intrinsic Fe XIV brightness, individual artifact attribution, or independent validator using another telescope is established. Reconstructing signed residuals from full Level-2 parents would additionally test the model assumption. [Instrument design](https://doi.org/10.1051/0004-6361/202557026) · [First results on coronal dynamics](https://arxiv.org/abs/2511.01679)

## Reproduce from the released products

The pixel protocol and sample were frozen at `73e12e5919b761905dc395f19bb13824e5074d50` before examining any brightness values. The primary script subsequently masked infinite saturation markers explicitly. The frozen result reproduced byte for byte. An independent implementation recomputed all 20 annuli from the same 30 files using adjacent-pixel averages for half-pixel perturbations. Integer counts agreed exactly, and fractional K-brightness quantiles agreed within numerical rounding. [Independent replay output](peer-replay.json)

The exact archived API response bytes are [orbit 795](orbit795.json.gz) and [orbit 772](orbit772.json.gz). The API URLs and full FITS manifest are in [sources.json](sources.json) and [fits-inputs.json](fits-inputs.json). Each of the 30 FITS arrays is a public L3 v03 file from the Royal Observatory of Belgium ASPIICS Science Centre. The full input set is 503,798,400 bytes and is retained outside Git on D. Pin Python 3.14.0, NumPy 2.5.5, Astropy 8.0.2 and SciPy 1.17.1. From this folder, download each hash-matching FITS in the manifest into RAW_DIR and run:

```text
python select_pairs.py . NEW_SELECTION.json
python audit_headers.py . NEW_LINEAGE.csv
python measure.py RAW_DIR fits-inputs.json NEW_MEASUREMENTS.json
python intensity.py RAW_DIR fits-inputs.json NEW_INTENSITY.json
python peer_replay.py --selection selection.json --raw RAW_DIR --manifest fits-inputs.json --output NEW_PEER.json --compare NEW_MEASUREMENTS.json
```

The scripts stop on input or lineage mismatch. [SHA-256 manifest](SHA256SUMS) · [Arithmetic claims](claims.json) · [Royal Observatory ASPIICS archive information](https://www.sidc.be/proba-3/aspiics-data)

ASPIICS observations are provided by the ESA Proba-3 mission and the Royal Observatory of Belgium's ASPIICS Science Centre.
