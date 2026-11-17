# WD 1856 b and two binary controls

The blind box search ranked the known clock first in 40 of 41 SPOC products for WD 1856+534 b, 8 of 8 for RR Cae and 7 of 7 for TIC 60040774. All 56 products were measured. The planet was missed in Sector 26 at 120 s. Sectors and cadences of the same star share the astrophysical signal. These counts cannot estimate discovery completeness across white dwarfs.

| System | Latest sector with 120-s product | Published period (d) | Blind 120-s period (d) | Blind 200-s result in that sector |
|---|---:|---:|---:|---|
| WD 1856+534 b | 86 | 1.407939217 | 1.40796224 | Missed in the 1.5-pixel-radius aperture |
| RR Cae, WD + M dwarf | 97 | 0.3037036366 | 0.30370800 | Rank 1 at 0.30370032 d |
| TIC 60040774, WD + M6.5 | 108 | 0.404807 | 0.40482473 | Rank 1 at 0.40481591 d |

The 20-s Sector 86 data recovered the planet at 1.40794276 d. At the independently searched 120-s ephemeris, 17 supported events in the 200-s FFI had a median aperture flux deficit of 0.69%. A Gaussian model of Gaia neighbours assigned 1.29% of the aperture light to the WD at width 0.85 pixel, varying from 0.99% to 1.75% over widths 0.65–1.05 pixels. Near this degree of blending, the nominal dilution-corrected depth is highly model dependent. The difference-image fit lies 0.22 pixel from the propagated WD position, with event-bootstrap scatter of 0.93 and 0.47 pixel in x and y. The FFI positions alone do not establish the source. Both binary controls gave stronger difference images, with offsets of 0.068 pixel (RR Cae) and 0.025 pixel (TIC 60040774). These are free-Gaussian fits and the resampling scatter omits PRF-model and astrometric systematics.

Each control was chosen from the literature and analyzed outside the blind survey denominator. The inventory enumerates the MAST public S3 target directories for sectors 1–110 on 2026-11-17 and records listing hashes. The search used all the 120-s and fast light-curve files found there, with zero quality flags, finite measurements, and uncorrected simple-aperture flux. The search had no ephemeris input. A returned period had to drift by at most one fitted full duration over the visit and agree in phase at the visit midpoint within one duration. For TIC 60040774 after Sector 33 we used only the drift test because the published period uncertainty loses phase over the intervening years. Its coarse historical ephemeris predicts an eclipse about 22 minutes from the Sector 108 measurement. The primary period checks, sources and matches are in the script and results. A raw file hash and the eight-column digest used by the repository's earlier transit survey accompany every light-curve measurement.

At the latest available 120-s control sector we also extracted a 9 × 9 TESScut stamp, projected the Gaia DR3 coordinates to the observation epoch, subtracted the per-cadence edge median and summed pixels within 1.5 pixels of the WD. The 200-s detector ran independently. For each predicted event from the 120-s ephemeris, the difference image subtracts the eclipse mean from the mean of symmetric flanks 1–2.5 durations from center. At least one quality-zero sample in the center and on each flank is required. A single circular Gaussian with a constant sky is fit to the signed mean difference image. The median aperture depth and photocenter shifts summarize individual events. SPOC and FFI exposures overlap in time and are not independent observing epochs.

[WD 1856 ephemeris](https://doi.org/10.1093/mnras/stad766) · [RR Cae ephemeris](https://academic.oup.com/mnras/article/376/2/919/1077105) · [TIC 60040774 binary](https://doi.org/10.1093/mnras/stac2197) · [Prior polluted-WD TESS transit survey](https://arxiv.org/abs/2407.21743) · [2026 known-debris census](https://arxiv.org/abs/2609.08923)

## Reproduction

Run from this controls directory, with its parent on PYTHONPATH for the shared reduction and search modules. The existing repository helper paths refer to adjacent results directories. Cache space for the products is roughly 0.5 GB.

1. Place '..' on PYTHONPATH.
2. Run 'python query.py' to refresh the holdings (optional).
3. Run 'python positive_controls.py control_inventory.json cache control_results.json --previous-fetch ../../tess-additional-transits/fetch.py'.
4. Run 'python control_pixels.py control_results.json . cache control_pixel_results.json --pixel-tools ../../tess-asymmetric-dips/pixels.py'.

The supplied inventory freezes the sector, TIC and product selection. Rerunning its first command updates archive timestamps and listing hashes. The committed Gaia cone replies preserve the measured field. The output retains the full difference matrices, fitted offsets, bootstrap results, two extraction results and input hashes.

## WD 1856 b in custom short-cadence apertures

The planet's published clock ranked first in blind searches of both Sector 86 target-pixel files after a noise-based small-aperture selection. In 20-second pixels the period was 1.40794276 d with 11 supported events and a midpoint phase offset of 0.22 minutes from the literature. In 120-second pixels it was 1.40793894 d with 10 supported events and a 0.005-minute offset. A fixed 1.5-pixel-radius aperture on each stamp missed the clock. These simultaneous reductions test the pipeline on a known positive control and carry no survey occurrence weight. [Blind results, input hashes and selected pixel coordinates](planet_target_pixel_results.json)

The short-cadence FLUX pixels already have SPOC background subtraction. A propagated Gaia coordinate sets the target location in each product's aperture WCS. Each alternative aperture sums a subset of pixels within two pixels of the predicted position. The subset maximizes expected target signal divided by the quadrature sum of pixel uncertainties, using a 0.85-pixel-width Gaussian integrated over each pixel and the median archived pixel variances on every 100th quality-zero exposure. The selection uses no ephemeris or in-transit flux. The five-pixel fast aperture returned a 7.14% observed box depth and the six-pixel two-minute aperture returned 4.72%. Neighbour models give nominal WD-light fractions of 7.88% and 5.13%, respectively. Varying the Gaussian width from 0.65 to 1.05 pixels spans 3.41–26.34% and 2.53–15.45%, respectively. Blending limits interpretation of the depth. The known planet remains outside the survey crowding cut.

Run from the parent result directory with the same Python dependencies as the primary survey.

```sh
PYTHONPATH=. python controls/planet_target_pixels.py controls/control_inventory.json controls/control-gaia-267574918.tsv cache controls/planet_target_pixel_results.json
```
