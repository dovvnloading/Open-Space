# Two later-sector phase predictions are not recovered in held-out TESS images

The phases chosen from the start of the later visits did not yield a secure recurrence of the 2018 dimmings in the held-out images of either asteroid. For (29489), 9 of 19 test cycles were fainter in the predicted window. For (3177), 7 of 17 were fainter. Their median depths were −0.0054 and −0.0030 mag, respectively. These measurements constrain a persistent signal at the tested later-season phases and amplitudes. Viewing geometry, intrinsically evolving events, and the (3177) spin/cadence alias remain possible explanations of the 2018 measurements.

| Asteroid / sector | Retrieved FFIs | Usable difference fluxes, entire visit | Held-out flux points | Qualifying cycles | Fainter in window | Median depth (mag) | Cycle-bootstrap 95% interval (mag) |
|---|---:|---:|---:|---:|---:|---:|---:|
| (29489), Sector 20 | 926 | 916 | 452 | 19 | 9 | −0.0054 | −0.0074 to +0.0682 |
| (3177), Sector 39 | 1,177 | 936 | 593 | 17 | 7 | −0.0030 | −0.0108 to +0.0139 |

The intervals are descriptive and assume cycles are exchangeable. The alternative fixed-radius aperture measured −0.0024 mag (9/19) and +0.0119 mag (9/17). Its corresponding intervals were −0.0326 to +0.0672 mag and −0.0367 to +0.0318 mag. The primary extraction's per-point residual scatter was 0.082 and 0.043 mag. Both sources had measurable rotational signals in the phase-training intervals. Artificial 0.050- and 0.100-mag events injected into the already-extracted test fluxes increased the measured depth by 0.051 and 0.098 mag for (29489), and 0.058 and 0.114 mag for (3177). This tests attenuation by detrending and scoring on the actual samples. It does not test loss upstream in the image extraction. The frozen phase ranked 41st of 64 circular offsets for (29489) and 46th for (3177), where the deepest offset ranks first. Adjacent offsets share many observations.

[Registered-test results](holdout.csv) · [Individual cycle depths](events.csv) · [Phase and injection sensitivity](sensitivity.csv)

## Observation and extraction

These are new reductions of Sector 20 and Sector 39 TESS full-frame images from the [MAST mission archive](https://archive.stsci.edu/missions-and-data/tess) accessed from the [public MAST S3 mirror](https://registry.opendata.aws/mast/). The Minor Planet Center already holds TESS station C57 positions and magnitudes from the same visits. The neighboring external-check report examined those catalog magnitudes before we opened these image cutouts. This experiment tests separately extracted image flux on a time interval held out from our own phase selection.

We interpolated C57 astrometric positions to each exposure midpoint in UTC, projected the position with each image's WCS, and extracted a 13-by-13-pixel square from the calibrated flux image. The accompanying acquisition program and manifest identify the archived exposures. The derived flux tables retain the midpoint, fractional pixel position, image-difference scatter, two measured fluxes, reference count and quality value. We retained files with chip coverage and with a C57 position within 0.065 days. The manifest identifies the downloaded frames exactly.

Fixed, same-detector reference stars register each cutout against its immediate neighbor. Integer translations from −4 to +4 pixels are compared using the median absolute pixel difference after subtracting the median difference. Gaps longer than 0.065 days start a new segment. For each science stamp, we use neighboring stamps in that segment where the asteroid is displaced between 2.5 and 12 pixels. We mask a 2.6-pixel radius around the asteroid in each reference, align the reference, subtract its median sky shift relative to the science image, and take the pixel-by-pixel median reference. At least three reference images, sufficient finite reference coverage, and zero quality flags are required. After reference subtraction, a residual background is estimated outside a three-pixel radius. The primary flux is the least-squares amplitude of a circular Gaussian with fixed 0.8-pixel standard deviation in a 2.8-pixel radius centered on the predicted location. A simple 1.65-pixel-radius sum on the same difference image supplies a second extraction. Both extractions share pixels and their errors are correlated.

## Frozen test

The Sector 3 event periods and widths were fixed. Only event phase in the later sector was fitted to an initial batch of image flux, separately for each asteroid. Fixed phases were committed in [12ff722](https://github.com/dovvnloading/Open-Space/commit/12ff72230b675a6f5bb136856e2f0eed62709182). The image and scoring code were committed in [f346799](https://github.com/dovvnloading/Open-Space/commit/f3467991fd63913a34b38952f0b4ee76967e4c4c), before inspecting the held-out image flux. A gap of one day separates the phase-training endpoint from the first permitted scoring observation. No event parameter is fitted in the test interval. The 0.8-pixel Gaussian is the primary result. The aperture is a sensitivity check.

Within the scoring interval, the program requires positive extraction flux, discards the highest 3% of image-difference scatter, removes a running median within ±15 hours, and fits three harmonics at the discovery primary-spin period. A four-scaled-MAD cut removes outliers before one refit. Each cycle needs at least two in-window points and four far from event phase (more than 0.27 cycles). The reported depth is the median of individual cycle depths. Bootstrap intervals resample cycles 10,000 times, treating them as exchangeable. The observing cadence and stellar backgrounds can correlate neighboring cycles.

The archived phase-training photometry (`training_29489.csv.gz` and `training_3177.csv.gz`) preserves the exact input snapshot, including reference-image availability at its boundary. Sector 20 contributed 408 training rows and Sector 39 contributed 285. The growing acquisition did not change either registered phase.

From the directory containing this report, run the following. The extraction command uses the output directory from the FFI retrieval.

```sh
python extract_photometry.py /path/to/cutouts --out /tmp/asteroid-flux
python replay_photometry.py --out /tmp/asteroid-check
python phase_scramble.py --out /tmp/asteroid-check
```

The replay defaults to the two compressed measurement tables next to the script. If the registered code lives elsewhere, supply `--locked-code /path/to/registered/code` to these commands.
