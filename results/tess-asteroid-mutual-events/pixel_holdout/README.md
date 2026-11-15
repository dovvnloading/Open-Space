# Four later visits test whether the 2018 dimmings recur

We extracted 13 × 13 moving-target pixel stamps from 3,412 calibrated TESS full-frame images. Each later visit tests the **period and duration measured in 2018**, allowing the event phase to be set using early data from the later visit. This avoids extrapolating an uncertain ephemeris over years. Every depth compares the in-window and out-of-window magnitudes in one cycle. Positive values mean the target grew fainter at the predicted time.

| Asteroid | Sector | Held-out measurements | Fainter / covered cycles | Median depth (mag) | Whole-cycle bootstrap 95% interval (mag) |
|---|---:|---:|---:|---:|---:|
| (29489) | 20 | 452 | 9/19 | −0.0054 | −0.0074 to +0.0682 |
| (3177) | 39 | 593 | 7/17 | −0.0030 | −0.0108 to +0.0139 |
| (29489) | 73 | 114 | 3/5 | +0.0632 | -0.4199 to +0.2243 |
| (3177) | 22 | 290 | 4/8 | +0.0014 | -0.0071 to +0.0177 |

The two most precise later visits give median depths near zero in the fixed-phase tests. In sector 20 the upper end of the interval is below WQ’s 2018 depth of 0.128 mag. In sector 39 the upper end is below Chillicothe’s 2018 depth of 0.0477 mag. Injections into the **extracted photometry** verify that the scoring and detrending can recover fixed-phase dips at the original scale. These bootstrap intervals describe the sampled cycles and are particularly imprecise in the small sector 73 sample. [Individual cycles, frozen phases, artificial dips and comparison aperture](FIRST_TWO_VISITS.md) · [Two additional visits](other_visits.csv)

## Independence and observing geometry

For sectors 20 and 39, the original periods, the trained phases, and test cutoffs were committed before inspecting the held-out pixel fluxes. The reference-image extraction and scoring code were also frozen. For sectors 73 and 22, the dates, phase-fitting rule, period, window width and scoring procedure were committed before reading the derived test fluxes. The latter two use a 0.5-day separation between the calendar training cutoff and start of testing. Each side is reduced separately to avoid shared reference images. Previously inspected MPC magnitudes cover these same visits. These are held-out tests of new image reductions, and all TESS reductions within a given sector share exposures. [Phase checkpoint](https://github.com/dovvnloading/Open-Space/commit/12ff72230b675a6f5bb136856e2f0eed62709182) · [Frozen extraction and score](https://github.com/dovvnloading/Open-Space/commit/f3467991fd63913a34b38952f0b4ee76967e4c4c) · [Additional date splits](https://github.com/dovvnloading/Open-Space/commit/b334a32f46579c8d3646afd6dee9657ef647b784)

Geometric vectors from JPL Horizons quantify how different the observing directions were. From WQ, the direction to Earth differed from the 2018 visit by 97.7° in sector 20 and 76.2° in sector 73. For Chillicothe the corresponding angles were 148.5° in sector 22 and 122.1° in sector 39. Directions to the Sun also changed. The calculation uses the geocenter as a proxy for TESS at the representative Julian dates in the table. A binary pole and mutual orbit would be needed to predict where eclipses should disappear or return. These angles describe the geometry and do not independently identify a companion. [Vectors and dates](geometry.csv) · [Recompute angles and source URLs](geometry.py)

## Replay

The target track uses MPC observatory C57 astrometry. The 80-column uppercase satellite-observation lines contain target positions. The lowercase continuations contain spacecraft coordinates. The independent-checks folder records the original responses, hashes and extraction source. For these stamps we linearly interpolate observed right ascension and declination, project them through each FFI’s celestial WCS, and read the 13 × 13 stamp around the nearest pixel. We use the mean of the primary science image’s UTC start/end timestamps. We preserved fractional-pixel positions, the quality masks, calibrated flux pixel values and every accepted FFI name. The sample includes frames with nonzero quality flags, which the difference-photometry stage excludes. The images were selected near reported MPC positions and are a partial sample of the full sectors. [MPC observations and provenance](../external_checks/) · [MAST full-frame image archive](https://heasarc.gsfc.nasa.gov/docs/tess/data_products.html) · [Frozen FFI manifest and numerical checksums](ffi_manifest.csv.gz)

From this directory use Python, NumPy and SciPy for the offline tests. Reconstructing the pixels also needs Astropy. The replay used NumPy 2.3.5, SciPy 1.17.0 and Astropy 8.0.2. Replays of the raw archived images under other WCS library versions can differ slightly in the last digits of the derived subpixel coordinates. The original reduced measurements are included for offline analysis. [Frozen observing track](mpc_ephemeris.csv.gz)

```sh
python replay_photometry.py --out /tmp/tess-pixel-replay
python phase_scramble.py --out /tmp/tess-pixel-replay
python replay_other_visits.py
python geometry.py
```

To replay from the calibrated pixels, run the commands below. The first command verifies the numeric payload of each stamp against its SHA-256. The FFI filenames in the manifest uniquely identify the AWS archive products.

```sh
python retrieve.py /tmp/tess-asteroid-stamps
python extract_photometry.py /tmp/tess-asteroid-stamps --out /tmp/tess-first-two
python rebuild_auxiliary_photometry.py /tmp/tess-asteroid-stamps --out /tmp/tess-other-visits.csv.gz
```
