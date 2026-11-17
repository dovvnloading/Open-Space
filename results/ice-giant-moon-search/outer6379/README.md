# Uranus, 155,000–649,000 km, 2025 February 2

**None of 246 first-pass positive orbital templates placed a localized source on the track in both visits.** The 12 native NIRCam exposures cover a larger field than the inner-moon search. Of 447,559 tested circular, equatorial tracks, 240,362 had usable positions in both chronological halves of the F150W2 sequence and both short F150W frames. Outside 584,000 km, 2,587 of 91,397 tracks passed these geometry cuts and zero met the initial brightness threshold. The declining detector coverage makes the outer edge particularly incomplete. This exploratory screen has no measured limiting magnitude. [Per-radius counts](coverage.csv) · [Threshold crossings](crossings.csv) · [Individual control tests](vetting.csv)

The first ten frames and two short frames come from different program 6379 visits on the same date. Independently measured disk centers, the native WCS, and 17 emitted-time samples per exposure define each orbit trail. The radial grid spacing is 2,000 km and the azimuth spacing is at most 1,400 km in the orbital plane. Known moon paths, detector edges, pixels flagged beyond a corrected jump, the inner 285 pixels and the immediate planet/ring region are excluded. Local 128-pixel tiles normalize the difference of two Gaussian smooths. Thresholds require a summed score of 7 in the ten deep exposures, 3 in each half, 5 in the short pair and 2 in each short exposure. These correlated scores are screening statistics. The same cuts produced 104 inverted crossings. [Geometry and masks](screen.py) · [Inputs with full-file SHA-256](inputs.csv)

At each surviving track, the second test compares the predicted position to 48 displaced paths at 6, 10 and 15 pixels and allows nine common offsets over ±2 pixels. It requires local contrast of 5 scaled MADs in each visit and 2 in each deep half. The test was repeated after subtracting separately estimated fixed-sky and planet-fixed backgrounds from other exposures separated by at least 46 minutes. No track passed even the original localized-contrast test. All three comparisons were simultaneously evaluable for 241 of the 246 cells. Several neighboring positive cells sample the same structure. [Diagnostics](vetting.json) · [Stationary-background implementation](vet.py)

## Reproduction

Download the 12 original CAL files in the input list using the MAST URIs. Rebuild the February visit 001 ephemeris from the [deep-visit observables](../jwst6379_9482/f150w2_6379_observables.csv) and visit 002 from [short-visit observables](short_observables.csv). The latter uses a short URA184 kernel excerpt for body 75052 and public Horizons for the other bodies. Raw response and kernel hashes are in [the receipt](short_sources.csv). With both ephemeris NPZ files and the CAL directory:

```sh
python outer6379/screen.py CAL_DIR FEB_VISIT_001.npz FEB_VISIT_002.npz
python outer6379/vet.py CAL_DIR FEB_VISIT_001.npz FEB_VISIT_002.npz
```
