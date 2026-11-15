# XMM-Newton imaging of the 2018 M51 ULX-1 brightening

The three EPIC detectors locate the low-to-high transition at M51 ULX-1 (13h29m39.94s, +47°12′36.6″). Its companion ULX-2 lies 8.42″ to the northwest. Di Stefano et al. studied the distinct M51-ULS-1, at a different sky position, in the same observation. The parent report contains the Chandra timing and relevant literature.

In ObsID 0824450901, profiling a common sharp transition over three independent cameras yields geocentric barycenter-corrected MJD(TDB) 58252.13450. The profile envelope at ΔC ≤ 3.84 is [58252.13375, 58252.14425]. This is a conditional sharp-step localization and lacks coverage calibration for a changepoint. A smooth transition model gives midpoint 58252.1402 and a 10–90% transition time of 0.076 d (6.6 ks). Its profiled midpoint grid spans 58252.128–.152 at ΔC ≤ 3.84. Both timing envelopes are conditional on their assumed light-curve shape. The MOS1, MOS2 and pn best splits are 58252.13250, .14500 and .14425. At the common split, observed pn source rates rise from 1.85 ± 0.67 to 15.67 ± 1.15 counts/ks in 0.9–7 keV. Source amplitudes and background can change independently on either side. Freeing ULX-1 individually improves Cash by 55.5, 47.3 and 106.5 in MOS1, MOS2 and pn, respectively. Labels describe a brightening. An eclipse egress is one interpretation.

The next pn visit covers 58264.72976, the first 12.59526-day return or second 6.29763-day return from this transition. In a 22.03 ks interval before that instant, ULX-1 gives 15.53 ± 1.32 counts/ks. Afterward it gives 11.50 ± 1.40 counts/ks over 16.51 ks. The 2018 May low-to-high transition did not repeat at that cycle count. MOS was in partial-window mode in the three subsequent 2018 visits and the target had no local MOS event coverage there.

## Reproduction

Python dependencies are NumPy, Pandas and SciPy. Rebuilding the compact inputs also requires Astropy. Run from this directory.

```sh
OPENBLAS_NUM_THREADS=1 python analyze.py
OPENBLAS_NUM_THREADS=1 python egress.py
OPENBLAS_NUM_THREADS=1 python transition_checks.py
OPENBLAS_NUM_THREADS=1 python ramp_sensitivity.py
OPENBLAS_NUM_THREADS=1 python robustness.py
OPENBLAS_NUM_THREADS=1 python recurrence.py
```

The inputs are lightly rounded per-event offsets, energies and barycenter-corrected times, CCD- and flare-screened good-time intervals, and individual reference translations. FITS identifiers, retrieval URLs and SHA-256 hashes are in `xmm-products.csv`. To regenerate the compact inputs from raw products, first obtain the [Chandra reference image, ObsID 13813](https://heasarc.gsfc.nasa.gov/FTP/chandra/data/byobsid/3/13813/primary/acisf13813N002_evt2.fits.gz) (SHA-256 `c11ef86ae0ff9db0605d198df73804293c9c65cc9ac1ba334de3ad9520ed0155`), then run this command.

```sh
python rebuild_inputs.py --raw-cache ./archive-cache --chandra-evt2 /path/to/acisf13813N002_evt2.fits.gz --output ./regenerated --manifest xmm-products.csv
```

## Data selection and models

The archive screening used public HEASARC XMM PPS EPIC event lists. Select 0.5–7 keV FLAG=0 events with pn PATTERN ≤ 4 and MOS PATTERN ≤ 12. Project events onto a tangent plane centered on ULX-1. Apply the detector-specific median sky translation from isolated Chandra field sources within 310″; individual offsets are retained in each references file. Retain the source CCD's standard GTIs, intersected with 500 s bins whose 2.5–7 keV counts in the 120–400″ annulus lie below 1.8 times the positive-bin lower quartile. The compact dataset includes events within 38″ of ULX-1 before registration. Time conversion adds the interpolated Earth–barycenter light time. The omitted spacecraft light time is below a second.

Fit unbinned local photon positions in a 38″ square with independent amplitudes at the two fixed Chandra positions. Both use a circular King profile with core radius 4.5″ and exponent 1.5, plus a local positive bilinear background. Exposure is the integrated duration of filtered source-CCD GTIs (not corrected for instrument dead time). Per-visit, single-detector rate ratios avoid cross-camera throughput differences. Counts are expected counts in the fit square. Quoted count errors come from the inverse observed information matrix and are approximate close to a zero-amplitude boundary. The empirical PSF approximation is tested by varying the PSF parameters, extraction size, energy band and registration, alongside the counts in a 3″ aperture compared with three equal-radius azimuthal controls around the neighbor. These tests are in the sensitivity and azimuth tables.

`xmm-alias-rates.csv` folds the XMM data on integer cycle counts between the Chandra egresses at 56180.216165685095 and 59373.11583257826. Ingress 56191.61960343367 sets the putative width. Each contact is trimmed by 0.08 d. Rows require ≥1 ks and ≥25 photons in the fit square. Columns labeled `in` and `out` refer to predicted eclipse and bright intervals under that particular clock. The `delta_c` column compares adding versus removing ULX-1 in the given interval. A source can remain faintly visible in eclipse. The empirical test is whether predicted deep intervals instead have bright-state rates. Other cycle counts and variable eclipse geometry remain possible.
