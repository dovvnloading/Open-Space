# A long optical flare on UX Arietis during two MAXI X-ray scans

On 2021 October 22, UX Ari brightened **9.95%** above its local optical baseline. Its TESS-band equivalent duration was **1,604 seconds**, about **92 times** that of the UX Ari flare tabulated by Xu et al. in this same TESS sector. MAXI measured excess soft X-rays in two successive passes during the optical decay. The amplitudes, time-resolved fluxes and difference images support a large flare in the unresolved UX Ari system, but cannot identify which component flared. The MAXI alert itself had only a reported 2.0σ single-alert significance; these data establish a timed cross-band association, not a precisely determined X-ray peak or physical lag. The published X-ray alert named UX Ari. The additional evidence is the optical event, its scale, and the cross-band association. [MAXI 2021 alert](https://maxi.riken.jp/alert/novae/9509734289/9509734289.htm) · [Xu et al. catalog query](https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2bA/699/A322/tablea1&TIC=29701408)

A second pairing captures the 2018 August 16 outburst of AT Microscopii. TESS measured an optical peak **9 min 56 s before the MAXI alert**, and the source was still about 73% brighter than its local baseline at the alert time in the corrected light curve. Günther et al. already modeled this optical event as a two-component flare, outburst 43. Their first component peaks 1.39 minutes from the maximum TESS cadence used here. The measurement added here links that published optical outburst to a specific X-ray exposure. [MAXI 2018 alert](https://maxi.riken.jp/alert/novae/8346898070/8346898070.htm) · [Günther et al. catalog query](https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/AJ/159/60/table1&TIC=441398770&-out.max=500)

| Observable | UX Ari, 2021-10-22 | AT Mic, 2018-08-16 |
|---|---:|---:|
| TESS optical peak, spacecraft UTC | 20:39:56 | 19:12:40 |
| Peak fractional excess in corrected TESS photometry | 9.95% | 186.16% |
| Peak range across raw/corrected flux and baseline fits | 9.33–10.22% | 128.76–186.94% |
| TESS equivalent duration across measured window | 1,604 s (1,401–1,624) | 1,717 s (1,192–1,738) |
| MAXI scan center, ISS UTC | 21:52:59 | 19:22:06 |
| TESS excess interpolated at alert time¹ | 7.12% | 73.25% |
| MAXI 2–4 keV photon flux (photons cm⁻² s⁻¹) | 0.218 ± 0.073 | 0.170 ± 0.055 |
| MAXI 4–10 keV photon flux (same units) | 0.115 ± 0.055 | 0.077 ± 0.040 |
| Optical timing control, sampled sector time reaching the alert score | 87 / 15,148 (0.574%) | 7 / 16,778 (0.042%) |
| Difference-image centroid offset from quiescent source | 0.14 TESS pixel | 0.39 TESS pixel |

¹ Linear interpolation of two-minute integrations to the published alert time. MAXI's trigger timestamps and scan centers differ by seconds. The round-off on single-cadence optical peak positions is small relative to their approximately one-minute sampling uncertainty.

## How exceptional was the UX Ari flare?

Xu et al. list one sector-44 flare for TIC 29701408, peaking at BJD−2457000 = 2501.154, about 9.2 days before this event. Their table has three entries for UX Ari across all sectors. Integrating their sector-44 event between the published endpoints gives a 17.46-second equivalent duration. Applying its published energy of 3.53 × 10³⁵ erg with the same bolometric conversion to the October 22 event gives **3.24 × 10³⁷ erg**. The baseline and photometry variants span **1.84–5.01 × 10³⁷ erg**. This estimate assumes both flares have the same conversion between TESS counts and bolometric radiation. A changing flare spectrum, a different contributing component, or flux outside the integration window introduces additional uncertainty. The directly measured quantity is the TESS-band equivalent duration. The 92-times comparison is against one catalogued flare in the same sector, not against all flares from this system; the quoted bolometric energy is a conditional conversion, not a directly observed energy measurement. [Xu et al. (2025)](https://doi.org/10.1051/0004-6361/202451962)

The October 22 rise and decay sit on rotational modulation. At the same phase one orbit earlier, one orbit later and two orbits later, the local maxima in the same eight-hour window are 1.61%, 3.37% and 0.84%. The 2–4 keV scan following the alert gives 0.157 ± 0.069 photons cm⁻² s⁻¹, about 93 minutes later. The two-scan inverse-variance score against zero flux is 3.68. None of 217 other eligible neighboring scan pairs in the local window reaches that score. The two-scan interval was inspected after looking at the event. The statistic describes the local comparison and does not supply a survey-wide false-alarm probability. [UX Ari orbital period](https://www.aanda.org/articles/aa/pdf/2001/18/aah2495.pdf) · [MAXI scan data](https://maxi.riken.jp/pubdata/v7.7l/J0326+287/glcscan_lcbg_gsc0124578.csv)

## Selection and controls

We compared named active-star entries in the MAXI nova-alert listing with MAST TESSCut sky coverage and broad sector dates. An exploratory roster of 29 names matched 118 entries in the index from the start of TESS science. Four candidate overlaps emerged from the sector dates parsed in this screen. This roster and date screen are incomplete, especially in later sectors. The entries with good TESS data bracketing the alert were UX Ari and AT Mic. For the 2020 May 1 SS Boo alert, the nearest good cadence was 10.56 hours away. For the 2021 December 3 σ Gem alert, it was 21.22 hours away. Neither gap constrains an optical counterpart. [MAXI index](https://maxi.riken.jp/alert/novae/index.html) · [TESS sector notes](https://heasarc.gsfc.nasa.gov/docs/tess/data_release_notes.html)

Use SPOC two-minute cadence data with QUALITY zero and finite PDCSAP flux. Subtract the supplied barycentric time correction and 69.184 seconds to obtain spacecraft UTC to millisecond-scale time-standard accuracy. TESS-to-ISS light travel remains in the comparison at up to roughly two seconds. Detrend UX Ari against the two flanks 6–12 hours from the alert and integrate from −5 to +3 hours. Use 1.5–4 hours and −0.45 to +1.35 hours for the shorter AT Mic profile. The primary fits use a line. Sensitivity fits use linear and quadratic backgrounds, flanks widened in radius by a factor of 1.33, and both raw SAP and dilution-corrected PDCSAP flux. Interpolate trapezoidally across gaps shorter than 246 seconds only. The sensitivity ranges describe model choices, not confidence intervals. AT Mic's large raw/corrected amplitude difference warrants quoting both flux products.

For a whole-sector optical control, divide each corrected light curve by its 12-hour running median (at least 200 samples) and take a 20-minute running median of the resulting excess (at least seven samples). Require eight or more good cadences in 0.08-day side windows centered 0.2 days on either side. Compare the good cadence nearest each external X-ray alert to all eligible good cadences. Adjacent cadences and flares are correlated. The fractions describe time occupied by equal-or-higher scores under this rule. They carry no correction for selecting stellar alerts or exploring multiple targets.

The archived MAXI case reports quote individual alert significances of 2.0σ (UX Ari) and 2.9σ (AT Mic). Our scan checks use the source-position 2–4 keV monitoring product. These are reductions of the same MAXI photons. At UX Ari, 0/233 comparison scans exceed the trigger scan's flux/error ratio. At AT Mic, 1/363 do. Comparison scans have quality flag zero, finite positive flux error, and lie 0.5–16 days from the trigger. The sparse X-ray exposures (48 and 51 s) cannot locate the X-ray maximum or measure a physical peak-to-peak delay. [MAXI scan format and binning](https://maxi.riken.jp/top/lc_readme.txt)

The peak-minus-pre-event images use all calibrated target stamp pixels. Flux-weighted centroids include pixels exceeding 2.5% of each image maximum. The AT Mic residual is compact, brightens multiple adjacent pixels, and falls within 0.39 pixel of the quiescent photocenter during a sector with documented pointing problems. The UX Ari residual lies within 0.14 pixel. The TESS pixel scale is about 21 arcseconds. Neither measurement distinguishes individual members of the stellar systems. [TESS sector 1 instrument notes](https://heasarc.gsfc.nasa.gov/docs/tess/sector1_summary.html)

## Reproduction

With Python 3, NumPy 2 and pandas, run the following commands in this directory. The first command downloads the source products, checks their identities and rebuilds the optical windows, difference images, MAXI scan snapshots and measurements. The second computes the reference-flare comparison and same-orbital-phase controls.

```sh
python analyze.py --out rebuilt
python compare_uxari.py --out rebuilt
```

`measurements.json` records the source URLs, optical-file SHA-256 hashes, target-pixel byte ranges and hashes, and MAXI snapshot hashes. The pixel frames and light curves agree on cadence timestamps to within 10⁻⁷ day. The MAXI snapshot columns are MJD center/start/stop, flux and one-sigma error for each of 2–20, 2–4, 4–10 and 10–20 keV, followed by quality. Flux is measured in photons cm⁻² s⁻¹. The optical CSV files contain calibrated electrons s⁻¹, spacecraft UTC JD−2457000, and local fractional excess. Difference-image values are electrons s⁻¹, with rows increasing in detector y. Multiple publications and catalogs report flares on both systems. Searches by date, target, TIC, and instruments did not locate a previous joint measurement of these two specific associations. Broader publication priority remains open.

This research used MAXI data provided by RIKEN, JAXA and the MAXI team, and NASA TESS observations archived at MAST.
