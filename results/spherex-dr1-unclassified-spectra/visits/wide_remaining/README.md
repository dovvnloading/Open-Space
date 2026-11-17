# Native exposures at every original broad-field trigger

**None of the fitted peak residuals reproduced with the same sign at diagnostic score 3 in the 25 independently measurable QR3 windows.** Across 52 targets, DR1 supported all three bands for 42 and QR3 for 25, with 21 measurable in both. Six targets lacked a usable local contrast in both releases. Five DR1 and four QR3 windows had at least two clean dates in each of the three bins. No complete DR1 window reproduced a same-sign fitted residual at score 3 either.

| Archived cutout outcome | DR1 | QR3 |
|---|---:|---:|
| Selected cutouts at all 52 peaks | 415 | 230 |
| Clean source, blank apertures and science wavelength | 244 | 144 |
| Source aperture flags / coverage | 126 | 60 |
| Too few eligible blanks | 35 | 21 |
| Science wavelength outside target bin | 10 | 5 |
| Unresolved read failures | 0 | 0 |

One target displays a repeating local emission feature with a negative fitted residual. UGC 08856 NED01 has a +47% coadd local contrast at 3.383 µm, +41% in native DR1 and +79% in held-out QR3. The signed native diagnostic scores are +3.77 and +12.93. Its spectroscopic redshift 0.031062 places the familiar 3.29 µm aromatic emission band at 3.392 µm. The positive local flux contrasts are consistent with that identification. QR3 has one usable date per bin; source blending, non-simultaneity and the sparse sampling limit the comparison. The negative template residual is not evidence for a newly identified absorption band. [Catalog redshift and counterparts](../../wide/catalogs/ned_matches.csv)

The full broad-field model triggered on 52 of 59 eligible sources. The main follow-up list already included the original peak channel for twelve of those 52 sources. This extension investigated the other forty source–channel pairs at the exact original peaks, including sources previously measured at a different peer-adjusted wavelength. The `targets.csv` in this directory lists the forty. The combined `all_raw_trigger_checks.csv` preserves both releases for every original peak, the routing to its source-level exposure receipts and the reason an individual measurement failed when the masks or coverage were insufficient.

The independent observing set is QR3. DR1 native exposures also test original coadd construction. For the forty additional peaks, the selection considers up to three distinct MJD days per release and wavelength bin, using the same footprint index, detector-calibration grid, frozen wavelength windows, selection rules, native pixel flags and complete-aperture requirements as the main follow-up. The original twelve overlaps retain their main selection of up to six days per bin. When a release/bin has no usable measurement after the initial pass, the quality-only supplement tries the remaining eligible footprints nearest the fixed target wavelength, including same-day alternatives, until that bin has one clean measurement or the options run out. Supplemental receipt rows are marked `missing_bin_alternate`. One blue, core and red measurement can produce an interpretable local contrast. The dates in those three bins need not coincide. `all_bins_2plus` requires two or more retained dates separately in all three bins.

`same_sign_ge3` compares the native local-continuum contrast sign with the **fitted template residual** sign. `local_sign_ge3` additionally checks the **unmodeled coadd local-contrast** sign from `wide/coadd_checks.csv` at the same wavelength. A model can overpredict the strength of an emission band, producing a negative residual alongside a positive local feature. Scores are diagnostic and are not trials-corrected significance levels. The two coadd and native local-continuum estimators have different band windows. Wavelength coverage can differ across releases. An unmeasurable release or epoch does not supply a negative detection.

The raw per-image CSVs give URLs, hashes, dates, nominal selection and recalibrated science wavelengths, aperture coverage and flags, sky errors and signed fluxes; source JSONs aggregate the finite measurements. Inputs and individual follow-up criteria are described in [the native exposure methods](../README.md). From the result directory, after preparing the wide-field footprint and science-calibration caches as described there, replay the archived image choices and verify their hashes first.

```bash
python visits/remaining_wide.py /tmp/spherex-wide-index /tmp/spherex-native-wcs --replay
python visits/summarize_wide_triggers.py
```

To make a fresh selection against the live archive, use a separate checkout. The following commands write new source-level receipts in that checkout.

```bash
python visits/remaining_wide.py /tmp/spherex-wide-index /tmp/spherex-native-wcs
python visits/remaining_wide.py /tmp/spherex-wide-index /tmp/spherex-native-wcs --retry-errors
python visits/complete_bins.py visits/wide_remaining/targets.csv wide/sources.csv /tmp/spherex-wide-index visits/wide visits/wide_remaining /tmp/spherex-native-wcs
python visits/remaining_wide.py /tmp/spherex-wide-index /tmp/spherex-native-wcs --retry-errors
python visits/summarize_wide_triggers.py
```

The output summary reports figures for all 52 original peak locations. Use it separately from the two-field native follow-up counts because twelve pairs occur in both.
