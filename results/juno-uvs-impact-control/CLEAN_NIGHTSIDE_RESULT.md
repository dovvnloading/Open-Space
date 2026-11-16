# Metadata-selected nightside pilot completed: no trigger

After the failed dayside/overflow pilot, housekeeping was screened **before new photon inspection** in all five remaining PJ74 products. All had only sparse overflow-clear periods. Long clear blocks associated with planned dark masks were excluded from science selection. The earliest sustained 20-second clear interval under a normal mask was selected: S02 P74OBS, July 11, 2025, 13:56:37.207–13:56:57.207 UTC.

## Measured outcome

- Retrieved 9,787 photons from frozen rows 12,431,751–12,441,537 inclusive, using 822,108 source bytes.
- Native photon times span 0.9264–19.9993 seconds within the selected window. This is not a claim of continuous exposure for 20 seconds.
- 3,016 photons have valid on-disc geometry. **All are nightside.** Incidence-angle 1st/50th/99th percentiles are 137.92/154.23/156.80 degrees.
- 404 nightside photons fall in the 170–200 nm band.
- The unchanged trigger threshold of 20 longwave photons in a 20-ms bin was never reached: the observed maximum was 5. There are **zero candidate triggers**, not a new discovery.
- Twenty conditional-label diagnostics also yield zero triggers, with maxima of 4–6. Labels were shuffled within 100-ms strata, 16-pixel slit-coordinate segments and night/day status. This preserves local count/scan structure better than global time shuffling, but breaks the exact detector-X/wavelength relation and is **not a calibrated physical false-alarm probability**.

## Checks performed rather than deferred

Eleven housekeeping samples bracket the window. All have overflow=0, checksum error=0, high voltage on, limiting=0 and safety inactive. Acquisition records give frame-quality value 1. We preserve science-quality bytes without inventing a bit interpretation.

The actual acquisition tag is 0xca06, matching `edgemsk1` in this file's mask table. No photon falls in its masked detector rectangles. Actual tags take precedence over planned wall-clock transition times.

Detector STIM is active. Ancillary data locate its two stable detector centers at (84,67) and (1989,170). Sensitivity boxes of half-width 0, 4, 8 and 16 pixels contain **no nightside science photons**, so these STIM centers do not explain or alter the nightside count result. The largest boxes contain 364 photons outside that nightside set.

“Telemetry-clear” is not proof of perfectly lossless calibrated exposure. This small negative pilot establishes that the selection and controls execute on real nightside data; it gives **no impact-rate limit, occurrence rate, survey completeness or discovery claim**.

## Reproduce

With NumPy, Astropy and curl, from this folder:

```sh
python code/screen_pj74_housekeeping.py
python code/extract_clean_interval.py
python code/fetch_clean_aux.py
python code/conditioned_pilot.py
python code/clean_instrument_check.py
```

The pinned photon header is committed. `clean-screen/clean-evidence.json` preserves metadata selection, measurements, conditional diagnostics, instrument checks and hashes. The metadata screen reads 12-MB tails because 5-MB tails did not reach housekeeping. Source binaries are deliberately not committed.

Next scientifically meaningful expansion: apply the fixed metadata/geometry qualification across independent rotations and acquisitions, construct response-aware injections, and validate an exposure-conditioned null. More unqualified count triggers would not be scientific progress.
