# Measured positive control, not a discovery

A real calibrated photon slice recovers the published 2020 April10 impact-like flash from Giles et al. (2021), https://doi.org/10.1029/2020GL091797. This is an independent computation on the **same observation**, not independent observational confirmation.

- Product: `UVS_S01_639769969_2020101_P26OBS_V01.FIT`, PJ26.
- Native photon list:14,476,358 rows,84 bytes each, beginning at byte4,320,000.
- Frozen slice: rows12,332,139 through12,444,380 inclusive;112,242 photons within ±2 seconds of12:57:10 UTC.
- Target gate: published53°N,200°W with ±3° geometry window. Geometry is the archive's400-km reference shell, not the paper's inferred225-km altitude.
- Peak longwave1-ms bin:66 photons at12:57:10.2835 UTC (bin center).
- Within80ms of that peak:2,885 photons in target gate, including993 at170–200nm and1,066 at130–160nm.
- Median latitude/longitude:53.0023°N,200.0847°W.

This reproduces the spatial/time localization and presence of long-wavelength photons. It does **not** reproduce the paper's calibrated temperature, altitude, mass or impact rate. Raw band counts are not a calibrated spectrum. Nearby off-window counts have different scan exposure and must not be used as a significance denominator. Masks, full housekeeping, independent scans, synthetic-injection recovery and survey-wide false-alarm controls remain required.

## Replay

Install NumPy and Astropy from the standard package registry; curl is required. From this folder:

```sh
python code/extract_control.py
python code/summarize_control.py
```

Default retrieval uses one approximately9.43MB HTTP byte range, checks Content-Range/length, validates time bounds and monotonicity, and records its SHA-256. `--rediscover-bounds` independently rederives bounds with additional requests. Cached source bytes and `control-photons.npz` are intentionally not committed. `acquisition.json` identifies source bytes; `control-summary.json` preserves measured results.

## New search scope and freshness

The actual archive lists PJ74 July11 2025 photon observations beyond the guide's old index throughPJ70. The PJ74 S01 label covers13:17:09–13:32:17 UTC and was modifiedDecember12 2025. Do not confuse the October20 2026 **derived** release with a new photon-data release. The new scientific lane is a frozen post2020 search for additional localized continuum-rich transients, not rebranding this known event as new. No such new event has yet been established.
