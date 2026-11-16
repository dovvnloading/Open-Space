# Newer photon pilot: trigger falsified by scan background and overflow

**No new impact established.** This is an executed, prospectively specified pilot on previously unexamined July 11, 2025 photons, not merely another known-event recovery. The frozen rule failed as a physical-transient detector in this observing regime.

## Data and measured results

First 500,000 rows of PJ74 S01, chosen mechanically before inspection. Retrieved 42 MB from the 1.288 GB FITS product, verified byte range and SHA-256. Native photon span is 2025-07-11 13:17:07.667 UTC through 7.652944 seconds later. Product label start instead refers to housekeeping coverage; do not substitute it for native photon time.

- 396,489 on-disc photons; 94,367 on-disc photons at 170–200 nm.
- 235 twenty-ms bins pass the 20-longwave-photon count cut; 81 also pass the fixed spatial-cell concentration cut.
- Those triggers advance through consecutive longitude cells with the scan. They are not 81 independent events.
- Every on-disc longwave photon is dayside; incidence-angle 1st/50th/99th percentiles are 55.14/56.96/62.14 degrees. Ordinary reflected-light and scan-exposure explanations have not been removed.
- 55 off-disc bins also exceed the count cut, reaching 375 longwave photons in one bin.
- All 20 time-shuffle trials produce zero localized triggers. This does **not** make real triggers significant: shuffling destroys scan/sky association. The null is unsuitable for physical discovery inference.

The frozen plan accidentally calls its latitude grid planetographic. Implementation uses the archive's explicitly planetocentric 400-km geometry; this terminology correction does not change a threshold or selection.

## Instrument checks actually executed

Parsed all nine FITS HDUs independently instead of assuming control-file layouts. Retrieved acquisition records, housekeeping and mask information.

- Four housekeeping samples overlap the photon pilot. All show `SCI_OVERFLOW=1`. The PDS label and instrument SIS define this as overflow during high-speed science transfer. This disqualifies completeness/rate claims unless loss behavior is modeled; it does not prove every individual photon is invalid.
- All four show hardware hot-pixel masking enabled, high-voltage status 1, no voltage limiting and checksum error 0.
- Planned mask `edgemsk1` covers the entire pilot; no photon lies inside its nonzero inclusive detector rectangles. Declared unmasked detector-area fraction: 0.64100647. This is a detector-area fraction, **not** sky exposure or detection completeness.
- Acquisition tag 51712 matches mask tag 0xca00. Frame-quality values 243/244 and science-quality bytes 246/248 are preserved, not interpreted as Boolean quality flags. Inspected SIS field descriptions do not supply their bit decoding.
- The product has no housekeeping sample preceding the first photons. Planned masks are not an independently verified history of actual per-photon mask state.

## Exact barrier to the next scientific claim

The data transfer itself is working. The barrier is scientific qualification: need a clean, preferably nightside acquisition; explicitly decoded acquisition/science-quality bits; actual scan exposure/neighboring-rotation background; and detector-frame injection/recovery with telemetry losses represented. A repeat of this naive trigger on more files would mainly accumulate scan/background candidates. No upper limit or impact rate follows from this pilot.

## Reproduction

From this folder, with NumPy, Astropy and curl:

```sh
python code/fetch_pj74_headers.py
python code/fetch_pj74_photons.py
python code/fetch_pj74_aux.py
python code/blind_pilot.py
python code/instrument_diagnostic.py
```

`newer-pj74/evidence.json` contains measurements, all triggers, shuffle counts, mask/HK checks and byte receipts. Binary data are re-fetched, not committed. The code uses no browser tests or hosted inference.

Primary references: [archive product label](https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DATA/ORBIT-74/UVS_S01_805506529_2025192_P74OBS_V01.xml), [UVS SIS](https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DOCUMENT/12029-RDRDP_SIS-01.PDF).
