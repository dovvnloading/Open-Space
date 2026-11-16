# New independent lane: Mars Express VMC 2025-H2 global atmospheric transient search

Status: **prospective screening plan only; no discovery claim yet**. Started 2026-11-15. This directory is deliberately new: it does not extend the existing Mars craters, impact timing, or spacecraft-geometry investigations. It does not use the SPHEREx moving-object lane in open PR #102.

## Source and scientific question

Primary source: ESA Mars Express Visual Monitoring Camera (VMC), calibrated RDR Extended Mission 9 `MEX-M-VMC-3-RDR-EXT9-V1.0`, copied at NASA PDS Geosciences Node:
https://pds-geosciences.wustl.edu/mex/mex-m-vmc-3-rdr-ext9-v1/mexvmc_l3_0009/

The PDS catalog reports an **October 16, 2026** release extending coverage to **December 17, 2025**:
https://pds-geosciences.wustl.edu/missions/mars_express/vmc.htm

Question: Does the newly released **July–December 2025** VMC imaging contain an unusually large, spatially coherent and temporally rapid Martian cloud or dust event, independently supported by consecutive frames rather than calibration or viewing artifacts?

**Prioritize** events with plausible broad atmospheric significance; do not call generic frame-to-frame brightness changes discoveries. The initial observing unit is a dated VMC calibrated frame with supporting labels; a potential event requires multiple independent acquisition times.

## Frozen eligibility and falsification plan

1. Inventory the post-2025-06-30 calibrated VMC frames (latest newly released half-year), retaining image IDs, image dimensions, acquisition times, calibration units, filter and exposure metadata where available. Record eligible and excluded counts.
2. Examine adjacent independent acquisition sequences, controlling for disk position, illumination, off-planet background, Bayer pattern, saturation, exposure, and navigation. Do not difference unregistered images.
3. Define candidates using spatially coherent changes occupying a materially important portion of the Mars disk, not isolated pixels. Compare against matched same-sequence negative controls and another imaging epoch; report all inspected denominators.
4. For any strong candidate, search an independent instrument (e.g., MRO MARCI / Mars Express other frames) at comparable time and location; distinguish clouds from dust, geometry, and surface albedo. No physical or discovery claim without cross-sensor checks.
5. Archive compact reproducible measurement tables, PDS file URLs, SHA-256 hashes, analysis code and null tests. Document negative and incomplete analyses explicitly.

## Critical limitations

VMC images are often broad, color-Bayer, changing-geometry views. Uncalibrated brightness variation is **not** optical depth or a confirmed atmospheric event. A late data-release date is not a recent observation date. Positive-looking image differences are screening-only unless repeated in independently obtained data.

This initial checkpoint exists to preserve the unique lane and preregister controls before outcome selection.
