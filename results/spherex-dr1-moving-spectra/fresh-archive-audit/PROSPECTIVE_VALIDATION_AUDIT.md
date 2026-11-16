# Prospective-validation audit and independent optical context

## No untouched SPHEREx D5 epoch in the available query
A fresh MOST query requested March 1, 2025 through November 16, 2026 for C/2025 R2. The archive clamps the interval to April 11, 2025 through October 16, 2026 at 12:00 UTC. It returns the same 348 intersections, including the same 57 D5 exposures. Comparison by observation ID and detector finds **zero new D5 records**.

All 57 D5 records across 17 dates had already been attempted in the full screen. Pixel-level photometry was obtained for 56. Therefore none of these dates can honestly be relabeled as a prospective untouched epoch. The January 2 family comparisons remain adaptively reused later-epoch tests, as disclosed in FINAL_STATUS.md.

The sole failed download is January 1 at 03:20:00.416. The three other January 1 exposures contain only one valid inclusive aperture (4.376803 microns); the others fail core validity. One valid datum cannot fix both continuum and molecular amplitude before predicting a new channel. Retrying that missing exposure would not supply an independent-epoch shape test. No new pixel retrieval or post-hoc selection was performed for this audit.

A genuine prospective test must wait for new release data or use an independent instrument. Freeze a response convention, physical-template family/temperature rule, quality mask, extraction, uncertainty model and success/failure thresholds before inspecting such pixels. Do not claim that selecting a new subset of already processed dates restores blinding.

## Independent archive attempt
The [current LCO developer documentation](https://developers.lco.global/) uses `target_name`, not the legacy `OBJECT` query parameter used earlier. We corrected the query and tried both C/2025 R2 and 2025R2 plus public=true. Both returned estimated count16 and an empty result list; a proposal/date query returned count0. These results neither identify 16 accessible SWAN exposures nor establish that no private/other-alias exposures exist. No authentication or outreach was attempted. The saved JSON documents the responses.

## Independent optical observations actually retrieved
The public [NEOx MPC-format download](https://neoexchange.lco.global/target/156360/measurements/mpc/download/) provides 2,494 lines, including 131 reported magnitude records between December 20, 2025 and January 5, 2026. These are aggregated observations from multiple observatories, not all LCO measurements. Dates, filters and observatory codes are retained in `independent-optical-context.json`.

For example, site L27 reports three r-band magnitudes on January 2 (17.8–18.0, median17.8) and three on January 4 (17.7–17.8, median17.8). Other sites/filters span much wider values; aperture definitions, zero points and uncertainties are not provided in these MPC records. These entries establish contemporaneous optical observation, but cannot deliver calibrated dust evolution, gas/dust ratios, narrowband infrared validation or novel chemistry. No physical fading/outburst claim is made.

## Conclusion
There is no valid pristine temporal shape test available in the currently retrieved D5 archive. Independent optical records exist and were retrieved, but no independent contemporaneous IR spectrum has been verified. The candidate remains exploratory. This is an explicit evidence limit, not a successful discovery-validation result.
