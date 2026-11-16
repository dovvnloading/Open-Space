# Mission 05 investigation log

## 2026-11-16 · Opening record, before inspection of reach measurements

- Base commit `3be0a9dabc756bb5e24df4ff4d4a81cc4e288f94`. Examined 121 advertised branch tips and 4,266 distinct file paths across their trees. Searched paths for SWOT, Hydrocron, Zambezi, hydrology, hydraulics, river discharge, width and anomaly variants. No matching path. Default-branch code search for SWOT and issue searches for SWOT, river discharge, hydrology and Zambezi returned no matches. Reviewed the six open PR titles and descriptions (#114, #108, #106, #105, #104, #102). No overlapping investigation surfaced. This establishes repository scope as of the check.
- Unit of observation is one SWORD reach, pass and UTC acquisition. Adjacent reaches on the same overpass share errors. Satellite discharge and its input width, slope and height also share errors.
- Primary source is NASA/CNES SWOT L4 Sword of Science River Discharge Version 3. Initial target is the upper Zambezi, where published river-authority observations and the documented 2024 drought allow a physical cross-check. A second basin may provide a comparison. Reach selection, all assessed candidates and exclusions will be logged.
- Test for changes in stage-conditioned width or slope or a seasonal discharge departure on repeated, quality-filtered passes. Report effect sizes, observational support, multiple testing and persistence limits. Compare with independent observations, and verify gauge use before interpreting any agreement with discharge. The known regional drought is a positive control.

[PO.DAAC release](https://podaac.jpl.nasa.gov/node/2037)
