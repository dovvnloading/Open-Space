# New investigation: orbital scale in a public exoplanet catalogue

**Status: preregistered screening plan, not a physical discovery claim.** This investigation is independent of the existing Open-Space target campaigns. The upstream Open Exoplanet Catalogue table is a 2023-10-15 snapshot and must not be mistaken for a current census.

## Locked input

- Open Exoplanet Catalogue (OEC) generated CSV: `https://raw.githubusercontent.com/OpenExoplanetCatalogue/oec_tables/bf1efad9e4bdc9f0887fe392b77308c3585e263f/comma_separated/open_exoplanet_catalogue.txt`
- Table source commit `bf1efad9e4bdc9f0887fe392b77308c3585e263f`; source Git blob `db832c9d70e953ca6e1851122e1d46578c3b5d58`.
- Exact available fields and units in `https://github.com/OpenExoplanetCatalogue/oec_tables/blob/bf1efad9e4bdc9f0887fe392b77308c3585e263f/FIELDS.md`.
- Screening eligibility set **before** selecting follow-up objects: confirmed, nonbinary (`binaryflag=0`) with positive orbital period (days), semimajor axis (AU), and host mass (solar masses). No missing values imputed.

## Primary question and independent falsification

Can the reported orbital period, semimajor axis, and host mass satisfy Kepler's third law? The expected semimajor axis is `a_K = ((M_star + M_p) * (P/365.25)^2)^(1/3)` AU (mass in solar masses, planet mass when reliable). Compute `Q = a_reported^3 / ((P/365.25)^2 * M_star)` for all eligible rows. Use the conservative massless approximation to flag gross mismatches, and explicitly check massive-companion exclusions.

### Provisional results and checks to complete

Initial deterministic screening: 5,414 source rows, 2,275 eligible. **70** had `Q > 1.5` or `Q < 1/1.5`; this is **not** 70 confirmed errors. Kepler-445 b/c/d each imply ~16 solar masses for a host listed at 0.18; Kepler-446 b/c/d each imply ~16 solar masses for a host listed at 0.22. Possible shared origin: converting `a/R_star` as though the stellar radius were one solar radius. Test this mechanism explicitly against the measured `a/R_star` and stellar radius for each planet in the 2015 discovery paper, and against NASA Exoplanet Archive per-planet solutions. Screen the full table for other coherent groups and errors in source mass, period, axis and conversion. Retain all nonconfirmed flags as candidates.

**Falsifiers:** If the raw source catalogue differs from generated CSV, isolate an export issue; if an independent stellar/planet model justifies the quoted axes, retract the correction; if an alternate star/planet association explains the discrepancy, do not call it a unit conversion. Compare with positive control planets and test against numerical tolerance and binary exclusions. Publication priority is unknown.

## Deliverables

A standalone standard-library replay script, pinned-data provenance and bounded eligibility counts, machine-readable flagged sample, independently corroborated high-impact corrections, unit tests, and a limitations section. Preserve this early commit as a branch checkpoint even if the hypothesis fails.
