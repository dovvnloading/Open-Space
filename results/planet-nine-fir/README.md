# Follow-up of published far-infrared distant-planet candidates

**The Phan et al. pair has no independent moving-source recovery in the public images searched so far. Two DES nights constrain bright members of the assumed 500–700 au family. The result remains conditional on the historic astrometry and the planet's optical brightness.**

## Completed measurements

The IRAS catalog independently gives 0.2396 Jy at 60 μm and 0.5164 Jy at 100 μm. The published single-scan AKARI 90 μm detection is 0.27 Jy. The public all-scan map gives −0.058 Jy at the reported position. The original positive scans use different processing and the public coadd combines detections and misses. A W4-bright galaxy 32.86 arcsec from the IRAS position is a specific stationary counterpart worth testing. [Historical catalogs, arrays and photometry](historical/README.md)

The nominal positions imply retrograde ecliptic inclinations of 115.03–117.33° on the reconstructed grid. At 00:00 UTC on 2026-11-17 its geocentric footprint spans RA 34.9358–35.0616°, Dec −49.8341 to −49.7047°, with sky rates 0.197–0.278 arcsec per hour. Ranges describe correlated sky positions. [Assumptions, sensitivity and validations](orbits/ORBITS.md) · [Monthly ephemeris](orbits/ephemeris.csv)

![Twenty-eight possible paths illustrate the annual parallax loops. Outlines mark two dates for the nominal historical positions.](orbits/search_region.png)

| Data | Independent measurement |
|---|---|
| DES, 2016-12-30/31 | No detection above the moving-source threshold. Real-image injections measure a steep loss of completeness between r = 23 and 24 AB. [Search results](archives/difference_result.json) |
| NEOWISE W2, 2023–2024 | No repeated signal among 24,118 eligible templates of 24,132. Out of 256 separate test tracks, 240 were recovered at W2 = 14 Vega and 35 at W2 = 15. [Images, masks and calibration](wise/README.md) |
| Rubin | Zero associated public ALeRCE alerts in the query rectangle; no image-level magnitude or exclusion fraction. [Receipts](rubin/README.md) |
| DECam follow-up, 2026 | Twelve 900 s r exposures labeled for the candidate are in the archive, with science access scheduled for 2028. [Archive dates and identifiers](archives/dedicated_2026.csv) |

## Other candidates

We propagated both Chen AKARI candidates, the Rowan-Robinson IRAS proposal and the 20 strongest published ACT trajectories. [Source table and monthly predictions](other/README.md) · [WISE field-by-field screening](wise/compact/field_summary.csv)

## Pending or uncalibrated at this checkpoint

Final Euclid, SPHEREx, historical scan timing, optical sensitivity to uncertain centroids and whole-family exclusions. No posterior measure is supplied by two historical positions. Nominal survey depth, alert absence and footprint overlap do not furnish injection limits.
