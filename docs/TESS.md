# Search for periodic dimming in TESS photometry

The transit worker searches original two-minute SPOC light-curve products from NASA's Transiting Exoplanet Survey Satellite. It selects stars from archive metadata, measures a period in one sector, and asks whether the dimming returns at its predicted times in another sector. Every input FITS file is pinned by MAST URI and SHA-256. The results record rejected hypotheses, coverage, depth estimates and individual event measurements. [MAST TESS collection](https://archive.stsci.edu/missions-and-data/tess) · [Two-minute product inventories](https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html) · [Worker](../open_space/tess.py)

## Sample and controls

MAST lists 7,074 TICs with standard two-minute light curves in both Sectors 1 and 2. Sort their TIC identifiers by a deterministic hash, then take a contiguous page. The pilot uses the first 24 targets. A separate named control, WASP-18 b, uses Sectors 2 and 3 with exactly the same search and fixed-epoch validation. The frequency grid, durations, thresholds and veto rules were specified before examining pilot light curves. The 96-target extension takes the next disjoint page under the same rules. This sample covers the short-cadence target list and a particular sector overlap. It cannot give an all-sky planet or variable-star rate. [Search rules](TESS_PROTOCOL.md) · [Pilot source IDs and digests](../data/tess-plan.json) · [Extension source IDs and digests](../data/tess-expansion-plan.json)

Across the pilot and named control we verified 50 products totaling 101,047,680 bytes and retained 892,290 good raw cadences. The blind control search measured WASP-18 b at 0.94156265 days, within 0.012% of the 0.9414522-day catalog period. Its observed depths were 0.97% and 0.94% in the two sectors. Twenty-seven first-sector and nineteen second-sector events had sufficient local coverage and supported the dimming. A positive control validates this strong signal and time convention at one point in parameter space. [Exo.MAST planet measurements](https://exo.mast.stsci.edu/exomast_planet.html?planet=WASP18b) · [Pilot measurements](../results/tess-pilot.json)

Two of the 24 exploratory stars produced a strong recurring dip, each accompanied by a second dip near the opposite orbital phase. TIC 234495456 measured 0.76084 days and TIC 234338669 measured 0.91902 days. Both were flagged in both sectors. None of the initial exploratory targets cleared all primary rules. A box search can select a harmonic of a stellar variability period. The extended readout includes a held-out model comparison for that problem.

Injecting two-hour boxes into the unmodified flux of the first exploratory star, applying the masks and detrending again, recovered 7 of 9 signals with the right period, phase and two-sector disposition. At depths of 0.5% and 1%, all three periods were recovered. At 0.1%, only the 1.7-day signal passed. The 3.1-day injection recovered the period but failed the validation rule, and the 5.3-day injection yielded another period. A single background, fixed phases and three depths cannot define survey completeness. The full search is also run with flux sign inverted for each star; real periodic stellar variations can produce detections in both directions. [Injection outcomes and per-target inversion control](../results/tess-pilot.json)

## Reproduce and tile the archive

Use the archive extra for Astropy, NumPy, FITS and Box Least Squares. From a checkout on a machine with a local cache and a connection to MAST:

```sh
python -m pip install -e '.[archive]'
python -m open_space.tess run --plan data/tess-plan.json --cache runs/tess-fits \
  --output runs/tess-pilot.json --workers 4 --injection-controls
python -m open_space.tess run --plan data/tess-expansion-plan.json --cache runs/tess-fits \
  --output runs/tess-expansion.json --workers 4
```

Plan a new page by changing the sector pair and the offset. Keep the seed and inventories fixed to partition a pair without overlaps. Different sector pairs can share targets. Archive inputs are read as data; no downloaded shell scripts run. Source acquisition and CPU analysis are separate steps for large jobs, with safe concurrent publication of matching cached files and checksums validated on replay. Use the discovery job coordinator for resumable per-target work, or use these standalone commands for one batch.

```sh
python -m open_space.tess plan --sectors 1 2 --offset 120 --limit 96 \
  --inventories runs/tess-inventories --output runs/tess-page-unpinned.json
python -m open_space.tess fetch --plan runs/tess-page-unpinned.json --cache runs/tess-fits \
  --workers 4 --output runs/tess-page.json
OPENBLAS_NUM_THREADS=1 python -m open_space.tess run --plan runs/tess-page.json \
  --cache runs/tess-fits --workers 8 --output runs/tess-results.json
```

## What the worker measures

- The worker checks TIC and sector in the FITS header, the catalog URI, time units and time-system reference. It rejects every cadence with a nonzero quality flag, nonpositive uncertainty or nonfinite measurement. It combines samples into weighted 10-minute bins with at least three good cadences. Gaps exceeding 0.35 days split the data. Trend estimates use local medians over ±0.65 days at anchors separated by at most 0.2 days. Short segments and nonpositive trends supply no coverage. Per-bin errors are floored by an empirical successive-difference scatter.
- Weighted Box Least Squares searches 0.75–8 days and 1, 2 and 4 hour durations in the discovery sector. The strongest peak supplies the only ephemeris tested in the later sector. Alternatives separated by a distinct Fourier resolution element are retained, including aliases. These box widths, the trend window and the period range make this particular detector insensitive to many long-period, single-event, shallow or very short transits. [Astropy BLS implementation](https://docs.astropy.org/en/stable/timeseries/bls.html)
- Each covered event needs at least half a duration's expected samples and at least five comparison samples with observations on both sides. The reported depth is the difference between the weighted local baseline and core flux. When individual depths disagree more than expected, the combined uncertainty increases. These signal-to-error ratios use approximate noise estimates for ranking. Correlated stellar noise and the many periods scanned prevent interpreting them as Gaussian discovery significances.
- The two-sector screen requires discovery and validation depth/error ≥7 and ≥5 respectively and ≥3 supported events in both sectors. Deep or inconsistent depths, significant odd/even differences, or a strong half-phase secondary are vetoes. A qualifying box still requires comparison against other periodic shapes, image localization, contamination checks, catalog searches and independent observations. The TESS target position names an aperture; it cannot identify the physical source of a dip on its own. [TESS science data product descriptions](https://heasarc.gsfc.nasa.gov/docs/tess/documentation.html)

The follow-up morphology script fits a fixed centered box, a continuous sinusoid and a four-harmonic periodic shape to the discovery sector. Each model is evaluated at its frozen period and coefficients in the later sector, allowing a scalar zero point for that sector. Compare weighted squared-residual improvements relative to a constant. The Fourier coefficients, harmonic and related periods diagnose broad or multiple-wave variability. The model check is conditional on passing the original dip screen, uses the same held-out data as the triggering calculation and does not provide independent confirmation or a physical classification.

```sh
python -m open_space.tess_morphology --pair data/tess-plan.json runs/tess-pilot.json \
  --pair data/tess-expansion-plan.json runs/tess-expansion.json \
  --cache runs/tess-fits --output runs/tess-shapes.json
```

## What the next 96 stars actually yielded

The extension was chosen after the initial pilot, from target metadata alone, and used the original unchanged thresholds. All 96 were screened. Ninety stayed below the two-sector rule, two reached the signal requirement with vetoes, and four passed the initial photometric rules. Across both exploratory pages the totals are 112 below threshold, four with a veto and four passing out of 120 targets. With the separately identified known-planet control, the two searches used 242 unique FITS products totaling 489,225,600 bytes and 4,356,057 retained two-minute cadences. Three extension stars also passed the same raw screen on *inverted* flux, demonstrating the value of checks on periodic brightness variations. [Extension measurements](../results/tess-expansion.json) · [Predictive shape comparisons](../results/tess-shapes.json)

| Exploratory target | Blind box period (days) | Evidence that changes its interpretation |
| --- | ---: | --- |
| TIC 50380257 | 1.34813 | Exact target and position match in the published TESS-EBs catalog at 1.348285 days, a 0.0115% difference. The 1% deep signal and weaker half-phase feature are a rediscovery of a documented eclipsing-binary signal. A weak secondary below the screen's 20% veto survived that automated threshold. |
| TIC 425975034 | 1.11475 | A third harmonic near 0.3716 d dominates the phase curve. Difference images in both sectors locate the signal roughly 38 arcsec away from the TIC, near the cataloged eclipsing binary BB Tuc. |
| TIC 279511712 | 1.65195 | A sinusoid fit in Sector 1 predicts Sector 2 nearly three times better than a box by weighted improvement over flat. The aperture is centered on the known spectroscopic multiple HD 53921. Attribution among its components requires more evidence. |
| TIC 31780319 | 1.31943 | A single sinusoid at this period accounts for almost the entire improvement delivered by the four-harmonic model on the later sector. The deficit lies near the target in the archival pixels. The photometry supports near-target variability, with its physical mechanism unresolved. |

We queried the frozen MAST TESS-EBs v1.0 catalog for each of the eight exploratory objects that reached the two-sector signal requirement. Its 4,584 rows include one exact counterpart, TIC 50380257. For the other seven there was no match at the TIC or within 90 arcsec in this catalog and epoch convention. A miss here cannot establish novelty: BB Tuc is an instructive counterexample. The public 2018 Planet Hunters sector-1 notes also record the approximately 1.348-day signal of TIC 50380257. The physical and literature classifications are post-screen follow-ups; the underlying blind records retain their original machine disposition. [EB table checks](../results/tess-eb-check.json) · [MAST eclipsing-binary catalog and citation](https://archive.stsci.edu/hlsp/tess-ebs) · [2018 primary observing notes](https://talk.planethunters.org/boards/BPH0000008/discussions/DPH0001l8z/) · [Pixel localization](PIXELS.md)

```sh
python -m open_space.tess_eb_catalog runs/tess-pilot.json runs/tess-expansion.json \
  --cache runs/tess-catalogs --output runs/tess-eb-check.json
```

These stars supply an empirical triage workload for short-period SPOC searches. None currently supports a claim of a new planet. Other period ranges, single-transit searches, diluted events, eclipse timing and less common variable-star shapes remain outside the transit screen.

## Three archival observations nearly two years later

We selected three unresolved variable apertures from the two-sector analysis. Before inspecting later flux, we chose the earliest available two-minute product in Sectors 27–29 for each target, using MAST inventory membership. We independently searched each later sector with both Box Least Squares and a 0.2–8 day, weighted floating-mean sinusoid. We also applied the original box period, duration and phase without refitting. The later-sector tests address recurrence in a selected sample. [Selection and digests](../data/tess-third-plan.json) · [Measurements](../results/tess-third-sector.json)

| TIC | Later sector and elapsed time | Independently strongest sinusoid, early → late | Result at the original box phase |
| --- | --- | --- | --- |
| 234495456 | S28, 736.6 days | 0.380551 → 0.380710 days; 81.8% → 62.6% weighted variance explained | 27 of 27 covered events support a dip; combined depth/error 35.7. The independent box period in both sectors is twice the dominant sinusoidal period. |
| 31780319 | S27, 711.0 days | 1.323085 → 1.320816 days; 78.1% → 86.7% | Depth/error 0.11. Variability recurs on a different camera and CCD, but the original box fails to predict the later phase. |
| 234338669 | S27, 711.0 days | 0.459137 → 0.229690 days; 48.1% → 7.4% | No supported events; depth/error −1.4. A weaker half-period harmonic dominates the later search. The original box ranks third among distinct later peaks. |

These results distinguish repeated frequencies from a persistent ephemeris. A period estimate from one 28-day observing window can accumulate appreciable phase error in two years. Evolving modulation, blending and period aliases remain possible. Three additional files brought the light-curve total to 245 products and 4,403,604 retained two-minute cadences.

```sh
python -m open_space.tess_reobserve plan --tics 234495456 234338669 31780319 \
  --sectors 27 28 29 --inventories runs/tess-inventories --cache runs/tess-fits \
  --output runs/tess-third-plan.json
python -m open_space.tess_reobserve run --later-plan runs/tess-third-plan.json \
  --cache runs/tess-fits --output runs/tess-third-sector.json
```
