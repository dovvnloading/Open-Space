# GWTC lensing census: runner result

**Disposition: no validated strong-lensing discovery.** This is a *partial* posterior search plus an all-pair summary pre-screen; no published detection claim.

- Eligible events by version-pinned catalogue: {'GWTC-2.1': 34, 'GWTC-3': 23, 'GWTC-4.0': 84}.
- All-pair mass/chi_eff/time-delay **proxy**: 9870 pairs. Approximate SIS population, actual GWOSC event cadence, 1,500 catalogue-stratified shuffles.
- Proxy best: GW190720_000836 / GW190728_064510; proxy family-wise null-tail fraction = 0.6376. **Not** a posterior-overlap false-alarm probability.
- Downloaded and MD5-validated 10 distinct released PE files; 0 unavailable. 5 preselected pairs screened with full PE samples. Download SHA256 and sampling labels in posterior_receipts.json.
- LVK O4a control: GW231113_200417 / GW231223_202619; published joint log10 Bayes factor -3.3, **not lensed**. Earlier O3a controls likewise negative joint odds.
- No phase-coherent Morse image-type likelihood, joint GWOSC strain fit or foreground galaxy crossmatch: **pending**, not 'passed'.
- A mass in the **source frame** depends on luminosity distance/redshift and so is not a lensing invariant; detector-frame masses are the gate.

## Screened posterior pairs

| Event 1 | Event 2 | log10 approximate overlap B | status |
| --- | --- | ---: | --- |
| GW231113_200417 | GW231223_202619 | 5.94 | not discovery |
| GW190720_000836 | GW190728_064510 | 4.49 | not discovery |
| GW191103_012549 | GW191204_171526 | 3.42 | not discovery |
| GW190421_213856 | GW190910_112807 | 3.33 | not discovery |
| GW191129_134029 | GW191216_213338 | -57.40 | not discovery |

## Constraints and receipts

See result.json, all_pairs_proxy.tsv.gz, proxy_background.json, time_delay_population.json, gwosc_event_freeze.json, release_manifest.json, posterior_pairs.json, posterior_receipts.json and download_failures.json (if any).
Posterior overlap uses symmetric KDE ratios to an approximate common PE prior, and sky HEALPix overlap; no lens selection bias or full-population PE FAP has been evaluated. The summary pre-screen cannot supply sky localization, image phase or foreground lens identification.
The controls' joint Bayes factors are published **external** results, not our recalculated joint evidence; a large overlap B cannot override their negative joint-analysis conclusions.
No [REMARKABLE] label per docs/REMARKABLE_FINDINGS.md.
