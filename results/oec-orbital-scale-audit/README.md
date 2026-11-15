# Independent orbital-scale audit: a repeated solar-radius error and three unrelated field transcriptions

**Research contribution | 15 November 2026 | Source-record corrections, not a claim to discover new planets.**

This is a new investigation, not an extension of an existing Open-Space target campaign. The repo's open pull requests and 227 `results/` entries were checked before selecting it. The [preregistered plan](PLAN.md) predates these follow-up conclusions. It seeks falsifiable inconsistencies in *reported orbital and stellar properties* and separates anomalies inherited from the underlying literature from mistakes in this particular catalogue.

## Main result: six planets in two systems carry one identical unit-conversion error

The **live, 15 November 2026 Open Exoplanet Catalogue (OEC) XML** still reports the following semimajor axes. The primary discovery paper gives the transit-fit `a/R_star` values. Multiplying those *dimensionless* values by the **nominal solar radius** (0.00465046726 AU) reproduces the six published OEC axes to within 1.5%. But the hosts are only **0.21 and 0.24 solar radii**. Multiplying by the *actual quoted stellar radius* instead yields the following central-value estimates:

| Planet | Period (days) | OEC `a` (AU) | Original `a/Rstar` | Star radius (Rsun) | From original paper (AU) | OEC axis too large by |
|---|---:|---:|---:|---:|---:|---:|
| Kepler-445 b | 2.984151 | 0.102 | 21.94 | 0.21 | 0.02143 | 4.76× |
| Kepler-445 c | 4.871229 | 0.141 | 30.21 | 0.21 | 0.02950 | 4.78× |
| Kepler-445 d | 8.152750 | 0.200 | 42.95 | 0.21 | 0.04194 | 4.77× |
| Kepler-446 b | 1.565409 | 0.067 | 14.20 | 0.24 | 0.01585 | 4.23× |
| Kepler-446 c | 3.036179 | 0.104 | 22.40 | 0.24 | 0.02500 | 4.16× |
| Kepler-446 d | 5.148921 | 0.147 | 31.60 | 0.24 | 0.03527 | 4.17× |

This is more specific than just identifying a mismatch with Kepler's third law: it **identifies a common transformation** that reproduces the erroneous values across both unrelated stars and all six planets. The paper's independent stellar temperatures and incident-flux estimates provide a falsification test:

| Planet | Incident flux in Muirhead et al. 2015 (Earth flux) | Recomputed with corrected `a` | Miscomputed with OEC `a` |
|---|---:|---:|---:|
| Kepler-445 b | 8.57 | 8.60 | 0.379 |
| Kepler-445 c | 4.52 | 4.53 | 0.199 |
| Kepler-445 d | 2.24 | 2.24 | 0.099 |
| Kepler-446 b | 26.22 | 26.30 | 1.472 |
| Kepler-446 c | 10.54 | 10.57 | 0.611 |
| Kepler-446 d | 5.30 | 5.31 | 0.306 |

Insolation scales as `1/a²`. The reported OEC axes would reduce incident flux by roughly **17–23 times** versus the source publication, substantially affecting habitability or irradiation-bin studies using these fields. The **same primary paper** (Muirhead et al. 2015, *ApJ* 801, 18, [DOI](https://doi.org/10.1088/0004-637X/801/1/18), Tables 3 and 4 and stellar properties) supplies both transit `a/Rstar` and incident flux; the irradiation calculation from the corrected axes is an additional, nonorbital cross-check, not independent telescope data. A separate [NASA Exoplanet Archive system record](https://exoplanetarchive.ipac.caltech.edu/overview/Kepler-445) confirms the source-model parameters and presents other, distinct later models. Recomputations and diagnostics are machine readable in `primary_checks.csv` and the checked-in `unit_conversion_checks.csv`.

**Important interpretation:** We do not declare the exact current best-fit orbital axes. The quoted stellar radii have uncertainties and the NASA archive contains later models. These are corrected values **on the OEC's adopted original 2015 stellar/transit solution**, not replacements for a new joint fit. Their uncertainty bars must be propagated in any upstream correction; the OEC's present `a` uncertainty must not simply be retained unchanged.

## Three further corroborated field errors

These are different error mechanisms, in the **same live 2026 upstream OEC XML**, rather than artifacts of the stale 2023 export.

| Object / field | Live OEC value | Correct same-solution value | Evidence and competing explanation |
|---|---:|---:|---|
| **HD 39392 b** / `semimajoraxis` | 0.11 AU | **1.08 ± 0.03 AU** | Source star mass 1.08 Msun, orbit 394.3 d; Wilson et al. 2016, [A&A DOI](https://doi.org/10.1051/0004-6361/201527581), and [exoplanet.eu model provenance](https://exoplanet.eu/catalog/hd_39392_b--2566/) list 1.08 AU. Same OEC ±0.03 AU error bar is attached to the wrong, one-decimal-place-short value. |
| **Kepler-441** / host `mass` | 0.0572 Msun | **0.572 (+0.049/−0.053) Msun** | [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/overview/Kepler-441) quotes the Torres et al. 2015 stellar solution at exactly 0.572 with the **same asymmetric uncertainties** as OEC (which erroneously moved only the central decimal). Other published stellar solutions agree in order of magnitude. This is a host-star correction, not a tenth planet. |
| **WASP-34 c** / `semimajoraxis` | 14.96 AU | **5.05 (+0.65/−0.46) AU** | The OEC has `mass*sin(i) = 14.96 Mjup`, **identical to its erroneous central `a = 14.96 AU`**; the older 4093-day circular-orbit solution at [exoplanet.eu](https://exoplanet.eu/catalog/wasp_34_c--10100/) reports **5.05 AU** with the **same a error bars** as OEC. This is a probable wrong-field copy, not a new companion. The 5.05 AU value describes the **historical** solution ([Knutson et al. 2014](https://doi.org/10.1088/0004-637X/785/2/126)), not the latest fitted orbit. As of the 2026-11-15 review, [exoplanet.eu](https://exoplanet.eu/catalog/wasp_34_c--10100/) displays a newer default value of **3.158 AU** (updated 2026-10-20) from a different model; the old model must not be silently substituted for that update. |

All nine field-level discrepancies across **five systems** were checked directly in the pinned source XML at `f4228675cfde38cdadc2f3dc93688edb11d89f77` (committed 2026-11-15 07:01 UTC). See `validated_cases.csv` for one row per correction, original XML links, independent comparison and qualification. Their prior identification by other researchers, if any, has **not** been established; publication priority is unknown.

## Separate discovery: inconsistent source *publication* values propagated to many catalogues

The 2019 *A&A* study of **K2-285** ([Palle et al. 2019](https://doi.org/10.1051/0004-6361/201834001)) reports four periods and semimajor axes for planets around one star of mass **0.830 ± 0.023 solar masses**. The values appear unchanged in OEC, [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/overview/K2-285) and [TEPCat](https://www.astro.keele.ac.uk/jkt/tepcat/planets/K2-285e.html). Yet taking the given orbital axes as *Keplerian* and the listed periods at face value implies **four different central dynamical host masses**:

| Planet | Period (days) | Published `a` (AU) | Mass implied by `a³/P²` (Msun) | Factor vs host mass 0.830 |
|---|---:|---:|---:|---:|
| K2-285 b | 3.471745 | 0.03817 | 0.616 | 0.742 |
| K2-285 c | 7.138048 | 0.08240 | 1.465 | 1.765 |
| K2-285 d | 10.45582 | 0.11780 | 1.995 | 2.403 |
| K2-285 e | 14.76289 | 0.18041 | 3.594 | 4.330 |

The relative mismatch is particularly large for e and cannot be explained by the reported small planetary masses. All four reported `a` values nevertheless agree closely with the paper's **separately fitted transit** `a/Rstar` ratios multiplied by its stellar radius **0.787 Rsun**. Thus this is **not** an OEC conversion mistake: it is an **internal tension between published transit-fit geometry and the stated shared stellar mass/periods**, propagated by downstream catalogues. It requires a review of transit fit, circularity, transit-duration systematics and any priors or parameter definitions before any physical explanation or replacement axis is claimed. Report as an **unresolved source-model consistency finding**, not four confirmed transcription errors. Exoplanet multiplicity, eccentricity and transit-fit systematics must be investigated. This contrasts with the strongly diagnosed 445/446 unit mistake. Inputs and numeric checks are in `k2_285_primary_checks.csv` and `source_model_checks.csv` (the latter produced by the same standard-library replay).

## Frozen screening denominator, and what 70 flags do *not* mean

A deterministic search over a **2023-10-15 OEC generated CSV**, not the live 2026 master, yielded:

| Quantity | Number |
|---|---:|
| CSV records inspected | **5,414** |
| Confirmed, nonbinary planet records with finite positive orbital period, axis, host mass | **2,275** |
| `Q > 1.5` or `Q < 1/1.5` | **70** |
| Factor disagreement > 2, > 3, > 5, > 10 | **38**, **24**, **15**, **10** |

`Q = (a/AU)^3 / [(P/365.25 days)^2 × (Mstar/Msun)]`, and the associated mass ratio is dimensionless. The solar mass of a massive companion is neglected in *screening* only (a ~14 Mjup planet contributes ~1.3%, not factors of 2–90). A single host's planets must share a dynamical host mass for consistent semimajor axes under a common two-body convention. A 1-year, 1-AU orbit around a 1-solar-mass star is a positive numerical control. **Seventy is not the number of catalogue errors**: source fits, alternate catalog releases, eccentric/transit geometry, unresolved binaries, planetary mass, imprecise or broad error bars, and stale source properties can flag valid rows. Not every flagged host has been independently validated. See generated `flags.csv` for the **complete list** and `summary.json` for the denominator and each exclusion.

### Reproduce (standard-library only, no live dependencies)

1. Download the frozen [OEC 2023 CSV](https://raw.githubusercontent.com/OpenExoplanetCatalogue/oec_tables/bf1efad9e4bdc9f0887fe392b77308c3585e263f/comma_separated/open_exoplanet_catalogue.txt) and keep those exact bytes. Git commit `bf1efad9e4bdc9f0887fe392b77308c3585e263f`, Git blob `db832c9d70e953ca6e1851122e1d46578c3b5d58`.
2. From repo root run:

   ```bash
   python results/oec-orbital-scale-audit/audit.py \
     --input /path/to/open_exoplanet_catalogue.txt \
     --primary-checks results/oec-orbital-scale-audit/primary_checks.csv \
     --source-model-checks results/oec-orbital-scale-audit/k2_285_primary_checks.csv \
     --out /tmp/oec-orbit-replay
   python -m unittest discover -s tests -p 'test_oec_orbital_scale_audit.py' -v
   ```

3. Inspect `/tmp/oec-orbit-replay/flags.csv`, `summary.json`, and `unit_conversion_checks.csv`. The audit **verifies the input's Git blob SHA-1 against the pinned 2023 snapshot by default** and fails closed on different bytes. It also reports the actual SHA-256 of its input and refuses an explicit `--expect-sha256` mismatch. To analyze fixtures or newer archives, pass `--allow-unpinned-input`; the resulting JSON then marks `pinned_input_verified: false` and does **not** attribute the input to the pinned OEC commit. It does not download live data or mutate existing campaigns. Its separate primary check uses independently published transit ratios, stellar temperatures and irradiation rather than fitting the orbital audit metric to the answer. **Never** apply the curated correction values across a different stellar-model epoch without checking compatibility.

### Evidence provenance and limits

- **Live source XML (pinned):** [`Kepler-445`](https://raw.githubusercontent.com/OpenExoplanetCatalogue/open_exoplanet_catalogue/f4228675cfde38cdadc2f3dc93688edb11d89f77/systems/Kepler-445.xml), [`Kepler-446`](https://raw.githubusercontent.com/OpenExoplanetCatalogue/open_exoplanet_catalogue/f4228675cfde38cdadc2f3dc93688edb11d89f77/systems/Kepler-446.xml), [`HD 39392`](https://raw.githubusercontent.com/OpenExoplanetCatalogue/open_exoplanet_catalogue/f4228675cfde38cdadc2f3dc93688edb11d89f77/systems/HD%2039392.xml), [`Kepler-441`](https://raw.githubusercontent.com/OpenExoplanetCatalogue/open_exoplanet_catalogue/f4228675cfde38cdadc2f3dc93688edb11d89f77/systems/Kepler-441.xml), [`WASP-34`](https://raw.githubusercontent.com/OpenExoplanetCatalogue/open_exoplanet_catalogue/f4228675cfde38cdadc2f3dc93688edb11d89f77/systems/WASP-34.xml), and [`K2-285`](https://raw.githubusercontent.com/OpenExoplanetCatalogue/open_exoplanet_catalogue/f4228675cfde38cdadc2f3dc93688edb11d89f77/systems/K2-285.xml). Git SHA prevents a nightly XML update from silently changing this result.
- **Independent original literature:** Muirhead et al. 2015 ApJ 801,18 DOI `10.1088/0004-637X/801/1/18`, Wilson et al. 2016 A&A 588,A144 DOI `10.1051/0004-6361/201527581`, Torres et al. 2015 stellar solutions via [NASA Archive](https://exoplanetarchive.ipac.caltech.edu/overview/Kepler-441), old Knutson et al. 2014 WASP-34 c 4093-day solution via [exoplanet.eu](https://exoplanet.eu/catalog/wasp_34_c--10100/), Palle et al. 2019 A&A 623,A41 DOI `10.1051/0004-6361/201834001`.
- **Review limits:** These are archival measurements and logical checks, not independent new photometry or a publication-priority search. We have not submitted a correction to the upstream OEC maintainers or refit the original data. The 2023 source table is stale relative to the 2026 live XML, and source-model changes can invalidate naive comparisons. Good next tests are to patch copies of the five affected XML systems, ask maintainers to assess original uncertainty propagation, and request an author-side review of K2-285 fit conventions. A full 2026 XML all-catalog audit remains unperformed.
