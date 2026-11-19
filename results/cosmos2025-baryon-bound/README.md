# COSMOS2025 z > 6 cosmic-baryon-bound test — investigation record

**Claim (provisional, 2026-11-18): No discovery-level violation is established.** This is a *validation status*, **not** a measured null. A positive or statistically meaningful null requires a frozen catalogue extract, mass–redshift posterior propagation, halo volume/extreme-value model, and individual astrophysical vetting. Do not cite this status as an astrophysical constraint.

## Question and prespecified decision
In the public COSMOS-Web COSMOS2025 galaxy catalogue (approx. 0.54 deg² JWST NIRCam/MIRI footprint), is there any galaxy with support at z > 6 for a stellar mass above \(f_b=0.16\) times the most massive halo expected in the survey's appropriate comoving redshift-bin volume? Recompute a Sheth–Tormen or Tinker halo mass function and treat halo maximum as a *distribution*, not a deterministic cutoff. Report both a conservative probability of exceeding a clearly defined percentile halo ceiling and a posterior-predictive population statistic, including selection, Eddington scatter, and a trials correction. A posterior tail above a 1-object expected halo mass is **not** sufficient by itself to claim ΛCDM failure.

**Freeze rule:** fix publication/release, source files and SHA-256, sky footprint, redshift bins, mass conventions/IMF, selection, photometry + flags, priors, model versions, HMF convention and cosmology *before* drawing conclusions. No point-estimate-only candidate cut.

## Ledger (observed, 2026-11-19 UTC)
| Check | Status | Evidence / missing prerequisite |
| --- | --- | --- |
| Scientific label rule | Completed | Read [remarkable findings](../../docs/REMARKABLE_FINDINGS.md); do not add [REMARKABLE] absent discovery-grade validation |
| Public release query + bytes and SHA-256 | **FAILED acquisition (two distinct checks)** | First CDS timed out on [Actions #44445639190](https://github.com/dovvnloading/Open-Space/actions/runs/44445639190), [receipt](run/receipt.json). An independent **local, no-Actions** [official/alternate-host probe](run/host_probe.json) found DNS failures on IAP, Caltech, IRSA, Zenodo, CDS, MAST and a GitHub control. No catalogue rows or bytes verified |
| Posterior-based z>6 candidate selection | Pending | Need joint individual \(p(M_*,z)\) or defensible likelihood-based posterior sampling, not just catalog best fits |
| Halo mass function and survey volume | Nominal geometry checked; HMF provisional and not independently replayed | Numerical job #44445915491 computed a model envelope, but exposed an exponentiation bug since fixed; [exact mathematical repair](HALO_PROVISIONAL.md). No accepted new-run result or survey footprint |
| Eddington bias, look-elsewhere, completeness | Pending | Need measurement likelihood and survey selection |
| Independent SED fits (EAZY + Prospector/Bagpipes) | Pending | Need individual multi-band fluxes, errors, covariances, apertures, priors and calibrated installations |
| AGN/line, dusty z~2–4, brown dwarf, blends | Pending | Need imaging/deblend and competing SED fits |
| MIRI, X-ray/radio, DJA/MAST spectroscopy | Pending | Need source coordinates, crossmatches, archival evidence |
| Independent referee review | Static review complete; signoff blocked | [Referee note](INDEPENDENT_REFEREE.md); arithmetic script committed, not executed against full data |

## Recorded work and replay commands

- **Completed: release inventory and protocol.** [Source register and validation contract](SOURCES_AND_PROTOCOL.md) identify official IAP v1.1 (2026-05-04), the older CDS tabular mirror, precise fixed-width columns and source-vetting gates.
- **Completed: independent nominal survey shell geometry.** [Nominal volumes](nominal_shell_volumes.json), for area=0.54 deg² and specified cosmology, were verified using separate D_c³/3 and ∫ D_c²/H(z) dz formulas. *Not footprint-masked.*
- **Committed:** [full CDS source/hash screen](screen.py), [photometry join](join_photometry.py), [ST99/Tinker08 envelope](halo_bound.py), [joint-posterior diagnostic](posterior_exceedance.py), and [numerical referee](referee.py). **Full catalogue scan failed** because three CDS ReadMe endpoints timed out; no source rows or hashes were acquired. A preliminary HMF run completed with an identified interpolation exponent bug; [receipt and algebraic recovery](HALO_PROVISIONAL.md) are retained for transparency.
- **Evidence status:** [Audit workflow](https://github.com/dovvnloading/Open-Space/actions/workflows/cosmos2025-baryon-bound.yml) generated [a committed failed-download receipt](run/receipt.json). No verified `columns.json`, high-z candidate row extract, photometry join or joint posterior exists. Remaining historical queued jobs are not evidence of success; only a valid full-source scan and independent referee count.

```bash
python results/cosmos2025-baryon-bound/screen.py --out results/cosmos2025-baryon-bound/run
python results/cosmos2025-baryon-bound/join_photometry.py --shortlist results/cosmos2025-baryon-bound/run/highz_triage_top250.csv --out results/cosmos2025-baryon-bound/run
python results/cosmos2025-baryon-bound/halo_bound.py --out results/cosmos2025-baryon-bound/run/halo_bounds.json
python results/cosmos2025-baryon-bound/referee.py --out results/cosmos2025-baryon-bound/run/referee_report.json
# posterior_exceedance.py requires ACTUAL independently fitted paired source posterior draws;
# do not substitute LePHARE marginal intervals.
```

**No catalogue rows or their downloaded SHA-256 hashes, direct post-fix halo quantiles, SED-refit posterior chains, high-z candidate classifications, or spectroscopic crosschecks are currently verified.** Do not treat a failed transfer as a null and do not infer a ΛCDM result from the provisional HMF table.

**Disposition after first attempted replay:** blocked by unavailable source data and source-level astrophysical checks. Draft PR #194 remains unmerged; no root catalog entry or discovery log addition.

## Second access attempt: alternate official routes from local execution container (2026-11-19 UTC)

**Result: BLOCKED AT DNS; no new observational evidence.** Re-run **python results/cosmos2025-baryon-bound/probe_catalog_sources.py --out results/cosmos2025-baryon-bound/run/host_probe.json** locally; no GitHub Actions involved. The [machine-readable result](run/host_probe.json) contains full URLs, exact exceptions and timestamps; the [probe code](probe_catalog_sources.py) runs in parallel and uses bounded HTTPS HEAD only after DNS resolution. Run timestamp: **2026-11-19T03:23:04.970475+00:00**.

- Nine URLs / seven distinct hostnames: official IAP catalogue page and human-access request form; official COSMOS-Web Caltech DR1 portal; IRSA paper bibliography and TAP discovery (neither proven to carry COSMOS2025); Zenodo records search (not a proven mirror); CDS catalogue; MAST bibliography (not a proven catalogue mirror); and GitHub as network control. Hosts: cosmos2025.iap.fr, cosmos.astro.caltech.edu, irsa.ipac.caltech.edu, zenodo.org, cdsarc.cds.unistra.fr, archive.stsci.edu, github.com. **All seven hostnames failed DNS resolution**.
- Every DNS request raised Python **gaierror, errno -3**, exact message **[Errno -3] Temporary failure in name resolution**; every HTTPS call recorded **not_attempted_dns_failed**. A separate local curl -I https://cosmos2025.iap.fr/ returned **curl: (6) Could not resolve host**. Direct TCP connect probes to public test endpoints 1.1.1.1:443, 8.8.8.8:53, 140.82.114.3:443 returned **[Errno 111] Connection refused**. These failures characterize **this local sandbox**, NOT publisher availability.
- The IAP catalogue page advertises **COSMOSWeb_mastercatalog_v1.1.fits**, but its download is reached through a request/captcha form, not an authenticated direct URL that the code can resolve. The catalogue paper explicitly names IAP and CDS; no official IRSA or Zenodo mirror was demonstrated.
- Local file SHA-256 immediately after execution: [probe code](probe_catalog_sources.py) **a11c743d6eb92e111e7fcd9f6b1764d6855e10f1fe55bde2169e2ccd7e00e5d0**; [output receipt](run/host_probe.json) **24780cda971f0fb8d48cd9742bf1e1eb7ec1259924be73ca992854f41bc7555f**. These are hashes of the probe and failure receipt, **not the missing catalogue**.

**Scientific stop condition reached:** A frozen catalogue could not be acquired here. No verified z>6 mass-posterior selection, independent SED fitting, per-object spectroscopy/interloper vetting, selection/Eddington corrections or exceedance probability was performed. The numerical halo envelope remains provisional. There is no measured null, no robust violation and no discovery label. Keep PR #194 draft and unmerged; do **not** add a root README result.

## Third access attempt: one dedicated path-filtered GitHub runner (2026-11-19 UTC)

[Actions run #44454134461](https://github.com/dovvnloading/Open-Space/actions/runs/44454134461) ran **one job** on an internet-connected GitHub runner, not the local sandbox. It recomputed the complete nominal-area six-bin Sheth–Tormen/Tinker halo envelope with current Colossus and passed the new numerical referee's independent shell-volume integration. This model scale still is **not a survey upper limit** or a COSMOS2025 observation.

The [complete machine-readable alternate-source receipt](run/noncds/mirror_receipt.json) records **nine exact routes** and responses: IAP catalogue and download-form URLs resolved in DNS but timed out at HTTPS; Caltech DR1 and the IRSA COSMOS overview both returned HTTP 200, but linked no verified full 784,016-row COSMOS2025 FITS; IRSA TAP tables returned HTTP 200 with no COSMOS2025 table; the paper's A&A URL returned HTTP 403; GitHub connectivity control returned HTTP 200. Both Zenodo search URLs returned **HTTP 400 because our original anonymous request asked for 30 results, above Zenodo's 25-result limit**. This client defect is recorded, **not passed off as a Zenodo outage or absence of a mirror**; the script has been corrected to 25 and the same single-job workflow is path-triggered again to retry discovery. [Zenodo's published limit](https://blog.zenodo.org/2025/11/25/2025-11-14-search-api-updates/).

The runner committed compact [source/host receipt](run/noncds/mirror_receipt.json), [full numerical halo output](run/noncds/halo_bounds.json), [version pins](run/noncds/numerical_versions.txt) and [independent referee output](run/noncds/independent_referee.json) without bulk FITS. Original catalogue bytes acquired: **zero**. [Referee's file-hash audit](run/noncds/independent_referee.json) verified the source receipt SHA-256 `2c85ec412bed846face13696cd0035af610ff6572e17503661781bf9ec2beecc` and halo output SHA-256 `6640731301df680bd862b5b7b83f5f0d5d21e4449e4a7b0ddd0fad4413712e8e`. These are hashes of **runner evidence**, not of the absent catalogue.

**Decision:** scientific test **inconclusive due to unacquired source data**; no joint posterior, no candidate refits, no per-object contaminant/NIRSpec checks, and no exceedance probability. Keep PR #194 **draft, unmerged**; no root catalogue finding and no [REMARKABLE] entry. The corrected Zenodo probe must not be called successful unless its own new committed receipt demonstrates real catalogue bytes and valid FITS rows.

## Final alternate-mirror outcome: corrected Zenodo run and locally replayed referee (2026-11-19 UTC)

The corrected one-job [GitHub Actions run #44454925847](https://github.com/dovvnloading/Open-Space/actions/runs/44454925847) **ran to completion and committed its full [mirror receipt](run/noncds/mirror_receipt.json)**. The anonymous Zenodo API queries now returned **HTTP 200**: 8 COSMOS2025-related records examined, 0 COSMOSWeb_mastercatalog records, and **no team-authored full 784,016-row mastercatalog linked**. IAP DNS resolved but the catalogue and download form both timed out. Caltech DR1, IRSA COSMOS overview and IRSA TAP returned HTTP 200; the IRSA TAP tables listed no COSMOS2025 table, and no direct full-release FITS was found. An A&A HTML endpoint returned HTTP 403. **Zero source bytes and zero verified rows**, so no catalogue SHA-256 exists. Do not confuse this controlled search with proof that no mirror exists anywhere.

The [new source receipt](run/noncds/mirror_receipt.json) has SHA-256 **472c644d622b785b6ab6d4780dab2e41f9cd590cb182a11111116bcba62a78df**. The freshly recomputed, entirely *model-only* [halo output](run/noncds/halo_bounds.json) has SHA-256 **d2ddfa4cba1250737e54b54ffc6f5b0d51dc5f66a3a409dffa9701285d5dc8f5**. These refer to the *second* run and supersede the earlier run's hashes. The released runner [ZIP artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44454925847) has SHA-256 **f4d517693bb5d679ee92a14d7d487d5407cb366bd3984352caa33b4338337e1b**.

**Independent local replay:** [audit script](audit/referee_runner_artifact.py), [machine summary](run/noncds/local_referee_audit.json). I unpacked the real GitHub artifact independently of the CI code; verified all source/hashing claims and recomputed both the six-shell comoving-volume integrals and 48 Poisson maximum-halo quantiles (two HMFs × six bins × four probabilities). Maximum relative differences were **8.01×10⁻⁶** in volumes and **9.11×10⁻¹⁵** in Poisson-tail interpolation. *Software defect recorded*: the second CI machine-referee script rejected a newly recomputed halo output against a **stale first-run SHA-256** stored in the overwritten mirror receipt, and the old referee JSON remained on disk. Its stale "passed" report does **not** validate run #44454925847. The independent local calculation above checks this run's actual output instead. The script's rollover hash logic must be repaired before a fresh successful run can claim machine-level referee signoff.

**Final claim: INCONCLUSIVE ACQUISITION FAILURE — NOT A MEASURED NULL, NOT A ΛCDM VIOLATION.** None of the candidate-posterior, two-code SED, AGN/dwarf/dust/blend, MIRI/X-ray/radio/NIRSpec or Eddington/selection/trials checks can be performed without the catalogue. **Keep PR #194 as a draft, do not merge, do not add a root result entry or [REMARKABLE] designation.**
