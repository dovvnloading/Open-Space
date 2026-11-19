# COSMOS2025 z > 6 cosmic-baryon-bound test — investigation record

**Claim (provisional, 2026-11-18): No discovery-level violation is established.** This is a *validation status*, **not** a measured null. A positive or statistically meaningful null requires a frozen catalogue extract, mass–redshift posterior propagation, halo volume/extreme-value model, and individual astrophysical vetting. Do not cite this status as an astrophysical constraint.

## Question and prespecified decision
In the public COSMOS-Web COSMOS2025 galaxy catalogue (approx. 0.54 deg² JWST NIRCam/MIRI footprint), is there any galaxy with support at z > 6 for a stellar mass above \(f_b=0.16\) times the most massive halo expected in the survey's appropriate comoving redshift-bin volume? Recompute a Sheth–Tormen or Tinker halo mass function and treat halo maximum as a *distribution*, not a deterministic cutoff. Report both a conservative probability of exceeding a clearly defined percentile halo ceiling and a posterior-predictive population statistic, including selection, Eddington scatter, and a trials correction. A posterior tail above a 1-object expected halo mass is **not** sufficient by itself to claim ΛCDM failure.

**Freeze rule:** fix publication/release, source files and SHA-256, sky footprint, redshift bins, mass conventions/IMF, selection, photometry + flags, priors, model versions, HMF convention and cosmology *before* drawing conclusions. No point-estimate-only candidate cut.

## Ledger (initial)
| Check | Status | Evidence / missing prerequisite |
| --- | --- | --- |
| Scientific label rule | Completed | Read [remarkable findings](../../docs/REMARKABLE_FINDINGS.md); do not add [REMARKABLE] absent discovery-grade validation |
| Public release query + bytes and SHA-256 | Pending execution | Freeze canonical COSMOS2025 release/table, machine query, download time, row count and exact hashes |
| Posterior-based z>6 candidate selection | Pending | Need joint individual \(p(M_*,z)\) or defensible likelihood-based posterior sampling, not just catalog best fits |
| Halo mass function and survey volume | Nominal geometry checked; HMF pending | Geometry in `nominal_shell_volumes.json`; two-HMF script committed, selection footprint/cosmic variance pending |
| Eddington bias, look-elsewhere, completeness | Pending | Need measurement likelihood and survey selection |
| Independent SED fits (EAZY + Prospector/Bagpipes) | Pending | Need individual multi-band fluxes, errors, covariances, apertures, priors and calibrated installations |
| AGN/line, dusty z~2–4, brown dwarf, blends | Pending | Need imaging/deblend and competing SED fits |
| MIRI, X-ray/radio, DJA/MAST spectroscopy | Pending | Need source coordinates, crossmatches, archival evidence |
| Independent referee review | Static review complete; signoff blocked | [Referee note](INDEPENDENT_REFEREE.md); arithmetic script committed, not executed against full data |

## Recorded work and replay commands

- **Completed: release inventory and protocol.** [Source register and validation contract](SOURCES_AND_PROTOCOL.md) identify official IAP v1.1 (2026-05-04), the older CDS tabular mirror, precise fixed-width columns and source-vetting gates.
- **Completed: independent nominal survey shell geometry.** [Nominal volumes](nominal_shell_volumes.json), for area=0.54 deg² and specified cosmology, were verified using separate D_c³/3 and ∫ D_c²/H(z) dz formulas. *Not footprint-masked.*
- **Committed, not yet executed successfully:** [complete 784,016-row CDS source/hash screen](screen.py); [all marginal high-z row indices and top-250 triage](screen.py); [full photometry source hash + coordinate/flux row join](join_photometry.py); [volume-integrated two-HMF envelope](halo_bound.py); [joint SED-draw posterior diagnostic](posterior_exceedance.py); [independent numerical referee](referee.py).
- **Queued:** [GitHub Actions research run](https://github.com/dovvnloading/Open-Space/actions/workflows/cosmos2025-baryon-bound.yml) should commit run/receipt.json, run/columns.json, run/highz_triage_all.csv.gz, run/highz_triage_top250.csv, run/photometry_receipt.json, run/highz_candidate_photometry.csv, run/halo_bounds.json, run/referee_report.json, and exact numerical package versions. **These files do not establish evidence until the job completes with valid source hashes, full row counts, and independent checks.**

```bash
python results/cosmos2025-baryon-bound/screen.py --out results/cosmos2025-baryon-bound/run
python results/cosmos2025-baryon-bound/join_photometry.py --shortlist results/cosmos2025-baryon-bound/run/highz_triage_top250.csv --out results/cosmos2025-baryon-bound/run
python results/cosmos2025-baryon-bound/halo_bound.py --out results/cosmos2025-baryon-bound/run/halo_bounds.json
python results/cosmos2025-baryon-bound/referee.py --out results/cosmos2025-baryon-bound/run/referee_report.json
# posterior_exceedance.py requires ACTUAL independently fitted paired source posterior draws;
# do not substitute LePHARE marginal intervals.
```

**No catalogue rows, their downloaded SHA-256 hashes, halo-model quantiles, SED-refit posterior chains, high-z candidate classifications, or spectroscopic crosschecks are currently verified. Never manufacture them.**

**Disposition at creation:** open investigation; not a detection; not a quantitative null.
