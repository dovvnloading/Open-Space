# COSMOS2025 z > 6 cosmic-baryon-bound test — investigation record

**Claim (provisional, 2026-11-18): No discovery-level violation is established.** This is a *validation status*, **not** a measured null. A positive or statistically meaningful null requires a frozen catalogue extract, mass–redshift posterior propagation, halo volume/extreme-value model, and individual astrophysical vetting. Do not cite this status as an astrophysical constraint.

## Question and prespecified decision
In the public COSMOS-Web COSMOS2025 galaxy catalogue (approx. 0.54 deg² JWST NIRCam/MIRI footprint), is there any galaxy with support at z > 6 for a stellar mass above \(f_b=0.16\) times the most massive halo expected in the survey's appropriate comoving redshift-bin volume? Recompute a Sheth–Tormen or Tinker halo mass function and treat halo maximum as a *distribution*, not a deterministic cutoff. Report both a conservative probability of exceeding a clearly defined percentile halo ceiling and a posterior-predictive population statistic, including selection, Eddington scatter, and a trials correction. A posterior tail above a 1-object expected halo mass is **not** sufficient by itself to claim ΛCDM failure.

**Freeze rule:** fix publication/release, source files and SHA-256, sky footprint, redshift bins, mass conventions/IMF, selection, photometry + flags, priors, model versions, HMF convention and cosmology *before* drawing conclusions. No point-estimate-only candidate cut.

## Ledger (initial)
| Check | Status | Evidence / missing prerequisite |
| --- | --- | --- |
| Scientific label rule | Completed | Read [remarkable findings](../../docs/REMARKABLE_FINDINGS.md); do not add [REMARKABLE] absent discovery-grade validation |
| Public release query + bytes and SHA-256 | Pending | Freeze canonical COSMOS2025 release/table, machine query, download time, row count and exact hashes |
| Posterior-based z>6 candidate selection | Pending | Need joint individual \(p(M_*,z)\) or defensible likelihood-based posterior sampling, not just catalog best fits |
| Halo mass function and survey volume | Pending | Need cosmology and footprint/mask, independent ST/Tinker check |
| Eddington bias, look-elsewhere, completeness | Pending | Need measurement likelihood and survey selection |
| Independent SED fits (EAZY + Prospector/Bagpipes) | Pending | Need individual multi-band fluxes, errors, covariances, apertures, priors and calibrated installations |
| AGN/line, dusty z~2–4, brown dwarf, blends | Pending | Need imaging/deblend and competing SED fits |
| MIRI, X-ray/radio, DJA/MAST spectroscopy | Pending | Need source coordinates, crossmatches, archival evidence |
| Independent referee review | Pending | Review outcome to be recorded after substantive tests |

## Reproducibility
Scripts, receipts, derived machine-readable summaries and exact pending versus completed checks will be added as evidence is obtained. Never manufacture masses, catalogue rows, hashes, photometry, posterior samples, SED refits or spectroscopy.

**Disposition at creation:** open investigation; not a detection; not a quantitative null.
