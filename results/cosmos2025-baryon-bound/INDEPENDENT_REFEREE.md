# Independent-referee pass — provisional, code-level only

**Disposition: BLOCKED.** Do not call this a measured null, a ΛCDM violation, a discovery or a vetted candidate. **Keep PR #194 draft; no [REMARKABLE] label or log entry.** This is a separate methodological review of the deposited specification and geometry, **not** an external astronomer's replication and not a signoff on unacquired data.

## Independently supported observations

- Nominal 0.54 deg² shell volumes checked using both analytical flat-space D_c³/3 and direct quadrature of dV/dz. Values agree to printed precision. For (6,7], the nominal volume is **4,717,016.070 Mpc³**; the five other values are in [nominal_shell_volumes.json](nominal_shell_volumes.json). These are *not* depth/footprint-mask volumes.
- CDS metadata specify **784,016** LePHARE and photometry rows and the fixed-width offsets used in the streaming code. Official IAP release **v1.1** is newer than the CDS revision and corrects MIRI-associated photometry in B5/B9/B10. Releases have not been demonstrated byte-equivalent.
- Shuntov et al. (A&A 704 A339), Sect. 7.1, report that their mass physical-parameter errors are obtained **at fixed photometric redshift**, without propagating its uncertainty. It is scientifically invalid to recover joint p(z, Mstar) by independently sampling separate point estimates or 68% intervals.
- The submitted scripts define a Poisson-distributed *maximum* halo rather than calling an expected-one threshold a hard law. Colossus Sheth99 FOF masses versus Tinker08 M200m masses require conversion or separate reporting.

## Blocking defects and missing evidence

1. **Catalogue acquisition failed.** [Actions #44445639190](https://github.com/dovvnloading/Open-Space/actions/runs/44445639190) timed out at the ReadMe stage on all three CDS source endpoints. The exact [committed failure receipt](run/receipt.json) is inspectable; no actual catalogue bytes, SHA256, full-row replay, source ID or measured flux has been obtained. This failure is not an observed null.
2. **No joint mass-redshift selection.** The exploratory high-z cut based on zfinal or the upper 68% photo-z bound cannot find every secondary high-z PDF mode. Galaxy/AGN Type 0/2 separation is a triage, and stars misclassified in the source or prior fit remain adversarial controls.
3. **No two independent SED posterior runs.** Fit EAZY plus Prospector/Bagpipes from calibrated multi-band source data; include AGN, high-EW nebular lines, dust/SFH/IMF redshift alternatives, covariance, model-prior sensitivity and sampler diagnostics.
4. **No individual source-vetting ledger.** Brown dwarfs, dusty z~2–4 interlopers, extreme lines, unresolved blends, little red dots, MIRI coverage and limits, Chandra and radio matches, and cuts from the actual images remain untested.
5. **No DJA/MAST/literature NIRSpec crossmatch.** For each actual coordinate, verify exposure, reduction, line ID and significance; no catalogue match is not proof of no spectrum.
6. **No Eddington/selection/trials-aware significance.** Survey depth masks, halo cosmic variance, lensing, mass systematic uncertainties and source/bin look-elsewhere effects must enter a posterior-predictive population comparison. An expected-one bound crossing is not a ΛCDM rejection.
7. **Provisional numerical HMF is not referee-cleared.** [Actions #44445915491](https://github.com/dovvnloading/Open-Space/actions/runs/44445915491) executed both HMF integrations, but the first threshold inversion used base-10 power for natural logs, producing unphysical max-halo masses. The exact [arithmetic repair receipt](HALO_PROVISIONAL.md) documents the mathematically recoverable values, the original SHA256 and corrective commits. Re-run the fixed code with independent Poisson tail tests before accepting them. Tinker08 z>6 extrapolation and FOF versus M200m mass definitions are additional systematics.

## Severity and software defect note

**Referee found and required correction of a blocking numerical bug,** and made the HMF code enforce mass-grid/monotonicity and Poisson-tail checks. A failed archive transfer and a preliminary numerical correction do not demonstrate source-level astrophysics. The manuscript claim remains unsupported.

## Decision

**DO NOT MERGE YET.** Do not mark ready, add a root README result, or alter docs/REMARKABLE_FINDINGS.md until the work yields a completed scientifically supported detection, quantitative null, or rigorously characterized inconclusive outcome. An unfinished data download is not a null result. Another referee should reproduce the frozen code/data hashes and independently sign off on all astrophysical vetting before escalation.
