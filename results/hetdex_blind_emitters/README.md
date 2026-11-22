# Rare emission lines in the two public HETDEX source catalogs

**No new emitter passed all confirmation checks** in 1,340,413 source observations. The blind scan found 541,205 fitted peaks at signal-to-noise of at least 7. The follow-up ledger records external-check attempts for 287 line or spatial hypotheses. Across all stages, 12,306 hypotheses remained unresolved. The releases share targets and exposures. Counts describe observations or feature hypotheses and can include repeated measurements of a galaxy.

## Search coverage

| Release | Source observations scanned | Distinct release source IDs | Fitted peaks | Raw shuffled spectra | Off-source apertures | Same-position stacks from distinct shots |
|---|---:|---:|---:|---:|---:|---:|
| SC1 | 232,650 | 223,641 | 141,783 | 28,928 | 588 | 7,039 |
| SC2 | 1,107,763 | 1,107,763 | 399,422 | 138,496 | 482 | 18,791 |

The source search uses published PSF-weighted extractions and catalog emission knots. Diffuse sources absent from these inputs are outside the candidate search. Off-source cube apertures serve as controls. Each public 1D spectrum combines fibres and dithers from one shot. A map from that shot tests spatial morphology.

## Leads and dispositions

| Channel | Single-visit screen, SC1 / SC2 | Excluded | Already published | Unresolved | Passed | External follow-up, SC1 / SC2 |
|---|---:|---:|---:|---:|---:|---:|
| He II | 1689 / 5098 | 6,926 | 4 | 183 | 0 | 32 / 32 |
| Extreme Lyα EW | 34 / 31 | 47 | 4 | 379 | 0 | 50 / 43 |
| Extended Lyα | spatial / spatial | 103,827 | 677 | 4,725 | 0 | 34 / 32 |
| Unassigned-line screen | 1879 / 5791 | 991 | 0 | 7,019 | 0 | 32 / 32 |

Spatial proposals, subthreshold helium–hydrogen line pairs, imaging-led EW measurements and independent-visit stacks enter after the single-visit screen. A hypothesis can appear in more than one channel. Ranked follow-up gives priority to independent visits or DESI counterparts. Rows outside external follow-up retain an unresolved disposition.

An isolated line anywhere on this wavelength grid admits a Lyα redshift. The unassigned-line channel retains strong features displaced from common transitions at the catalog redshift. An unmatched line at an independent host redshift means no coincidence in the retrieved lists within the adopted velocity window. Wavelength coverage and completeness vary across those lists. A projected source can also contribute a line. The 5540 Å air red edge cuts off He II at about z = 2.378. Both C III] components fall within VIRUS only below about z = 1.903 in the requested interval. At higher redshift the metal test needs a redder spectrum. The evidence distinguishes an out-of-band component from a measured upper limit. A candidate galaxy spectrum cannot on its own determine a Population III or direct-collapse origin.

| Unassigned-line follow-up | Audited hypotheses | With external redshift | Unmatched at that redshift in retrieved lists | Of unmatched, both external atlases returned transitions |
|---|---:|---:|---:|---:|
| SC1 | 32 | 0 | 0 | 0 |
| SC2 | 32 | 0 | 0 | 0 |

[All external follow-ups](evidence/followup.csv) · [Rows passing every gate](evidence/candidates.csv). Each release’s compressed vetting ledger contains IDs, comparisons, rejected alternatives, and pending checks. Fluxes are in erg s⁻¹ cm⁻², spectral densities per Å, wavelengths in observed air Å and EWs in rest-frame Å. Each trial redshift depends on the named identification.

## Controls and recovery

2,594 of 3,733 successfully injected, unblended primary lines were recovered (69.5%) over the specified mix of amplitudes and widths. This is a spectral recovery fraction conditional on an extracted source. It does not include the probability of entering the public catalog, sky coverage, spatial completeness, or success in external follow-up.

| Stream | Spectra | Fitted peaks | He II leads | Extreme EW leads | Off-list leads |
|---|---:|---:|---:|---:|---:|
| SC1 observed | 232,650 | 141,783 | 1,689 | 34 | 1,879 |
| SC1 shuffled | 28,928 | 49,442 | 299 | 5 | 2,751 |
| SC1 off-source | 588 | 0 | 0 | 0 | 0 |
| SC2 observed | 1,107,763 | 399,422 | 5,098 | 31 | 5,791 |
| SC2 shuffled | 138,496 | 3,054,031 | 62,563 | 1,144 | 274,334 |
| SC2 off-source | 482 | 0 | 0 | 0 | 0 |

| Release | Matched real spectra / mock spectra | Real / mock leads per 1,000 spectra | Matched mock passes / leads | Raw-shuffle sample spectra | Raw mock passes / leads |
|---|---:|---:|---:|---:|---:|
| SC1 | 3,459 / 3,459 | 14.5 / 237.6 | 0 / 822 | 1,024 | 0 / 103 |
| SC2 | 16,503 / 16,503 | 9.5 / 130.8 | 0 / 2,158 | 1,024 | 0 / 2,611 |

Each lead is a channel hypothesis. A peak can qualify for multiple channels. Passes count hypotheses that completed confirmation. External follow-up of control leads used the same ranked first-pass rules. The control rate is conditional on its construction and the matched population.

| Published He II positive controls | Matched of eight | Single-line screen recovered | Weaker paired-line screen recovered | Recovered by both | Recovered by either |
|---|---:|---:|---:|---:|---:|
| SC1 | 2 / 8 | 0 | 1 | 0 | 1 |
| SC2 | 7 / 8 | 0 | 3 | 0 | 3 |

Recovery counts use unique published objects within each release. Each screen used its own S/N and pair cuts. SC1 produced 11 weaker spectral pairs against 82.0 expected from displaced-wavelength controls before vetting. It also tested 7,039 wrong-position, distinct-shot stacks and 16 usable blank position/wavelength pairs for halo injection. SC2 produced 38 weaker spectral pairs against 223.5 expected from displaced-wavelength controls before vetting. It also tested 18,791 wrong-position, distinct-shot stacks and 15 usable blank position/wavelength pairs for halo injection.

| Release | Nominal injected S/N | Clean trials | Primary line recovered |
|---|---:|---:|---:|
| SC1 | 5 | 129 | 9 (7.0%) |
| SC1 | 7 | 137 | 68 (49.6%) |
| SC1 | 10 | 125 | 118 (94.4%) |
| SC1 | 15 | 129 | 128 (99.2%) |
| SC1 | 25 | 128 | 128 (100.0%) |
| SC2 | 5 | 634 | 31 (4.9%) |
| SC2 | 7 | 599 | 303 (50.6%) |
| SC2 | 10 | 627 | 595 (94.9%) |
| SC2 | 15 | 622 | 611 (98.2%) |
| SC2 | 25 | 603 | 603 (100.0%) |

The raw wavelength-by-wavelength shuffle mixes bright and faint continua and destroys adjacent-pixel covariance. Its very large positive counts make it a stress test. After observing this excess, we added a same-shot block shuffle that preserves source continua. Its rates apply to the matched population. Off-source apertures preserve cube interpolation and spatial noise. Empty sky has no catalog redshift, so its off-list test is undefined. The per-release tables break recovery down by input flux, S/N, wavelength, profile width and input rest EW. Wilson intervals are nominal per-trial binomial intervals. Six trials share each injected background, and the intervals omit that clustering. The paired-visit controls substitute unrelated sources observed in the same shots.

## Inputs, evidence and reproduction

Every unresolved narrow He II / Lyα pair meeting the flux-ratio and spectral cuts, and every unresolved EW lead with a qualifying lower bound, entered the supplementary positional search of independent spectral catalogs. The table distinguishes successful catalog queries from attempts. Full rechecks preserve the original ledger and each revised decision.

| Release | Eligible source / line hypotheses | Distinct positions | Rechecked hypotheses |
|---|---:|---:|---:|
| SC1 | 21 | 21 | 4 |
| SC2 | 35 | 35 | 3 |

| Supplementary spectroscopy | SPARCL batches with valid response / attempted | CDS DESI batches with valid response / attempted | Distinct retrieved spectrum files by checksum | Eligible hypotheses with an external line fit | With an agreeing warning-free redshift | With an earlier agreeing DESI visual redshift |
|---|---:|---:|---:|---:|---:|---:|
| SC1 | 0 / 1 | 0 / 1 | 4 | 4 | 2 | 0 |
| SC2 | 0 / 2 | 0 / 2 | 3 | 3 | 1 | 1 |

For supplementary redshifts, agreement requires |Δz| < 0.006. A single spectrum or position can contribute to several eligible hypotheses. Independent line fits and catalog redshift agreement are recorded separately.

| Native DESI DR1 files and single-cone catalog retry | Positions rechecked | Matched | No DESI DR1 counterpart | Partial or unavailable | Source hypotheses with new evidence re-evaluated |
|---|---:|---:|---:|---:|---:|
| SC1 | 21 | 4 | 17 | 0 | 4 |
| SC2 | 35 | 3 | 32 | 0 | 3 |

The retry includes failed spectroscopic catalog queries and unresolved spectrum retrievals. The receipts distinguish files searched, sky pixels with no published product, empty target cones and failed requests. A completed DESI DR1 lookup constrains that survey only.

[SC1 static catalog audit](evidence/sc1/focus/static-summary.json) · [SC2 static catalog audit](evidence/sc2/focus/static-summary.json).

SC1 had 162 two-knot systems whose geometry and wavelengths recurred in distinct shots. Of these, 31 met the local and prior-catalog selection for independent cube checks, 22 had measurable bridges in two shots, and 2 passed both knot-centered spatial tests twice.
SC2 had 609 two-knot systems whose geometry and wavelengths recurred in distinct shots. Of these, 24 met the local and prior-catalog selection for independent cube checks, 19 had measurable bridges in two shots, and 0 passed both knot-centered spatial tests twice.

The endpoint audits give the positions, trial identifications and dispositions of every system that passed both spatial tests twice. [SC1 endpoint measurements and source checks](evidence/sc1/blob-followup/systems.csv) · [SC2 endpoint measurements and source checks](evidence/sc2/blob-followup/systems.csv).

| Release / system | First knot RA, Dec (deg) | Second knot RA, Dec (deg) | Observed line (Å) | Trial Lyα z | Disposition | Evidence or remaining ambiguity |
|---|---|---|---:|---:|---|---|
| SC1 / 44 | 197.779510, 51.198437 | 197.779739, 51.199913 | 4765.4 | 2.92107 | unresolved | Legacy query unavailable at secondary knot; direct imaging checks unavailable; foreground interpretation unresolved; independent redshift unavailable |
| SC1 / 78 | 11.112092, -0.050852 | 11.113163, -0.049214 | 5422.5 | 3.46173 | unresolved | Legacy query unavailable at secondary knot; direct imaging checks unavailable; foreground interpretation unresolved; independent redshift unavailable; persistent detector/fibre response not excluded at both knots; primary knot has a competing counterpart redshift; association unresolved; secondary counterpart redshift disagrees; association needs review; secondary knot has a competing counterpart redshift; association unresolved |

| Check on ranked features | SC1 | SC2 |
|---|---:|---:|
| NED cone searches with a valid response | 148 | 139 |
| SIMBAD cone searches with a valid response | 148 | 139 |
| Legacy Surveys cone searches with a valid response | 0 | 0 |
| Direct HSC cutouts with measured photometry | 30 | 29 |
| Cube opened and target narrowband window tested | 145 | 139 |

SC1 could not retrieve the HETDEX-specific DESI follow-up catalog and the published HETDEX Lyα-nebula catalog during the initial source and control vetting. That limited novelty checks and the imaging shortlist for this release. The independent spectral and repeated-geometry results above report their own coverage.

For each tested cube, the extension statistic subtracts the median from disjoint offset-wavelength windows and divides by 1.4826 times their median absolute deviation (MAD). The numerator counts source-row hypotheses with an excess ≥5 times this scaled MAD. The denominator counts source rows with a measurable extension statistic in their primary cube. Multiple source rows can refer to the same physical system.

| Spatial test on source rows | SC1 excess / measurable source rows | SC2 excess / measurable source rows |
|---|---:|---:|
| Emission between two knots | 5 / 16 | 0 / 13 |
| Diffuse emission around one center | 0 / 14 | 0 / 8 |

The release summaries retain each unresolved reason, paired-metal upper limits, and exact ledger rows for the strongest remaining leads.

[SC1 evidence breakdown](evidence/sc1/evidence-summary.json) · [SC2 evidence breakdown](evidence/sc2/evidence-summary.json).

| Release input | Bytes | SHA-256 |
|---|---:|---|
| hetdex_sc1_spec_v3.2.fits | 4,864,757,760 | `52b1d5f437305822658dfca1d441ea1f18a9e64e847b53228983d4431a7b786a` |
| hetdex_sc1_detinfo_v3.2.fits | 141,235,200 | `1cac07886617d5b5a8481df9571f1b493ad3b5821ca119a3a264e6065435ce28` |
| hetdex_sc2_spec_v1.5.fits | 9,406,042,560 | `d1ed1ccd4745136776f904659f2e88dfe8ca64285978a18fbad6b33004fa650c` |
| hetdex_sc2_detinfo_v1.5.fits | 1,839,265,920 | `468696f44e86dc1a40fc96ef363c0381a7f0394f116f09b830745fd50edbbad3` |

SC1 · [Input checksums and mirrors](evidence/sc1/inputs-sc1.json) · [Ancillary and cube receipts](evidence/sc1/vet-inputs.json) · [Vetting ledger](evidence/sc1/vetted.jsonl.gz) · [Injections](evidence/sc1/injections.jsonl.gz) · [Completeness](evidence/sc1/completeness.csv) · [He II positive controls](evidence/sc1/positive-controls.json) · [Weaker paired lines](evidence/sc1/weak-heii-summary.json) · [Imaging-led EW measurements](evidence/sc1/ew-rescue-summary.json) · [Supplementary spectroscopy](evidence/sc1/focus/summary.json) · [Query receipts](evidence/sc1/focus/receipts.json.gz) · [Eligible spectroscopy ledger](evidence/sc1/focus/catalog.csv.gz) · [Repeated spatial systems](evidence/sc1/repeat-blobs/summary.json) · [Raw-shuffle follow-up](evidence/sc1/shuffle-vet/vet-summary.json) · [Continuum-preserving control](evidence/sc1/nuisance-summary.json) · [Matched control vetting](evidence/sc1/matched-shuffle-vet/vet-summary.json) · [Halo recovery](evidence/sc1/halo-injections-summary.json) · [Repeat controls](evidence/sc1/repeat-summary.json) · [Off-source spectra selection](evidence/sc1/offsource.jsonl.gz).

SC2 · [Input checksums and mirrors](evidence/sc2/inputs-sc2.json) · [Ancillary and cube receipts](evidence/sc2/vet-inputs.json) · [Vetting ledger](evidence/sc2/vetted.jsonl.gz) · [Injections](evidence/sc2/injections.jsonl.gz) · [Completeness](evidence/sc2/completeness.csv) · [He II positive controls](evidence/sc2/positive-controls.json) · [Weaker paired lines](evidence/sc2/weak-heii-summary.json) · [Imaging-led EW measurements](evidence/sc2/ew-rescue-summary.json) · [Supplementary spectroscopy](evidence/sc2/focus/summary.json) · [Query receipts](evidence/sc2/focus/receipts.json.gz) · [Eligible spectroscopy ledger](evidence/sc2/focus/catalog.csv.gz) · [Repeated spatial systems](evidence/sc2/repeat-blobs/summary.json) · [Raw-shuffle follow-up](evidence/sc2/shuffle-vet/vet-summary.json) · [Continuum-preserving control](evidence/sc2/nuisance-summary.json) · [Matched control vetting](evidence/sc2/matched-shuffle-vet/vet-summary.json) · [Halo recovery](evidence/sc2/halo-injections-summary.json) · [Repeat controls](evidence/sc2/repeat-summary.json) · [Off-source spectra selection](evidence/sc2/offsource.jsonl.gz).

The workflow pins Python dependencies, fetches both releases, checks the publisher MD5 values and hashes the received bytes. Run its **full** dispatch on the research branch, or use the same commands on a runner with network access and sufficient disk. Run from the repository root.

```sh
D=results/hetdex_blind_emitters
python -m pip install -r "$D/requirements.txt"
python "$D/check_pipeline.py"
for r in sc1 sc2; do
  python "$D/scan.py" "$r" scan --output "$D/evidence/$r" --cache "/tmp/hetdex-$r" --workers 3
  python "$D/scan.py" "$r" vet --output "$D/evidence/$r" --cache "/tmp/hetdex-$r"
  python "$D/focus.py" "$D/evidence/$r" "$D/evidence/$r/focus" --release "$r" --cache "/tmp/hetdex-$r" --vet
  python "$D/repeat_blobs.py" "$D/evidence/$r" "$D/evidence/$r/repeat-blobs" --release "$r" --cache "/tmp/hetdex-$r"
  python "$D/focus.py" "$D/evidence/$r" "$D/evidence/$r/focus" --static-retry --vet --release "$r" --cache "/tmp/hetdex-$r"
  python "$D/focus.py" "$D/evidence/$r" "$D/evidence/$r/focus" --replay-saved
  python "$D/blobs_followup.py" "$D/evidence/$r" --cache "/tmp/hetdex-$r"
done
python "$D/report.py"
python -m open_space.verify
python -m unittest discover -s tests
```

[Spectral method](SPECTRAL_METHODS.md) · [Spatial, imaging and literature checks](VETTING.md) · [HETDEX source catalog 1 and public data](https://doi.org/10.5281/zenodo.7448504) · [HETDEX source catalog 2 and public data](https://doi.org/10.5281/zenodo.19581262) · [HETDEX data products and associated catalogs](https://hetdex.org/data-results/) · [Prior HETDEX He II search](https://arxiv.org/abs/2609.16134).

These measurements use the public source catalogs of Mentuch Cooper et al. (2023, ApJ 943, 177) and Mentuch Cooper et al. (2026, ApJS 284, 67), obtained with VIRUS on the Hobby–Eberly Telescope. The release records list the HETDEX collaboration, telescope partners, funding and requested acknowledgements.
