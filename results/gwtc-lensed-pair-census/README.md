# No validated strongly lensed gravitational-wave pair in a bounded GWTC-2.1/3/4.0 screen

**Claim: null for this pre-registered catalogue-summary test, not a full-posterior lensing exclusion.** Of **141** eligible catalogued black-hole mergers across version-pinned GWTC-2.1 (34), GWTC-3 (23) and GWTC-4.0 O4a (84), a mass/spin/time-delay proxy searched all **9,870** unordered pairs. The top pair, **GW190720_000836 / GW190728_064510** (8.275 days apart), was beaten or tied in **956/1,500** catalogue-stratified shuffled-background maxima: add-one family-wise false-alarm fraction **957/1501 = 0.637575**. Thus **no significant pair in this proxy test**. Full posterior-overlap significance has **not** been calibrated, so this cannot rule out lensing in these catalogues. No discovery is claimed.

**Separate, more expensive follow-up:** a frozen five-pair shortlist (three leading proxy pairs plus two published LVK controls) subsampled **1,400** deterministic PE draws per event from **10** original Zenodo HDF5 files; all downloaded files passed declared MD5 and length, with independent original-file SHA256 and extracted-sample SHA256 recorded. A **355,238,670-byte** O4a LVK lensing archive was independently MD5-verified (262 members). The job and read-only referee replay both completed; the scientific limitations below remain.

## Prior-pinned data and calculation

| Release | Immutable Zenodo record | Eligible events with PE |
| --- | --- | ---: |
| GWTC-2.1 (v2) | [6513631](https://zenodo.org/records/6513631) | 34 |
| GWTC-3 (v3) | [22685054](https://zenodo.org/records/22685054) | 23 |
| GWTC-4.0 corrected O4a (v3) | [17602505](https://zenodo.org/records/17602505) | 84 |
| LVK published O4a lensing controls | [18163632](https://zenodo.org/records/18163632) | external comparison |

The exact event list and GWOSC input-page hashes are in [`outputs/gwosc_event_freeze.json`](outputs/gwosc_event_freeze.json) and [`outputs/eligible_events.json`](outputs/eligible_events.json). [`outputs/release_manifest.json`](outputs/release_manifest.json) pins *every* eligible Zenodo PE filename, version, URL, MD5 and byte length; [`outputs/posterior_receipts.json`](outputs/posterior_receipts.json) pins every downloaded original-file SHA256, selected waveform sample group, extraction seed, column list and extracted NPZ SHA256. Re-running against changed source bytes requires an explicit new freeze.

For all pairs, GWOSC 90% summaries give a **preselection score** using detector-frame chirp mass, effective aligned spin and a simplified SIS galaxy time-delay prior (SDSS early-type velocity dispersion function, magnification weighting and observed-cadence window); its fixed shuffled null calibrates **only that statistic**. A separate Haris-style posterior integral `∫ p₁(θ)p₂(θ)/π(θ) dθ` is approximated with 6D KDE component masses and spin tilts/magnitudes, HEALPix sky and inclination. Those posteriors have **nonmatched effective PE priors**; the reported overlap B is an **exploratory scale**, not reliable lensing odds. The sum of approximate log10 overlap and the lens delay log-likelihood ratio is also not a calibrated detection statistic.

## Five posterior-sample shortlist checks (not a global search)

| Event pair | Delay, days | Approx. log10 overlap B | + delay score (log10) | Independent status |
| --- | ---: | ---: | ---: | --- |
| GW231113_200417 / GW231223_202619 | 40.015 | +5.94 | +6.40 | Published LVK **joint** log10 B **−3.3**; no lens |
| GW190720_000836 / GW190728_064510 | 8.275 | +4.49 | +5.24 | Proxy global p=0.638; sky overlap is smoothing-sensitive |
| GW191103_012549 / GW191204_171526 | 31.659 | +3.42 | +3.93 | No phase/strain/lens-galaxy validation |
| GW190421_213856 / GW190910_112807 | 141.576 | +3.33 | +2.95 | Published O3 joint lens evidence negative |
| GW191129_134029 / GW191216_213338 | 17.329 | −57.40 | −56.74 | Six-dimensional intrinsic posterior mismatch |

**The high approximate overlap for the O4a control is a warning, not evidence.** LVK's stricter joint analysis rejected it (Hanabi log10 lens Bayes factor −3.3). Our *different* O4a summary pre-screen placed it **second** among its 3,486 O4a-only combinations; that is **not exact reproduction** of LVK's ranking or joint evidence. Two older published O3 control pairs require GW190424_180648, outside the selected file release; these controls remain untested. Consult [`outputs/lvk_reference.json`](outputs/lvk_reference.json) and [`outputs/lvk_archive_inventory.json`](outputs/lvk_archive_inventory.json) for the exact published archive receipt. The archived LVK joint values are external references, not recomputed by this pipeline.

All five tested pairwise detector-frame mass 90% marginal intervals overlap; this alone is weak evidence. **Source-frame masses are not lensing invariants**: lens magnification biases luminosity distance/redshift inference. Independent coarse equal-area sky rebinning shows the leading proxy pair lacks robust high-probability patch agreement; the inference depends strongly on sky binning and smoothing, so a sky BF is not independently validated.

## Completed and pending gates

| Gate | Disposition |
| --- | --- |
| Version-pinned Zenodo list, metadata hashes, 141-event selection | **Completed** |
| Full 9,870-pair mass/spin/delay proxy + 1,500 real-catalogue shuffles | **Completed; independent exact replay**, global proxy p=0.638 |
| Ten original PE files verified by Zenodo MD5/size, SHA256 receipts, five posterior pair scores | **Completed on a frozen shortlist**, not all-pair PE |
| LVK O4a release archive MD5 and published top negative control | **Archive verified; qualitative control reproduced, joint odds not** |
| Detector- versus source-frame mass conventions | **Checked at marginal level**; common prior and waveform systematics pending |
| Morse/image-type phase and coherent GWOSC strain joint likelihood | **Not performed** |
| Foreground EM lens candidates, spectroscopy, image geometry | **Not performed** |
| Full-catalogue posterior Bayes factors, injection/reweighting and trial-corrected FAP | **Not performed** |

## Run, referee and limitations

Run the single Actions job by pushing a change to `run.py`, `test_census.py` or its workflow on `research/gwtc-lensed-pair-census`. Outputs-only commits are not triggers. The job is `.github/workflows/gwtc-lensed-pair-census.yml`; source is [`run.py`](run.py), synthetic tests [`test_census.py`](test_census.py). The independent [`referee_replay.py`](referee_replay.py), executed **outside the computing job** on the deposited frozen outputs, reproduces every all-pair proxy score and background maximum with no imports from `run.py`; its [receipt](referee_checks.json) and [review](REFEREE.md) distinguish performed from pending work. [Runner summary](outputs/RESULT.md), [machine-readable result](outputs/result.json), [full proxy rows](outputs/all_pairs_proxy.tsv.gz), [null model](outputs/proxy_background.json), [lens-delay model](outputs/time_delay_population.json) and [five-pair posterior table](outputs/posterior_pairs.json) are deposited.

**Limits:** simplified SIS lens prior/observing window, GWOSC summary preselection, approximate PE prior and smoothing, no general detector selection efficiency, no global posteriors-derived false-alarm statistic, and no phase/strain/foreground-lens joint validation. This is a **calibrated null only for the proxy screen** and cannot set a lensing population rate. It does not alter `results/gwtc-successive-mergers`. The repository's [remarkable findings policy](../../docs/REMARKABLE_FINDINGS.md) bars a discovery label or log entry at this evidence level.

References: [Haris et al. 2018](https://arxiv.org/abs/1807.07062) · [LVK O4a lensing study](https://arxiv.org/abs/2512.16347) · [LVK public O4a data products](https://zenodo.org/records/18163632).
