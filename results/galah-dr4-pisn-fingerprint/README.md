# GALAH DR4: a coverage-limited search for a pair-instability-supernova fingerprint

**Result (18–19 November 2026): ZERO of 196 metal-poor GALAH DR4 stars with all five required unflagged elemental measurements passes the preregistered Na/Mg–Al/Mg–Co/Fe–Zn/Fe fingerprint.** No new PISN-descendant candidate is validated. This is a **measured-element screening null, not a population-wide PISN exclusion**: only 3/196 are at [Fe/H] ≤ −2 and none is at [Fe/H] ≤ −2.5; only 14/5,189 selected stars have all seven elements including Sr and Ba. PISN-versus-core-collapse/Type-Ia yield likelihoods and APOGEE/LAMOST spectroscopic cross-checks **were not completed**. No [REMARKABLE] finding is claimed.

## Frozen, completed search and counts

- **Primary source:** original GALAH DR4 `galah_dr4_allstar_240705.fits` (full published 917,588-star FITS; one star per record), not the derived Hugging Face mirror. HTTP 200 with **757,998,720 bytes**, Last-Modified `2024-07-10T08:28:59Z`; exact SHA-256 `97911d4436a42c77eb6d8c2116eecfb28fa927be83302aa66bbeb58d95889693`. The response URL, size, ETag, digest, query (read every row) and UTC run time are in [the machine-readable receipt](run/scan_summary.json). The raw FITS is not checked into git; its digest is the fixed source identifier.
- **Predeclared cuts:** [`criteria.json`](criteria.json) committed before access. `−5.5 ≤ [Fe/H] ≤ −1.5`; `flag_sp == 0`; GALAH CCD3 signal-to-noise `>30`. Each required `[X/Fe]` requires its own `flag_X_fe == 0`, a finite reported value, and `0 < e_X_fe ≤ 0.35 dex`. `flag_fe_h` is intentionally *not* used because DR4 documents a bug. Every upper limit, non-detection and missing value stays outside the detection-only gate.
- **Fingerprint gate:** `[Na/Mg] ≤ −0.70`, `[Al/Mg] ≤ −0.50`, `[Co/Fe] ≤ −0.25`, `[Zn/Fe] ≤ −0.20` (five separately detected elements: Na, Al, Mg, Co, Zn); the seven-element extension also requires measured `[Ba/Fe] ≤ 0` and `[Sr/Fe] ≤ 0`. These are **screening thresholds**, not yields or formal diagnostic likelihoods.

| Stage | Number |
| --- | ---: |
| All GALAH DR4 allstar rows acquired | 917,588 |
| In the frozen [Fe/H] interval (before other QC) | 11,181 |
| Metal-poor, `flag_sp=0`, CCD3 SNR >30 and finite [Fe/H] error | 5,189 |
| Na, Al and Mg all detected with individual QC | 217 |
| Na/Mg and Al/Mg both pass odd-even cuts | 9 |
| All five core elements detected with individual QC | **196** |
| All five detected and four fingerprint conditions met | **0** |
| All seven (including Sr and Ba) detected | **14** |
| All seven detected and conditions met | **0** |

Within the **5,189** quality-selected low-metallicity rows, individual usable detections are **Na 2,469; Al 217; Mg 3,582; Co 1,460; Zn 3,594; Ba 4,848; Sr 17**. In particular, **5,165 Sr values are flagged undetectable**, not low. The exact per-element flag breakdown appears in [`run/scan_summary.json`](run/scan_summary.json).

## Independent checks and abundance-space sensitivity

[`run/complete_five_rows.csv`](run/complete_five_rows.csv) contains **every one of the 196 measured five-element vectors**, including flagged-in-list `seven_gate` and uncertainties, for independent checking; [`run/five_gate_survivors.csv`](run/five_gate_survivors.csv) contains the header and **zero** detections. The [independent referee script](referee.py) re-evaluates every deposited detection and its uncertainty, exact gate, duplicate ID, [Na/Mg] and [Al/Mg] arithmetic, metallicity coverage, receipt SHA, and the disputed-star control without calling the screening code; its [report](referee_report.json) is produced by CI.

In the complete-five subset, `min [Co/Fe] = −0.1546` — no star meets the `−0.25` Co requirement (even though 9/196 meet the joint Na/Mg and Al/Mg cuts). `min [Zn/Fe] = −0.3049`, with five below −0.20. **Relaxing all four five-element thresholds by +0.10 or +0.20 dex still leaves zero simultaneous matches**, but +0.30 dex produces nine; threshold dependence is real. The most Co-depleted vector, `sobject_id=160422003001228`, `Gaia DR3=5401818786664214144`, `[Fe/H]=−2.2824`, has `[Na/Mg]=−0.7595`, `[Al/Mg]=−0.5671`, `[Co/Fe]=−0.1546`, `[Zn/Fe]=+0.0686`; it **fails** two critical conditions and is **not** a PISN candidate.

For the frozen artificial abundance vector ([Na, Mg, Al, Co, Zn, Ba, Sr]/Fe = [−1.2, 0, −1.0, −0.6, −0.55, −0.6, −0.6] dex), 256 Monte Carlo perturbations per star recover the **five-element** gate in 50,122/50,176 trials (**99.892%**) *conditional on all five measurements being available*, and the **seven-element** gate in 3,569/3,584 trials (**99.581%**) conditional on seven measurements. Multiplying by observed flag-eligible fractions yields descriptive coverage×recovery proxies **3.77%** (five) and **0.269%** (seven) of the 5,189 base stars. **They are not astrophysical completeness estimates:** abundance-dependent detection and flagging biases, metallicity selection, NLTE/systematics and unmeasured stars have no calibrated injection response. For 0/196, the exact 95% one-sided binomial upper bound is ≈1.52% **only in an identically selected, independently drawn complete-five population**; it cannot constrain the parent 917,588-star sample or the metal-poorest halo.

## Disputed control: LAMOST J1010+2358

[`j1010_xing_table1.json`](j1010_xing_table1.json) reproduces Xing et al. (2023, [Nature, doi:10.1038/s41586-023-06028-1](https://doi.org/10.1038/s41586-023-06028-1)) Table 1: [Fe/H] = −2.42, [Na/Fe] **< −2.02**, [Mg/Fe] = −0.66 ± 0.12, [Co/Fe] = −0.72 ± 0.12, [Zn/Fe] **< −0.25**, [Sr/Fe] **< −2.25**, [Ba/Fe] **< −1.37**. Therefore the tabulated central-ratio **limit** is [Na/Mg] **< −1.36**, not a measured Na detection. Aluminum is absent and zinc is only an upper limit, so **this published pattern cannot satisfy our strict five-detection gate**. [`control.py`](control.py) checks the published numbers and writes [`j1010_control_receipt.json`](j1010_control_receipt.json). This is a **literature-figure arithmetic control**, not independent recovery from a new LAMOST spectrum. Skúladóttir et al. ([2024](https://arxiv.org/abs/2404.19086)) and Jeena & Banerjee ([2024](https://doi.org/10.33232/001c.124113)) contested a pure-PISN explanation using revised chemistry and alternatives.

## Completed versus pending, and interpretation

**Completed:** prospectively committed thresholds; official full-catalogue acquisition with HTTP/byte/SHA-256 receipts; element-level flag gating and denominators; deposited all 196 complete-five rows and 0 survivor rows; Monte Carlo *measurement-conditional* sensitivity; contested J1010 published-data arithmetic control; an independently written arithmetic referee and its written limitations. All of these are locally replayable via the scripts and files here.

**Not completed:** quantitative fits to actual Heger & Woosley (2002) and Takahashi et al. (2018) PISN yield grids versus CCSN and Type Ia alternatives with common abundance covariance, censored limits, dilution, posterior/likelihood calibration; APOGEE DR17 and LAMOST individual-star remeasurement/zero-points; spectral line residuals, binarity/RV scatter, blending, SIMBAD/SAGA survivor vetting, or an independent full-FITS rescan. **There are no five-element survivors to vet**, but none of these tests can be counted as passed; in particular the published J1010 table does not substitute for a modern spectral control. No model likelihood ratio is reported. The proposed datasets, preconditions and likelihood specification are in [model comparison protocol](model_comparison_protocol.md).

**Meaning:** the frozen high-specificity fingerprint is absent **among the small subset with the required unflagged detections**. This is a reproducible abundance-screen measurement with substantial censoring and low-metallicity incompleteness, **not evidence that PISN descendants do not exist**. No [REMARKABLE] log entry was made.

## Reproduction

```sh
python -m pip install 'numpy>=1.26,<3' 'astropy>=6,<9'
python results/galah-dr4-pisn-fingerprint/screen.py \
  --cache /tmp/galah_dr4_allstar_240705.fits \
  --out results/galah-dr4-pisn-fingerprint/run
python results/galah-dr4-pisn-fingerprint/control.py
python results/galah-dr4-pisn-fingerprint/referee.py
```

This reacquires the **complete** 758 MB official input unless the local cache already exists; verify the fixed digest and the full 917,588 rows before interpreting any numerical changes. GitHub Actions [`galah-dr4-pisn.yml`](../../.github/workflows/galah-dr4-pisn.yml) independently executes the same frozen command and deposits receipts. A source mismatch or failed query is an error, **not** a null result.
