# Mars 2026 gravity-model audit: ~100 m reference-frame discrepancy explains 99% of large-scale differences

**Outcome: a new, quantitatively tested gravity-catalog alignment result, not a new physical gravity anomaly.**  
**Study date:** 2026-11-15. **Data release:** 2026-09-24, NASA PDS.  
**Scope:** Independent new result directory. No prior Open-Space investigation was modified.

## What we found

We acquired and parsed NASA/JPL's complete degree/order-120 gravity models **MRO120F (2020)** and **MRO120G (2026)**, using all 7,380 harmonic entries in each. Blindly subtracting their coefficients makes Mars's large-scale gravity appear to have shifted substantially. Fitting just harmonic degrees **2 through 10** yields an effective **6.06689-arcsecond longitude-frame rotation**, equivalent to **99.887 metres of arc distance at Mars's 3396 km reference radius**.

The result matters because the two archived models use **different definitions of Mars's prime meridian**. The MRO120F PDS label anchors the meridian to Viking 1; the MRO120G label says its meridian matches IAU 2011 to within five metres. The original 2026 investigators already published orientation and gravity-spectrum differences. **Our distinct addition is a frozen-fit, degree-held-out numerical measurement of the spurious gravity-model mismatch created by omitting a frame alignment.**

| Coefficient degrees | Rows | Naive squared-coefficient mismatch removed by frozen 6.067″ longitude correction | Equivalent radial-acceleration difference at 400 km, before → after |
| --- | ---: | ---: | ---: |
| 2–10 (fitting) | 63 | **98.55%** | 6.019e−8 → 7.055e−9 m/s² |
| 11–20 (not used to fit) | 165 | **99.08%** | 2.518e−8 → 2.295e−9 m/s² |
| 21–40 (not used to fit) | 630 | **88.45%** | 1.545e−8 → 3.453e−9 m/s² |
| 41–60 | 1,030 | 8.84% | 3.239e−9 → 2.854e−9 m/s² |
| 61–79 | 1,349 | 0.15% | 2.280e−9 → 2.276e−9 m/s² |
| 80–120 | 4,141 | ≈0% | 1.525e−9 → 1.525e−9 m/s² |

**Definition:** The percentage is 100 times [1 minus squared residual coefficient-norm divided by squared raw difference-norm], using fully normalized C/S coefficients. At harmonic order m, the fitted correction rotates complex C+iS by exp(i m Δlongitude). The acceleration numbers apply the degree-dependent *radial* gravitational acceleration factors to fully normalized coefficient power at 400 km above the archive's 3396 km reference radius. These are descriptive field-difference metrics, **not probabilities, independent-sample statistical significance, or the share of Mars's total gravity field explained**. The GM monopole is deliberately excluded from the harmonic table.

As a cross-check, independently fitting the **held-out coefficient degrees** estimates **6.03077″** at degrees 11–20 and **6.03769″** at degrees 21–40, close to the frozen 2–10 fit. We did **not** refit the rotation for the held-out percentage calculations. Excluding degree 2 from the fit gives **6.07173″**; giving each harmonic degree equal angular leverage gives **6.04427″**. These degree partitions are not independent spacecraft data; the underlying two model inversions share most tracking observations.

There is a critical exception: **the 80–120-degree mismatch is not explained by longitude rotation**. Its coefficient-band norm falls ~**1.50%** between releases. Both NASA labels specify a Kaula-type regularization starting at degree 80. The 2026 solution incorporates additional data, revised reference-frame assumptions and different force/ephemeris modeling. We **cannot attribute these remaining differences to physical changes in Martian crust or mantle** on the evidence presented.

We separately measure a **G−F system-GM revision of −0.000571598 km³/s²**, which produces a literal monopole radial-acceleration difference of **3.967e−8 m/s²** at 400 km. The PDS files call this Mars-system GM (Mars, Phobos and Deimos), and it cannot be eliminated by rotating the coefficients. This is a **parameter/model difference, not proof that Mars lost mass**.

## Why this is a useful new result

It quantitatively exposes a major confounder in direct comparisons of public Mars gravity products. If old and new coefficients are subtracted as though their longitude frames are identical, ~99% of the apparent 11–20 degree power is a false change. That is an especially consequential screening error when one is hunting for gravity anomalies or validating new gravity maps. The next analysis should spend its effort on the **aligned residual** and physical measurements, not on a misplaced prime meridian.

The new 2026 solution adds radio tracking through February 2025; MRO120F ended in August 2019. A full physical comparison would still require the *time-dependent three-dimensional* body-frame conversion, common GM, seasonal mass exchange/tide conventions, matched ephemerides and force models. Our constant longitude angle is only the dominant first-order alignment control. This investigation does **not** independently discover a Mars rotation-rate variation: Konopliv et al. (2026) have already reported that.

## Data coverage and exact provenance

**Eligibility:** all (degree,order) coefficients with degree 1–120 and order 0–degree, no exclusions or cherry-picked locations, **7,380 rows per model, 14,760 total**. Fit degrees 2–10 to omit degree-1 reference/translation terms; the remaining degrees are reported in disjoint bands including null and negative results. Both model labels are validated for 3396 km reference radius, full normalization and 7,380 coefficient rows.

| Input | NASA PDS table | Archived source SHA-256 |
| --- | --- | --- |
| JPL MRO120F | [2020 coefficient file](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120f_sha.tab) | ddd3de9c30d75879fe37aa17a1149e7c96c141c095962954cb7ea865a2c025b6 |
| JPL MRO120G | [2026 coefficient file](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120g_sha.tab) | b296318066d6e8b8ab5a3ca4878d99635f1cf0e977c39e3efb3f36bfcb06a569 |

[2020 PDS label](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120f_sha.lbl) · [2026 PDS label](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120g_sha.lbl) · [Release history (2026-09-24)](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/errata.txt).

[Machine-readable measurement and all limits](source-model-comparison.json) records exact PDS source URLs, table and label byte counts and SHA-256, 120-degree completeness, fitted alternative angles, all seven degree bands and GM separation. Only compact derived results and code are checked in; the original raw bytes remain at NASA PDS.

## Reproduce without model inference

Standard-library Python 3.10+; no package installation needed. The code verifies *actual source bytes* using the recorded SHA-256 digests and refuses silently altered files or caches.

    python -m unittest discover -s results/mars-gravity-mro120g-2026 -p test_analyze.py -v
    python results/mars-gravity-mro120g-2026/analyze.py --cache /tmp/open-space-mro120g-pds --output /tmp/mars-gravity-mro120g-summary.json
    python results/mars-gravity-mro120g-2026/gmm3_crosscheck.py --cache /tmp/open-space-mro120g-pds --output /tmp/mars-gmm3-control.json

The [GitHub Actions archive replay](../../.github/workflows/mars-mro120g-audit.yml) runs synthetic rotation, normalization and parser tests before downloading both original source tables and labels, then writes the full comparison and uploads a result artifact. These measurements were generated from the verified NASA files in GitHub Actions, **not simulated gravity data**. The synthetic data exist only to test the numerical algorithm.

## External-center control: NASA Goddard's independently inverted 2016 GMM-3

To stress-test whether the JPL F→G offset is merely a quirk of comparing two releases from one team, we independently acquired NASA/GSFC's [GMM-3 2016 SHADR table](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/gmm3_120_sha.tab) and [original label](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/gmm3_120_sha.lbl). That adds **7,380 separately modeled coefficient entries** to the experiment, or **22,140 across all three archived solutions**. GMM-3 is developed at Goddard rather than JPL, but **the models draw on overlapping MRO, Odyssey and MGS tracking data** and its label specifically warns that its reference frame is *NOT* the IAU system.

| Source → target | Fitted effective yaw, degrees 2–10 | Fraction of raw 11–20 squared mismatch removed using the frozen yaw |
| --- | ---: | ---: |
| Goddard GMM-3 (2016) → JPL MRO120F (2020) | **−6.37144″** | **96.34%** |
| Goddard GMM-3 (2016) → JPL MRO120G (2026) | **−0.30455″** | **3.62%** |
| JPL MRO120F (2020) → JPL MRO120G (2026) | **+6.06689″** | **99.08%** |

**Consequential cross-center observation:** On the same harmonics, the independently derived Goddard model is only about **0.305 arcseconds (~5.0 m equatorial arc distance)** out of effective yaw from JPL's *2026* solution, but **6.371 arcseconds (~104.9 m)** from JPL's *2020* field. The 2020 frame mismatch is therefore not confined to a comparison against 2026 from the same producer. The three pairwise best-fit angles close to numerical precision (~5.6e−8 arcseconds), but this is an algebraic consistency diagnostic, not independent scientific significance.

**Caveats:** Goddard's gravity coefficients are an independently *modeled* product with overlapping tracking observations. A single fit angle does not prove globally identical 3D frames, an absolute orientation, a correct prime meridian, or true mass change. In particular, the GMM-3 label explicitly says its 2006-based reference system is **not IAU**; the relative alignments above must not be mistaken for literal IAU coordinates. Different reference periods, regularization and estimated parameters influence these fit angles.

The replayable [Goddard/JPL three-model control](gmm3-control.json) contains *all* pairwise held-out statistics, Goddard model GM, exact 900,604-byte source receipt and verified SHA-256 values. [Reproduction code](gmm3_crosscheck.py) independently downloads original 2016 PDS files and checks their source bytes. Pinned GMM-3 table SHA-256: eb4913b1afb6682406e6a9dad5be7918a162fa8462473c9a2e7aae258d4c2c9c; label SHA-256: 097f9989fc87bc9adbf04be754734349160228de940f9a242164f3be1b80bb7d.

## Falsification, unresolved alternatives and the next major discovery test

**Supported interpretation:** this large-scale discrepancy is primarily a coordinate-frame effect; frame choices are explicitly documented in the archive labels. The 6.067″ angle is an *effective model-registration yaw*, not an independently measured physical rotation or absolute prime-meridian displacement at a dated epoch.

**Other explanations for the residual:** changes in observed tracking arcs, Mars's pole orientation, different constraints/regularization (especially ≥80), solar-system ephemeris, spacecraft force models, CO₂ seasonal mass exchange and the tide conventions. Two overlapping tracking-data inversions are **not independent observations**; treating their coefficient formal sigmas as independent for a p-value would be unsound.

**Next falsifiable test before claiming a new physical gravity phenomenon:** express MRO120F and G in the same complete, epoch-dependent 3D Mars-fixed frame, use identical GM and tide conventions, and propagate both solutions through **eligible radio-tracking arcs *not included in either fit***, including suitable post-February-2025 NASA PDS MRO/Odyssey data. Compare Doppler/range likelihood with identical tracking geometry and force modeling, using a predeclared independent holdout. If the effect vanishes under those controls, it is not a physical discovery. No new mass anomaly, Mars mass loss or navigation-error rate is claimed here.

**Published prior work (do not rebrand as new):** [Konopliv, Park, Le Maistre, Yseboodt & Baland (2026), *Updated Mars Orientation and Gravity Field and Indication of Decadal Changes in Rotation Rate*, DOI 10.3847/PSJ/ae9ba7](https://doi.org/10.3847/PSJ/ae9ba7). In particular, [Figure 4 explicitly compares the models' prime meridians](https://www.astroexplorer.org/details/psjae9ba7f4) and [Figure 8 compares their gravity power spectra](https://www.astroexplorer.org/details/psjae9ba7f8). This PR contributes the *degree-held-out numerical frame-confounder measurement*, not the first identification of those two models or their already published orientation revision.
