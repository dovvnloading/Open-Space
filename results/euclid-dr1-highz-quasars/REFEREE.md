# Independent referee assessment — Euclid DR1-F high-z QSO pilot

**Verdict:** ACCEPT as an **auditable bounded source-catalogue screen and negative *discovery validation***; REJECT any assertion of a novel z>7 quasar, a calibrated quasar-versus-dwarf probability, or a statistical null for the full DR1 release. No `[REMARKABLE]` classification.

## Checks reproduced separately from the primary selector

The standalone [referee replay](referee.py) reads deposited `s3_output/broad_rows.csv` and applies the signed VIS/J/H error, Y−J, J−H, ellipticity and image-flag cuts without importing the main selector. On the networked single Actions job it reproduced **37,443 broad input rows → 290 strict proxy rows**, and the exact five highest-J-SNR ranked IDs, verified `freeze.json`, `run_s3.py`, the downloaded HATS metadata and output SHA-256 manifests, and reconciled all eight publisher Parquet partition row counts.

The reference quasar **EUCL J125308.55+705432.3, published z≈7.69**, was independently recovered in the strict sample at **0.051219 arcsec**. Controls at published z≈7.77 and z≈7.61 were not in this fixed DR1-F sample, hence they cannot certify completeness or failure. The nominal sum of cone areas, 2.014698 deg², reflects separately rounded cone values. The Fornax deep probe intersects no HATS partitions. The **1,900 deg² released footprint is not 1,900 deg² searched**.

## Major scientific objections

1. **Unphysical shortlist colours:** all five top-J-SNR rows have J−H roughly −2.3 to −2.7 mag; **41 of 290** strict rows have J−H<−0.5. The original frozen query set only an upper limit. The after-the-fact −0.5 threshold and resulting **249** rows are diagnostics, **not preregistered selection or a validated candidate set**.
2. **Pointlike not proven:** Euclid `POINT_LIKE_PROB` is NaN for NIR-only objects. Low `ELLIPTICITY` is roundness rather than NISP PSF compactness. It does not separate dwarfs from compact galaxies or cosmic rays.
3. **Probabilities lack validity:** the displayed toy `P(QSO | QSO or dwarf)` is based on simple top-hat quasar tracks, fabricated L/T/Y colour centres and an arbitrary prior. A numeric value near 1 is **not** evidence of a quasar. Physically calibrated posterior probabilities are unavailable.
4. **Archival follow-up unfinished:** owing to TAP 502/timeouts, public W1/W2, proper motions, independent optical detections, radio/X-ray and spectra archives were **not checked**. A missing catalogue flux is not a confirmed non-detection.
5. **Selection function incomplete:** eight HATS partitions contain sources for fixed cones, not a characterized imaged area; neither injection completeness, known z>6.5 quasar recovery in the entire footprint, nor VIS ghosts/persistence/cosmic-ray vetoes were demonstrated. Broad-row CSV and S3 object ETags permit reproducible checking of the compact sample, but **S3 multipart ETags are not full Parquet SHA-256**.

**Recommendation:** merge the precisely scoped *method + measurements + limitations* as a non-remarkable, explicitly unfinished DR1-F pilot. No new spectroscopic redshift or quasar candidate can be responsibly published from this evidence alone. Resume only with independent image vetting, realistic synthetic filters/dwarfs, WISE and spectra checks, a measured footprint, and full positive-control recovery.

External grounding: [IRSA DR1 product identifiers](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_dr1.html), [IRSA HATS archive](https://irsa.ipac.caltech.edu/euclid-archive-documentation/data-products/), [Yang et al. 2026](https://doi.org/10.1051/0004-6361/202658883), [remarkable bar](../../docs/REMARKABLE_FINDINGS.md).
