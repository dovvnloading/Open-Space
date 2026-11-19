# Referee of the control-first Euclid DR1-F reselection

**Pre-computation decision: HOLD the new population and candidate claims.** Source-data per-object audit passes: 24 matched quasars, 14 matched dwarfs; exact initial 3/24 recovery; overlapping losses 20 due to Y−J, 10 due to |VIS S/N|, plus flags/ellipticity/VIS_DET. The frozen control-tuned relaxed cut yields **21/24 apparent matched-quasar recall and 5/14 matched-dwarf retention**. A nonrepresentative, label-tuned control sample offers *no held-out completeness or field-purity calibration*.

**Release gates:** full 768-partition relaxed-source reread and hash reconciliation; statistically bounded top-list register; SIA-to-S3 FITS pixel receipts (not just metadata); CatWISE/unWISE and public spectrum/literature status with error-vs-nondetection distinctions. The [machine referee](referee.py) checks the source manifests and independently replays all matched controls. If follow-up services or real pixels fail, report that failure: do **not** claim image clearance, spectral exclusion, calibrated posterior or discovery.

**Scientific evidence bar:** only a *new, independently verified spectrum and redshift* could support a discovery. Photometric ranks remain follow-up targets. Acceptable archival result, even with incomplete pixels, is strictly **measured/inconclusive with explicit open gates**.
