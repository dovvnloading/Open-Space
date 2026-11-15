# Fourth-epoch VLASS survey of ultracool dwarfs

The strongest new association from this survey is [WISE J0458+6434 AB](../wise-j0458-radio/README.md). A weak positive in the fourth VLASS visit prompted a check of earlier visits at the system's predicted positions. Emission follows the binary across the sky. The earlier-visit measurements, translated-position controls, source images and freely fitted stationary-source comparison are in the counterpart report.

## Sample and selection

We used the 3,890 rows in the Main table of UltracoolSheet 2.1.0. The published spectral-type index had to be at least 6 (M6 or later), with an optical or infrared type and without a photSpT flag. Of 3,500 such rows, 3,487 had both coordinates at 2025 July 3 and both proper-motion components. Intersecting those positions with CADC VLASS epoch 4.1 Quick Look Stokes I footprints yielded 1,526 catalog rows. The inventory contained 16,410 tile rows. A field needed a finite center, a polygon and a tt0 intensity FITS product. Separate catalog entries may refer to members of one radio beam. [UltracoolSheet](https://zenodo.org/records/15802304) · [VLASS](https://science.nrao.edu/vlass) · [CADC archive](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/vlass/)

Each position was propagated with the catalog proper motion to the image DATE-OBS. The fourth-epoch screening value is the brightest pixel center within 1.5″ of that position, after subtracting the median background. The noise scale is 1.4826 times the median absolute deviation in the 8–35″ annulus. Eight offset positions at 20″ use the same aperture. Off-image controls are excluded. At overlapping tile boundaries the search takes covering tiles in center-distance order and tries the next one if the initial cutout lacks enough science pixels. Annual parallax and binary orbital motion are omitted from this initial small-aperture screen.

An early distance-ordered batch contained the known L 726-8 radio system and two weaker features on very late-type targets. The two weak signals were chosen for earlier-epoch examination after their fourth-epoch values were seen. WISE J0458+6434 AB was the only one of these two with corroborating motion. The T8 dwarf 2MASS J07290002−3954043 had a 3.92 noise-scale peak in the fourth epoch. The 2020 and 2023 images yielded plane-and-beam fitted amplitudes of −0.21 and +0.13 mJy/beam (local scatter 0.28 and 0.11). These exploratory follow-ups were conducted while the rest of the survey was running. The field ranks in the WISE report apply to that source's field. They do not calibrate the selection of all possible targets or alternative tests. The positive-control system appears as G 272-61A in the Sydney Radio Star Catalogue. [Catalog](https://radiostars.org/data/)

## Coverage and sensitivity

All 1,526 covered catalog positions yielded a measurement. The median local noise was 0.140 mJy/beam. The six-noise-scale screen selected one catalog position and 0 of 12,202 offset measurements. The 3.681 peak-to-noise value at WISE J0458+6434 AB was equaled or exceeded by 5 catalog targets and 12 offset measurements. Its single-visit value motivated the examination of the earlier observations. The broad search measured 67 T8-or-later catalog positions.

Adding a synthetic beam-shaped source at ten times local noise to the *actual* image pixels of the T/Y sample recovered 2,700 of 2,704 injections at the six-noise-scale cut. The injections probe pixel sampling, beam shape and local backgrounds. They omit astrometric errors, intrinsic variability and calibration changes. [Injection replay](completeness.py) · [Injection counts](injections.json)

The parallax subsample in the JSON gives optional sensitivity benchmarks at spectral luminosities of 10¹⁴ and 10¹⁵ erg s⁻¹ Hz⁻¹. Eligibility requires a T or Y classification, a parallax signal-to-noise ratio of at least 10 and an expected on-axis flux density at least ten times the measured local noise. The simple upper bounds count threshold crossings of six noise scales in one visit, merge catalog positions within 4″, and assume independent detections of equally probable continuously bright systems. These are conditional binomial bounds. They do not model variability, optical selection, radio burst duty cycle or the sky area outside these images.

## Reproduction

Install Python with NumPy, pandas, SciPy and Astropy. These programs fetch public files from Zenodo and CADC and save their caches in the directory supplied. Replaying a live inventory can change ordering or products as the archive updates. The committed ledger fixes the object names, data product URIs, measured positions, offsets and observed values for this run. CADC sometimes regenerates FITS header comments on cutout delivery. The byte-level cutout hash identifies the payload actually measured. Fetch its original parent image to compare the science pixels. The original-pixel extraction adjusts the WCS reference pixel and preserves the beam metadata.

```sh
python prepare.py --cache /tmp/vlass4-ucd
python survey.py --cache /tmp/vlass4-ucd --workers 6
python finalize.py --cache /tmp/vlass4-ucd --out /tmp/vlass4-replay
python completeness.py --cache /tmp/vlass4-ucd
```

- [Machine-readable cutout ledger](cutout-ledger.csv.gz)
- [Coverage, control counts and sensitivity](coverage.json)
- [Input release and inventory digests](inputs.json)
- [Original branch checkpoint](search-plan.md)
