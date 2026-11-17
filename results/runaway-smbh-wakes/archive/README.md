# HST and JWST deep-field search

No new host-connected wake passed this archival imaging screen. One HST line passed the independent F606W and F814W profile cuts when each band fitted its own geometry. It failed the same test on a fixed sky axis in both follow-up ACS bands. The known RBH-1 feature was found in both discovery filters, but the full selection missed it. These measurements do not support a population upper limit.

## Completed

- The blind HST screen used five abutting, non-overlapping 1024 × 1024 windows in the 0.06 arcsec/pixel XDF mosaic, centered on the image and its four cardinal neighbors. All 5,242,880 pixels had positive weight in each band, covering 5.24288 arcmin² geometrically. Boundaries truncate lines. The JWST/JADES DR2 NIRCam F150W window covered a further band on 0.26205 arcmin² **inside** the central HST window. The observations are archival. The target selection for RBH-1 and the additional western follow-up stamp contribute zero blind area.
- We searched angles and centers within these windows using 4, 6 and 8 arcsec templates of 0.16 arcsec FWHM, a score threshold of 6, minimum measured length 3 arcsec and aspect ratio above 20. The HST search returned 637 F606W and 629 F814W raw shapes. There were 104 cross-band pairs satisfying the detector's segment-compatibility test in both directions. One pair also passed continuity in both fitted profiles. The JADES search returned 38 raw shapes, one passing continuity. Three raw shapes matched HST F814W. None of those pairs passed continuity in both bands. We inspected the central HST matched stamps, every JADES raw stamp, and the western field's sole bilateral continuity match. See [field totals](tables/field_totals.csv), [all triggers](tables/archive_triggers.csv), [matched pairs](tables/archive_pairs.csv) and [JADES stamps](tables/jades_raw.jpg).
- The western HST pair (F606W index 59, F814W index 74) is centered at 53.18577725°, −27.78513879°. F606W measured 4.44 arcsec long and an aspect ratio of 21.1. Fixing the F606W geometry across bands gives 6/8 positive-against-both-flanks bins in F606W, 3/8 in F775W, 3/8 in F814W, 2/8 in F850LP, and 2/8 in a targeted JWST/JADES F150W cutout. The originally independent F814W trigger passed only on its own fitted centerline. The stamps show separated patches and no compelling continuous structure pointing from a distinct host. DR10 and DR11 provide a shallow comparison. See [follow-up images](tables/western_followup.jpg) and [fixed-axis measurements](tables/western_followup.json).
- The nearest UVUDF non-star centroid to the western F606W endpoints is source 6118, offset by 1.067 arcsec. The independent F814W fit has a different nearest source (6833 at 1.153 arcsec). Its BPZ redshift is 3.38 with odds 0.463. Source 6732 along the fitted line has BPZ 4.24 (interval 3.28–4.53). Source 10275, 1.92 arcsec past the opposite endpoint and 0.07 arcsec off-axis, has BPZ 2.06 (1.91–2.81). These values and their uncertainties do not establish a common redshift. The independent MUSE HUDF DR2 catalogue identifies UVUDF 6732 as MUSE 4025, with z=4.04923, confidence 2, 0.11 arcsec off the axis and 2.59 arcsec from endpoint zero along it. Source 8020 has tentative z=2.51612 (confidence 1), 1.98 arcsec beyond the opposite endpoint and 0.04 arcsec off-axis. The catalogue supplies one redshift within the narrow on-segment strip. Individual source redshifts cannot be assigned to diffuse emission. See [25 MUSE cone matches, input hashes and coordinates](muse_matches.json) and [MUSE HUDF DR2, Bacon et al.](https://doi.org/10.1051/0004-6361/202244187). VLASS epoch-1 components and the RACS-low DR1 source catalogue each returned zero rows within 15 arcsec. This provides a shallow radio counterpart check. See [VLASS receipt](catalog_responses.json) and [RACS receipt](racs_cone.json). The catalogue table records the nearest non-star even when the association is distant or incidental. The complete UVUDF catalogue contains 9,969 sources, 42 flagged stars. See [UVUDF data and attribution, Teplitz et al. and Rafelski et al.](https://archive.stsci.edu/hlsp/uvudf).
- Eight fixed random sites/orientations per length in the central HST window were reused across five brightnesses (120 injections). Width was fixed at 0.16 arcsec. The XDF ACS F606W adopted AB zero point is 26.51. The numerator counts new matches passing the single-band continuity test. No endpoint galaxy was injected. These small conditional samples omit cross-band and host completeness. See [all injection trials](tables/xdf_inj.json) and [XDF calibration](https://archive.stsci.edu/prepds/xdf/).

| Length | 24 AB mag/arcsec² | 25 | 26 | 27 | 28 |
|---|---:|---:|---:|---:|---:|
| 4 arcsec | 7/8 | 6/8 | 4/8 | 2/8 | 0/8 |
| 6 arcsec | 7/8 | 7/8 | 4/8 | 3/8 | 2/8 |
| 8 arcsec | 8/8 | 7/8 | 6/8 | 6/8 | 0/8 |

## RBH-1 positive control

HST GO-16912 ACS/WFC products jety02010 and jety02020 are distinct F606W and F814W stacks. We chose a 19.2-arcsec field containing the published host and searched all centers and orientations within it, using 4, 6 and 8 arcsec filter lengths and widths of 0.12 and 0.20 arcsec. Comparison used the published host at 02h41m45.43s, −08°20′55.4″ and the southeastern direction along the published slit (147° east of north). The detector recovered the feature in both bands. Measured lengths reached the profile search limit. The adopted 1 arcsec endpoint-centroid rule missed the true alignment in both bands. See [measurements](tables/rbh1_measured.json) and [the original publication](https://arxiv.org/abs/2302.04888).

| Band | Rank among raw triggers | Length | Aspect ratio | Axis difference | Centerline offset | Closest endpoint to host | Continuity |
|---|---:|---:|---:|---:|---:|---:|---|
| F606W | 3 / 9 | 5.60″ | 26.1 | 0.45° | 0.006″ | 2.53″ | Passed, 7/8 bins, maximum flux fraction 0.217 |
| F814W | 4 / 8 | 9.20″ | 33.2 | 0.34° | 0.074″ | 1.14″ | Failed, 8/8 bins, maximum flux fraction 0.516 |

Forced grz measurements in Legacy DR10 and DR11 used the central 2–7 arcsec of the published axis, parallel flank subtraction, and robust scatter from 32 parallel displaced apertures. The largest positive contrast was 1.34 times the local scatter. This illustrates the ground survey's limited sensitivity for this object. The results test recovery of morphology only. Van Dokkum et al. interpret the JWST/NIRSpec data as a bow shock. Sánchez Almeida et al. argue that the line ratios also fit a star-forming disk. See [JWST bow-shock paper](https://arxiv.org/abs/2512.04166), [disk interpretation](https://arxiv.org/abs/2602.10682) and [archival data DOI](https://doi.org/10.17909/8d8q-r813).

## Reproduction

Use Python with NumPy, SciPy, Astropy and Matplotlib. Set DETECTOR to the directory containing wake_search.py and profile_screen.py. Run the commands from this directory. The MAST receipts include proposal response hashes, recorded query parameters and relevant rows. Full-image hashes and the native JADES science, uncertainty and weight slice receipts are in archive_inputs.json. The FITS files are regenerated from their public sources. The JADES source is approximately 11 GB. The program reads fixed byte ranges under the recorded object ETag condition and records both returned-range and cropped-image hashes. Retrieval at a future date may return updated data. Compare hashes before comparing outputs.

```sh
export DETECTOR=../detector
export PYTHONPATH="$DETECTOR"
python archive_inputs.py pixels --receipt new_inputs.json
python archive_screen.py "$DETECTOR" pixels raw --inject
python archive_audit.py pixels raw tables
python control_measure.py pixels raw
python mast_queries.py new_mast_receipts.json
python follow_inputs.py pixels new_west_inputs.json
python west_measure.py pixels raw tables
python catalog_follow.py new_query_receipts.json
python muse_match.py pixels raw new_muse_matches.json
```

The comparison plots were generated from calibrated pixels. Yellow circles show the ends of the fitted line. The input pointings, catalog releases and source links are recorded in [image receipts](archive_inputs.json), [catalogue receipt](host_catalog_receipt.json) and [MAST query receipts](mast_query_receipts.json).

## Pending

The search covers part of GOODS-S and a targeted positive-control field. No archival population completeness, remaining deep-field survey area, or upper limit has been measured. DESI DR1 has no covering tile for the 10.8 arcsec cone around the western feature (NESTED NSIDE 64 pixel 35935). [DESI coverage receipt](../desi/direct_cones.csv). Single-exposure persistence, quantitative disk/stream/arc model comparison, and spatially resolved spectroscopy were not established for the western feature. A true host-aligned wake can fail the endpoint criterion or the brightness-dominance veto demonstrated above.
