# Public JWST archive inventory

At the release cutoff **2026-11-17 06:51:31 UTC**, MAST returned 5,299 coronagraphic filesets over all targets, 4,192 normal-imaging filesets labelled `Star` with at least 300 s effective exposure, and 1,008 normal-imaging filesets from the supplementary program list. The queries overlap. After removing combined associations and deduplicating by MAST archive ID, the inventory contains **9,056 distinct exposure filesets, 2,653 visit/optical-element groups, and 240 programs**. A fileset can bundle both NIRCam channels. Each exposure contributes its integration time once.

The full [visit inventory](archive/visits.csv) preserves program ID and title, target and coordinates, observation and visit numbers, earliest and latest exposure and release, mode, optical elements, and an example dataset. Rows are sorted by most recent public release. The three compressed CSVs retain all matching public archive rows, including combined associations, individual fileset identifiers and rights. `release_last` is the date on the selected records at retrieval, which can differ across products in an observation. Publication dates and proposal numbers do not determine release order.

## Coverage and target selection

* The coronagraph query spans `MIR_4QPM`, `MIR_LYOT`, and `NRC_CORON`, without a target-name or depth cut. It includes 1,200, 57 and 3,831 distinct primary exposures respectively. Science, calibration, PSF, and background pointings are retained in the full table.
* The ordinary-imaging query covers `MIR_IMAGE` and `NRC_IMAGE` with the archive's proposer target class `Star` and a per-fileset effective exposure of at least 300 seconds. A supplementary query lifts both cuts in programs 1184, 1911, 3964, 4403, 4857, 6122, 7564, 7833, 8142, 8581 and 9079. This captures the short individual dithers in the MEAD survey. These cuts miss some short stacked observations, misclassified targets, off-axis stars, and stars outside the defined programs. Eclipse and parallel imaging may match the archive cuts, so each image still needs field, subarray and coronagraph checks.
* [Nearby fields](archive/nearby-fields.csv) lists **792 visit/filter groups from 70 programs**, with **331 distinct nearest Gaia source IDs**. A geometric host match requires a position within 5 arcsec of a Gaia Catalogue of Nearby Stars source with median distance below 100 pc. Eight groups have more than one candidate host in the matching radius. The table omits COM/CAL programs and records explicitly flagged as PSF references or background. Blank role flags do not establish a science target. For very bright stars absent from Gaia, tight binaries, catalogue omissions or offsets, use the full visit inventory. Moving-group membership is unverified unless established separately.
* The Gaia comparison uses sky positions propagated linearly from J2016.0 to the observation epoch, and also tests J2016.0 and J2000.0 to accommodate heterogeneous pointing coordinates. The table reports the smallest positional residual, tested epoch and number of distinct matches. Perspective acceleration and observer parallax are omitted at this field-selection stage. A pointing match identifies a potential host only. [Gaia Collaboration 2021](https://doi.org/10.1051/0004-6361/202039498)

## Prior searches

[Program–paper crosswalk](archive/program-papers.csv) contains 438 links to 334 distinct publications returned by 240 individual JWST program queries to MAST Papertrack. It returned records for 106 programs. This establishes a program association. A paper may use another target in a large survey, or use reference frames only. An empty receipt leaves target-specific prior-search status open. [Papertrack](https://mast.stsci.edu/papertrack)

For the four latest public coronagraphic targets, [SIMBAD's target bibliography](archive/go8826-host-bibliography.csv) returned 803 object–paper associations. The four targets in GO 8826 released on 2026-10-29 are all nearby (5.208–9.102 pc) and have older companion-search literature. The program-specific Papertrack query has zero records. SIMBAD links G 99-49, L 788-34 and Ross 619 to Oppenheimer et al.'s 8 pc coronagraphic survey, and V374 Peg to Durkan et al.'s Spitzer search. These associations identify prior coverage to inspect before assigning novelty. [Oppenheimer et al.](https://doi.org/10.1086/319941) · [Durkan et al.](https://doi.org/10.3847/0004-637X/824/1/58)

The [NASA Exoplanet Archive host search](archive/go8826-confirmed-planets.csv) identifies two previously reported RV planets around L 788-34, also known as GJ 4274. Their published semimajor axes are 0.0153 and 0.186 au (0.002 and 0.026 arcsec at the Gaia distance). The other three hosts returned no rows in the confirmed-planet table. This table does not cover all stellar or substellar companions. SIMBAD name/coordinate results, the explicit Gaia-ID lookup (important for fast-moving Ross 619) and exact queries are included. [Charpentier et al.](https://doi.org/10.1051/0004-6361/202557666)

## Reproduce

Run these commands from the directory containing the scripts (Python 3 with NumPy, SciPy and pandas 2.x). The GCNS TAP download is ~27 MB. Each query receipt records exact predicates, request time and response hashes. MAST's live outputs may change after reprocessing. The compressed CSVs serve as the fixed inputs to the visit and host tables. Table generation independently reconciles page counts and removes associations before summing exposure times.

```sh
python fetch_inventory.py coron nearby_programs stellar_imaging
python nearby_inventory.py --gcns /tmp/open-space-gcns.csv coron nearby_programs stellar_imaging
python program_papers.py
python latest_host_search.py
```

The main archive search, host-position crossmatch, per-program bibliography and four-target catalogue queries are complete. For hosts outside the selected reduction, image inspection, second-epoch verification, and target-specific literature adjudication remain pending. No entry in an inventory table constitutes a new companion or an occurrence-rate constraint.
