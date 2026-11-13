# Held-out galaxy sample

The second campaign drew 24 BOSS galaxies from three fields on separate plate visits. The launch sample consisted of 12 galaxies on plate 6497, MJD 56329. The selection used sky positions, plate quality, program, classification, redshift and observation identifiers. The observation IDs and catalog positions were fixed before downloading the held-out spectra. [Frozen scientific protocol](../protocols/HOLDOUT_V1.md) · [Target list and per-file checksums](../data/holdout-targets.csv)

## Metadata-only selection

Starting after plate 6497, scan good BOSS main-program DR17 visits in ascending plate and MJD order. Accept a visit only if its plate and MJD are both new and its center lies at least 8° from the launch center and each earlier accepted center. Apply the launch galaxy selection to its fibers: both class fields must say GALAXY, `sciencePrimary = 1`, `zWarning_noqso = 0`, `0.02 ≤ z_noqso ≤ 0.10`, `survey = programname = 'boss'` and `run2d = 'v5_13_2'`. Require at least six distinct eligible targets. Select at most the first eight by fiber ID. Exclude any position within 2 arcseconds of a launch or previously selected target. Stop at three accepted visits. The metadata queries were bounded to plate IDs 6497–6752 for `PlateX` and 6497–6560 for galaxy fibers; the third visit was reached by plate 6509. The script contains the two exact SQL statements. [Selection and retrieval code](../examples/fetch_holdout.py) · [SDSS catalog definitions](https://www.sdss4.org/dr17/spectro/catalogs/)

| Plate | MJD | Plate center RA, Dec (°) | Eligible before fiber limit | Selected | Center separation from launch (°) |
| ---: | ---: | --- | ---: | ---: | ---: |
| 6498 | 56565 | 359.5821, +31.4156 | 25 | 8 | 111.00 |
| 6503 | 56564 | 349.3870, +31.7280 | 9 | 8 | 105.74 |
| 6509 | 56486 | 338.3943, +31.1752 | 8 | 8 | 99.78 |

The minimum separation between holdout field centers is 8.6885°. The nearest two observed targets from different holdout fields are 7.1889° apart. The nearest launch–holdout target pair is 98.3952° apart. Even within a single selected visit, every two selected positions are separated by at least 164.7 arcseconds. BOSS plates subtend about 3° in diameter. The eight-fiber cap gives each field equal weight here, and ordering by fiber can favor one part of a plate or spectrograph. These data provide a test on new fields for the chosen galaxy population. [SDSS plate geometry](https://www.sdss4.org/dr17/algorithms/boss_tiling/)

## Check for shared exposures

SDSS can combine exposures obtained on more than one night into a single coadd. We inspected the primary FITS headers of all 36 spectra, including all 12 launch spectra. For each `EXPIDnnn` card, we extracted the science-frame number after the camera code, retaining leading zeros. Matching blue and red cameras observing the same integration count once. The union of science integrations for each visit is shown below. [SDSS pipeline and raw-frame indexing](https://www.sdss4.org/dr17/spectro/pipeline/) · [SDSS spectrum and plate identifiers](https://www.sdss4.org/dr17/spectro/spectro_basics/)

| Plate / MJD | Unique science-frame numbers, inclusive | Count |
| --- | --- | ---: |
| 6497 / 56329, launch | 00155803–00155809 | 7 |
| 6498 / 56565 | 00166857–00166859 | 3 |
| 6503 / 56564 | 00166761–00166764 | 4 |
| 6509 / 56486 | 00164755–00164758 | 4 |

These science-frame sets are disjoint across visits. Galaxies on the same plate still share observing conditions. All four fields use the same BOSS instrument family, reduction, survey redshift estimator and detector. Science-frame IDs establish exposure independence across the fields; they do not test common calibration errors. Sky-coordinate selection also means the comparison measures recurrence across different galaxies. [Target-level exposure audit](../data/holdout-targets.csv)

## Reproduce or resume

From the repository root, install the archive extra. The first command checks the live metadata selection against the published IDs and coordinates, then saves both raw catalog replies and the frozen five-column ID list. The fetch command reads the pinned roster directly and can run independently of SkyServer availability.

```sh
python -m pip install -e '.[archive]'
python -m examples.fetch_holdout select --output runs/holdout-catalog
python -m examples.fetch_holdout fetch --output runs/holdout-inputs
```

The fetch command saves original lite FITS files in `raw/`, verifies the source hash and exposure IDs, and verifies the normalized CSV and observation-manifest hashes. The manifest hash binds the redshift, target coordinates, mask convention and release used in the measurement. The 24 original files total 5,241,600 bytes. A rerun checks cached bytes and reconstructs missing files. An explicit FITS directory supports offline regeneration from prior downloads. It should contain filenames such as `6498-56565-0011.fits`.

```sh
python -m examples.fetch_holdout fetch --fits-dir /path/to/pinned/fits --output runs/replayed-inputs
```

| Audit item | SHA-256 for this experiment |
| --- | --- |
| Raw SkyServer plate response | `596fdc650639fd354e3a782c01e74f47e2cbea137119992d1f68fac8fd66ffb2` |
| Raw SkyServer galaxy response | `1efd464d0c216bf964ad32bb61fa75bc13fc843f5b70bca757c7f04f06a8d725` |
| Five-column ID/coordinate selection | `5ec2510cdb1f69b879706995cc8115ff910396a0f5d3661853a01a2f3b3626ae` |
| Released roster with source, derivative, manifest hashes and exposure IDs | `e36f133feb89003343093eba05a15f6f5adea9450a6e5c3c24934b36fd8e47f9` |

A survey service can change formatting on a later replay. The selection command compares the exact ID/coordinate serialization. Preserve the original declared hashes when reporting this experiment. [Data provenance and SDSS attribution](DATA.md)
