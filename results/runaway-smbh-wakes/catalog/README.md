# Supplementary catalogue screen

An independent pass selects entries from the public Q1 MER table with a VIS detection, 0.95 < ellipticity < 1 and semimajor axis at least 7.5 native pixels. The cut uses published catalogue moments. The pixel search can recover fragmented emission that this catalogue selection misses.

The [HATS reader](hats_screen.py) reconciles IRSA's XML partition listing to the publisher's row counts and applies the cuts directly to the public merged Parquet data. The receipt identifies the parent ETag and records SHA-256 for each exact range read or each decoded Arrow column read through the native S3 interface. The latter hash uses PyArrow 26.0.1. The supplementary axial test pairs selected entries with Gaia DR3 stars of G < 16 out to 120 arcsec and compares the result with the perpendicular direction. A geometric match alone is insufficient to establish a diffraction spike. The dedicated two-band cutouts target extragalactic catalogue entries with DET_QUALITY_FLAG = 0, an exploratory subset. The flag was retained in the full list.

The release's [catalogue dictionary](https://euclid.esac.esa.int/dr/q1/dpdd/merdpd/dpcards/mer_finalcatalog.html) and [morphology cookbook](https://euclid.esac.esa.int/dr/q1/dpdd/merdpd/mermorphologycookbook.html) describe the angle in different terms. The photometry measures both the stored angle counterclockwise from pixel +x and the perpendicular orientation. Gaia comparisons use the local TAN Jacobian. The empirical orientation check, the star alignments and the measured VIS/H profiles share the same catalogue selection.

Run the following from this directory. The Q1 tile and product lists are in the sibling euclid directory. Each script supports a cache outside the repository. The surface-brightness completeness measurements in the detector directory apply to the independent pixel search.

```sh
python hats_screen.py --cache D:/new-campaign/partitions --out D:/new-campaign/results
WAKE_GAIA_CACHE=/path/to/gaia python gaia_stars.py
python star_alignment.py --stars /path/to/gaia
python targeted.py --cache /path/to/vis-h-cutouts
```

A new full catalogue acquisition requires distinct, initially absent external D: cache and output directories. Existing partition JSON is never reused, and generated tables are written to the new output, preserving the deposited historical receipts. The compact selected rows and partition totals were checked offline; the full native catalogue columns were not independently replayed.
