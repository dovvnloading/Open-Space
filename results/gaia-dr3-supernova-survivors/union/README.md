# Combined spectroscopic catalogue limit

**Zero of 67,433 unique objects with qualifying Gaia astrometry have a measured heliocentric transverse speed above 1,000 km/s.** The threshold is a corrected parallax signal-to-noise ratio of at least five. Lowering that cut to three retains 75,684 objects and still yields zero. The largest point estimate at either cut is 704.56 km/s for Gaia DR3 3688712561723372672, classified as a DA white dwarf in the El-Badry et al. follow-up. [El-Badry et al. 2026, Table 4](https://arxiv.org/abs/2606.11293)

| Spectroscopic product | Unique IDs | Eligible at 3σ | Eligible at 5σ | Maximum at 5σ (km/s) |
|---|---:|---:|---:|---:|
| DESI DR1 visually inspected WD-candidate catalogue | 63,152 | 37,832 | 33,322 | 704.56 |
| SDSS-V DR19 classified SnowWhite entries | 10,364 | 8,541 | 7,843 | 486.90 |
| SDSS-V DR20 classified SnowWhite entries | 50,862 | 40,073 | 36,400 | 601.97 |
| LAMOST DR12 v2 low-resolution WD table | 13,850 | 9,972 | 8,956 | 457.86 |
| **Unique sources across the four products** | **113,166** | **75,684** | **67,433** | **704.56** |

The rows overlap. DESI and the combined SDSS catalogue share 11,196 source IDs. LAMOST shares 4,019 with their union and adds 9,831. DR19 and DR20 share 9,847. Every source is counted exactly once in the combined denominator. overlap.csv contains all 15 mutually exclusive membership combinations before and after selection. The mask adds 1 for DESI, 2 for DR19, 4 for DR20 and 8 for LAMOST. The optional membership CSV contains each source's mask and 3σ/5σ flags. Its sorted-content SHA-256 is recorded in summary.json.

Of the 113,166 distinct sources, 113,028 have a 5- or 6-parameter solution, finite parallax and proper motions, and positive parallax uncertainty. The individual Lindegren zero point is defined for 109,391. The quality criteria retain 95,082 of the calibrated solutions before any S/N gate. For each, the corrected parallax is catalogue parallax minus the source's zero point. We use its ratio to the reported parallax error and calculate 4.74047 times total proper motion divided by corrected parallax (mas/yr and mas). The five quality conditions require 0 < RUWE < 1.4, at least nine visibility periods, no duplicated-source flag, IPD multi-peak fraction ≤2%, and IPD goodness-of-fit harmonic amplitude <0.1.

Two of the twenty published control objects are in these WD-oriented spectral products. J1240+6710 passes both astrometric gates with a point estimate of 422.6 km/s. J1109+0001 falls below both at corrected S/N 1.32, even though its formal transverse speed is 1,026.1 km/s. None of the eleven D6 controls or US 708 appears in this product union. Crossmatches in other spectral products, including legacy survey spectra for some of these objects, are outside this denominator. DESI's parent explicitly includes programme contaminants. Neither the parent count nor the selected count is a census of confirmed white dwarfs.

The finite-catalogue result concerns **measured tangential speeds**. Under an additional hypothesis of independent, comparable draws from a population of eligible catalogue measurements, the one-sided exact binomial 95% upper fraction is 1 − 0.05^(1/N), or 44.4 per million at 5σ (39.6 per million at 3σ). This conditional calculation cannot constrain the fraction of physical supernova survivors without the survey targeting, class-assignment and kinematic selection functions. Objects dominated by radial motion can pass undetected at either threshold. Requiring noisy parallaxes to meet a measured S/N cut can further bias the sample.

## Reproduction and validation

First run the sibling desi/fetch.py and sdss/screen.py on the same cache, and spectra/lamost_parent.py for the verified LAMOST file. Then run:

```sh
python census.py /path/to/cache /path/to/lamost_wd.csv.gz --membership-out /path/to/membership.csv
```

The script verifies every previous Gaia batch hash, its full ID list and its returned columns. SDSS membership uses the classifications and exact IDs in the two already extracted FITS files, with archive SHA-256 and gzip CRC/size validation. The LAMOST archive is checked against its published-input receipt. The 45 new Gaia batches contain only the 9,831 missing LAMOST IDs. They have the same columns and ordering as the DESI batches. Their input-ID, exact ADQL and returned-data hashes are in lamost_gaia_receipt.json. Membership is checked against the published release-by-release counts from the upstream analysis. controls.csv keeps missing controls explicit. top10.csv gives the fastest ten eligible objects at both thresholds.
