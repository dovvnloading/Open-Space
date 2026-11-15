# HD 154615 before TESS

ASAS-3 observed HD 154615 from 2001 through 2009. Its catalog returns 628 V-band measurements in four fields. We shifted a two-dip template fixed from the TESS light curves across one full orbital cycle, fitting one magnitude offset per ASAS field and a common depth scale. In the first-listed aperture, 579 grade A/B measurements with magnitudes between 9 and 10 give a shift of **−50.9 minutes** relative to BTJD 1644.95 + 1.851076 × cycle. Independent chronological halves give −47.9 and −53.4 minutes. All five aperture point estimates lie between −57.0 and −47.9 minutes. MAG_0 also admits a disconnected alternative approximately half an orbit away at the nominal 95% threshold. The other four apertures favor a single phase interval. On 58 Gaia G-band measurements from 2014–2017, the same fixed photometric template gives **−30.2 minutes**. The main ASAS and Gaia profile intervals are −70.1 to −31.3 and −39.8 to −20.7 minutes. [Original ASAS measurements](inputs/asas_170827-4207.5.txt) · [Epoch shifts and aperture alternatives](early_phase.csv) · [Gaia epoch photometry](../inputs/gaia_epoch_5965875226751271040.csv)

These are best-fit mean shifts over each survey, conditional on the template and the displayed TESS period. A small change in that period has an increasingly large effect on the oldest residuals. The six-year TESS curvature reported above is invariant to that choice. Fits use a full-cycle search, positive eclipse depth and a profile sum-of-squares interval (Δχ² = 3.84 with one residual variance per fit). The ASAS HJD time stamps and Gaia BJD(TCB) were aligned by their published Julian-day origins, without converting between heliocentric and barycentric arrival times or time scales, which can shift epoch comparisons by about a minute. ASAS recorded these observations with 14-arcsecond pixels and multiple apertures. The historical series establish the persistence of the eclipses and give earlier timing constraints, with depth, blending and instrument-dependent systematics. [ASAS catalog and data credits](https://www.astrouw.edu.pl/asas/?page=aasc)

## Alias and cycle-count limits

The grid spans one complete orbital cycle (±0.925538 days) and is periodic at its boundary. The MAG_0 alternative corresponds to swapping the primary and secondary eclipse. The bounds in the CSV enclose both disconnected intervals when `confidence_segments=2`. No fit tests an absolute integer number of cycles across the long gaps. Moving the oldest ASAS series by a full cycle relative to the TESS epoch requires a period change of order 0.0005–0.0007 d across the corresponding baseline. Each half of the ASAS sample independently favors the same small phase displacement at the adopted period. An unknown long-period modulation or a sufficiently different past period can change the inferred cycle count. A period error of only 0.4 seconds would shift an epoch 15 years away by about 20 minutes.

The ASAS file header declares four field counts totaling 628, matching the included observations. The listed `median_time`, `earliest` and `latest` are Julian-day offsets from 2457000 using each survey’s original time convention. Run the analysis from this directory with Python and NumPy.

```sh
sha256sum -c inputs.sha256
python early_photometry.py
```

[Source queries and original file hashes](sources.json)
