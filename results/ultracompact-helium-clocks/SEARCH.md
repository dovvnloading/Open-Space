# Search for short photometric clocks in helium-accreting binaries

Select the 2026 Green catalogue of mass-transferring ultracompact binaries (Zenodo record 18842123). Search public ZTF photometry of systems lacking a tabulated measured period. Keep the input target denominator and each coverage failure. Candidate frequencies must survive temporal or passband holdout, the daily sampling window, and a target-specific check for previous reports. An outburst recurrence or superhump alone does not determine the orbital period.

The catalogue has 153 entries. This checkpoint records the research question before inspecting any light-curve periodograms.

[Catalog release](https://zenodo.org/records/18842123) · [IRSA archive](https://irsa.ipac.caltech.edu/Missions/ztf.html)

## Outcome of the initial ZTF screen

Of 51 missing-period entries, 49 were north of declination −32°. The initial archive pass obtained nine light-curve responses, including three zero-row responses. Three targets had at least 75 usable detections after quality selection. For each we fit trial sinusoids over 18–144 cycles/day on the first 65% of the photometry, then used their phases and amplitudes to predict the remainder. None of the eight strongest training peaks improved the held-out weighted residual variance. The other 40 target requests timed out on three attempts each. Of two additional known-period controls, AM CVn returned only 50 raw measurements and CR Boo timed out. The independently successful experiment used TESS images and spatially resolved ZTF alert photometry of the eruption. [Initial query URLs, data hashes, coverage and all peaks](ztf-screen.json) · [Measured result](README.md)
