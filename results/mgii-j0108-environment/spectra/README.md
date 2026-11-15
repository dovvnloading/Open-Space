# Repeat spectroscopy of the three background quasars

The DESI Mg II absorber catalog supplies the three redshifts. Holding each sightline's redshift fixed, we measured both lines directly in ten spectra, including seven earlier SDSS spectra. The 2010 October 6 SDSS plate alone recovers both lines in every sightline. The two 2001 spectra also recover both lines. Comparison to SDSS specifically tests whether the DESI detection repeats in different hardware at a different epoch.

| Sightline | DESI rest EW 2796 / 2803 (Å) | SDSS plate 4226 rest EW 2796 / 2803 (Å) |
|---|---|---|
| A | 1.36 ± 0.09 / 1.13 ± 0.09 | 1.42 ± 0.14 / 1.36 ± 0.14 |
| B | 1.68 ± 0.06 / 1.40 ± 0.06 | 1.54 ± 0.13 / 1.23 ± 0.14 |
| C | 1.40 ± 0.17 / 1.05 ± 0.18 | 1.27 ± 0.17 / 1.38 ± 0.16 |

For each doublet, the default estimator fits a weighted linear continuum within ±40 observed Å of its midpoint, excludes ±8 Å about all specified transitions and clips sidebands once at 4σ. It integrates fractional pixel widths within ±210 km/s of each predicted line center. An error propagates the pipeline inverse variance and the covariance of the fitted continuum parameters under a diagonal pixel-noise approximation. A line touching a masked pixel is omitted. Narrower and wider windows and a constant continuum provide sensitivity checks. These conditional signal-to-noise ratios describe the pixels at known redshifts. Resampling correlations, continuum systematics, spectral blending and survey selection affect their interpretation.

The A and B quasars lie at z > 3.2. Their shortest-wavelength iron-line windows fall in the Lyα forest and contain blends. Their strong isolated apparent iron features cannot establish an additional ionic identification. Sightline C places the iron windows redward of the quasar's Lyα line.

`python measure.py` regenerates `line_measurements.csv` using NumPy and the included flux, variance-weight and mask cutouts under `pixels/`. To regenerate the cutouts and compare their checksums with the committed files, run `python extract.py /path/to/cache` from this directory with Astropy and NumPy installed. The largest download is the 732 MB DESI coadd. SDSS plates 694 and 3736 observed A and B in October 2001 and January 2010, respectively. `sources.csv` records exact products and SHA-256 hashes. `visits.csv` records dates, coordinates and fiber IDs. DESI exposures shared across objects still use separate fibers.

[DESI DR1 spectroscopy](https://data.desi.lbl.gov/doc/organization/) · [DESI Mg II absorber catalog and data model](https://data.desi.lbl.gov/doc/releases/dr1/vac/mgii-absorber/) · [SDSS DR17 optical spectra](https://www.sdss4.org/dr17/spectro/spectro_access/)
