# Matching oxygen-line wings in DESI 39628534716566336

In a DESI DR1 spectrum of the galaxy at redshift 0.841615, the bright [O III] line has more flux in its wings than one Gaussian predicts after convolution with the published instrumental response. A profile fitted to that line also improves the fit to the weaker [O III] line 88 Å away in the observed spectrum. The previously published Open-Space analysis identified this galaxy and its projection near a lower-redshift galaxy, established the two apertures’ different observed widths, and fitted individual Gaussian profiles. This analysis tests the line shape against the instrument, transfers it between transitions, and checks sky fibers from the original exposure. It is a measurement of a known galaxy. Priority in the wider literature is unestablished. [Project investigation](https://github.com/dovvnloading/Open-Space/blob/5f8bc79e188ff1940d87a59e43fe2e96a9083059/docs/SDSS_SUPERPOSITION.md)

## Observations and model

The offset (northwest) and foreground-centered apertures are separated by about 3.86 arcsec. Their DESI target IDs are 39628534716566336 and 39628534716566366, respectively. The dark-program exposure 108405 and bright-program exposure 108404 were taken on observing night 20211111 through the same repositioned fiber 3291. Each coadd contains one exposure. The background-galaxy redshift used for all velocities is 0.8416150156358213. In the second aperture, the much brighter foreground galaxy is at redshift 0.421496. See provenance.json for archived products and their SHA-256 hashes. [DESI DR1](https://data.desi.lbl.gov/doc/releases/dr1/)

Use the directly downloaded Z-camera data, positive pipeline inverse variance, a zero pixel mask, and all qualifying pixels within ±680 km/s of the expected rest-vacuum wavelengths 5008.239 and 4960.295 Å. Model flux over each 0.8 Å pixel, then apply the target-, wavelength-, and aperture-specific DESI resolution matrix. Compare a single Gaussian against a shared-centroid sum of two Gaussians, both with a linear continuum. Nuisance amplitudes are fit by weighted least squares. The offset training fit returned positive flux for both components. The velocity centroid, both widths, and wing fraction from each strong-line model are frozen when predicting its companion. Each companion model has three freely fitted coefficients, one intensity and two for continuum. Both transitions were examined during exploratory analysis, so this is a cross-line consistency check, not a formally blinded observation. [Instrumental response convention](https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/tiles/GROUPTYPE/TILEID/GROUPID/coadd-SPECTROGRAPH-TILEID-GROUPID.html)

| Measured quantity | Result |
| --- | --- |
| Offset strong line, single Gaussian and two Gaussians | χ² = 133.65 and 51.49 for 53 pixels and 5 and 7 fitted parameters |
| Offset companion, frozen single and two-component shapes | χ² = 60.06 and 46.84 for 52 pixels and 3 coefficients in both fits |
| Model core and wing σ before instrumental convolution | 12.3 and 61.4 km/s |
| Wing fraction of total modeled 5008 flux | 20.4% |
| Strong-to-weak line flux ratio with transferred shape | 2.965 |
| Same profile transferred to foreground-centered 5008 line | χ² = 90.92 with free shift and normalization plus linear continuum, 46.33 with an additional 57.1 ± 8.5 km/s of Gaussian broadening |

The foreground-centered 4960 line alone gives only a 1.81 improvement in χ² on adding broadening. Across half-widths of 430, 520, 680, 820 km/s and constant, linear, and quadratic continua, all 12 companion-line comparisons favor the two-component prediction. Their improvements range from 3.93 to 13.92. The inferred strong-line wing fraction ranges from 19.9% to 20.5%. The finite-pixel core is near the sampling limit, so its fitted width is especially sensitive to modeling. The atomic intensity ratio of the pair is approximately 2.98. [Storey & Zeippen 2000](https://doi.org/10.1046/j.1365-8711.2000.03184.x)

There are 49 designated sky fibers from the same exposure and camera. Twelve have zero fiber status and pass the local pixel coverage requirements in both windows. Fitting continuum, core and wing at the target line shapes, their largest positive wing coefficients divided by their formal errors are 1.13 and 1.21. Passing this control disfavors additive blank-sky residuals. Calibration or line-spread errors shared by object light can survive this test. The two oxygen lines used the same instrument on the same exposure. The formal errors assume diagonal pixel noise and a fixed instrumental response. A second spectrograph with well-measured response or spatially resolved spectroscopy can test the source of the wings and the aperture-width difference. These data do not determine an outflow, a second gas component, or gravitational lensing.

The published FastSpecFit DR1 v3.0 catalog contains both oxygen-line measurements for the offset galaxy. Its OII_DOUBLET_RATIO is 0.7019, and it assigns OIII_5007_SIGMA = 20.37 km/s. This analysis adds a cross-line test of the profile of the cataloged doublet. [FastSpecFit catalog](https://data.desi.lbl.gov/doc/releases/dr1/vac/fastspecfit/)

## Reproduce

With Python 3.11 or later and NumPy, SciPy and Matplotlib installed, run the following command in this directory. It reads windows.npz and regenerates results.json and oxygen-profile.png offline.

```sh
python analysis.py
```

To independently regenerate the compact extraction from the public FITS files, install Astropy and run the following commands. The first command downloads approximately 405 MB, verifies exact bytes, and rebuilds the cutouts and source metadata. The cutouts contain the two target spectra and resolution matrices, and measured blank-sky fluxes, masks, and inverse variances. The plot directly displays the observed DESI flux and residuals. The curves and Gaussian widths are model-derived. Core and broad-component widths are Gaussian standard deviations; all flux ratios are dimensionless.

```sh
python extract.py --cache fits-cache
python analysis.py
```

## Data credits

The spectra, resolution matrices and sky-fiber products are from DESI DR1. This reanalysis and its compact data cutouts derive from those products. Use DESI's data license and acknowledgment when redistributing or publishing the measurements. [DR1 paper](https://doi.org/10.3847/1538-3881/ae4c43) · [DESI acknowledgments](https://data.desi.lbl.gov/doc/acknowledgments/) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

This research used data obtained with the Dark Energy Spectroscopic Instrument (DESI). DESI construction and operations is managed by the Lawrence Berkeley National Laboratory. This material is based upon work supported by the U.S. Department of Energy, Office of Science, Office of High-Energy Physics, under Contract No. DE–AC02–05CH11231, and by the National Energy Research Scientific Computing Center, a DOE Office of Science User Facility under the same contract. Additional support for DESI was provided by the U.S. National Science Foundation (NSF), Division of Astronomical Sciences under Contract No. AST-0950945 to the NSF’s National Optical-Infrared Astronomy Research Laboratory; the Science and Technology Facilities Council of the United Kingdom; the Gordon and Betty Moore Foundation; the Heising-Simons Foundation; the French Alternative Energies and Atomic Energy Commission (CEA); the National Council of Humanities, Science and Technology of Mexico (CONAHCYT); the Ministry of Science and Innovation of Spain (MICINN), and by the DESI Member Institutions: www.desi.lbl.gov/collaborating-institutions. The DESI collaboration is honored to be permitted to conduct scientific research on I’oligam Du’ag (Kitt Peak), a mountain with particular significance to the Tohono O’odham Nation. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the U.S. National Science Foundation, the U.S. Department of Energy, or any of the listed funding agencies.
