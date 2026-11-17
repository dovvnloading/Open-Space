# Spectroscopic check of the three recovered D6 controls

A common-velocity fit to the Ca II resonance doublet places D6-1 at +1,193 km/s in the heliocentric frame. The same fit yields +71 km/s for D6-2 and +108 km/s for D6-3. The latter two values are uncertain at the ~100 km/s scale. Correlation of the remaining 4,050–6,700 Å spectrum against D6-2 independently gives a D6-1 shift of +1,120 km/s and a D6-3 shift of −30 km/s. Thus the large D6-1 line-of-sight speed and the large contrast with both other stars are visible in the released fluxes.

| Control | Data | Ca II RV, km/s | Residual-bootstrap 16th–84th percentiles | ΔRV from D6-2 outside Ca, km/s |
|---|---|---:|---:|---:|
| D6-1 | SALT/RSS | +1,193 | +1,169 to +1,218 | +1,120 |
| D6-2 | NOT/ALFOSC | +71 | −6 to +170 | template |
| D6-3 | NOT/ALFOSC | +108 | +11 to +173 | −30 |

The inputs are the three reduced spectra released by Shen et al. Each file explicitly labels its air wavelengths as heliocentric. The wavelengths still contain the star's line-of-sight shift. We fit them as supplied, without applying Earth's velocity again. The paper's rest-frame overlay for D6-1 is a separate plotted representation. The reported resolution is R ≈ 900 for D6-1 and FWHM ≈ 15 Å for D6-2/3. [Shen et al., Section 4 and ancillary spectra](https://arxiv.org/abs/1804.11163)

## Fits and limitations

1. Fit the available pixels from 3,912 to 4,024 Å with a linear continuum and two Gaussian absorption profiles at NIST Ca II K and H laboratory wavelengths, 3,933.6614 and 3,968.4673 Å. Both lines share one relativistic Doppler factor and width. Depths vary independently. The common width can exceed the instrumental resolution. Multiple velocity starts search the bounded range −1,600 to +2,500 km/s. Each file provides only wavelength and flux, so pixels carry equal weight. Blockwise residual sign resampling estimates fit uncertainty on the scale of one instrumental resolution element. [NIST Ca II wavelengths](https://physics.nist.gov/PhysRefData/Handbook/Tables/calciumtable4.htm)
2. For a separate relative-shift test, normalize to a 40 Å Gaussian-smoothed continuum, approximately match SALT to the NOT resolution, then correlate D6-1 and D6-3 against the **observed D6-2 spectrum**. Exclude the Ca II fit region, the [O I] 5,577 Å sky feature, Na D, and the stronger red sky/telluric regions. The script gives the eight included bands. It scans velocities at 10 km/s spacing and weights each band's Pearson correlation equally. Omitting one band at a time gives +1,110 to +1,140 km/s for D6-1 and −30 km/s in every trial for D6-3. The grid and shared data prevent interpreting that zero spread as an error bar.
3. With the Ca-fit Doppler factor removed, measure signed equivalent widths relative to a local linear continuum in the Mg II 4,481 Å, O I 7,774 Å, Balmer and He I 6,678 Å regions. The oxygen region lies outside the SALT coverage. The D6-2/3 oxygen regions show 8.07/3.95 Å of integrated absorption, and the Mg-region indices span 2.57–3.15 Å across all three. None shows a strong Hα or Hβ absorption index with these windows. These are blend-sensitive band measurements, tabulated in feature_indices.csv. Low resolution, continuum placement, and the lack of a metal-rich atmosphere model preclude elemental abundance estimates or He upper limits from these indices. The cooler stars can hide substantial atmospheric He.

Wavelength-calibration errors, photospheric blends, unresolved interstellar calcium, and gravitational redshift are outside the bootstrap. The physical lower width bound is active in both NOT Ca fits. For D6-3 the doublet fit and the full-spectrum shift differ by about 70 km/s when placed on D6-2's Ca velocity scale. The +108 km/s measurement therefore does not supersede the published −20 ± 80 km/s. Published reference velocities in control/published_kinematics.csv are the comparison values for the orbit work.

## Reproduction

Run `OPENBLAS_NUM_THREADS=1 python measure.py` with Python, NumPy ≥2.0 and SciPy. Each download verifies the release file's wavelength-frame header and matches the committed SHA-256 receipt. A byte-identical local copy can be used with `--cache DIRECTORY`. Both paths produce line_measurements.csv and feature_indices.csv.
