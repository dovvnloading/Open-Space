# CFHT photometry of 5XMM J085618.5+380711

The deep public MegaPrime images cover the full XMM localization. The nearest ≥5σ optical peak appears in both g and r about 6.6 arcsec from the X-ray position. Both the native and locally boresight-corrected 3.204 arcsec regions lie clear of this feature. The r-band peak has coordinates (134.0794405, +38.1192415) deg. Its separation is 6.55 arcsec from the native centroid and 6.60 arcsec from the shifted centroid.

| CFHT filter | Header exposure | Largest aperture signal anywhere in either region | Nominal local 5σ sensitivity, AB | Conservative AB limit over both regions |
| --- | ---: | ---: | ---: | ---: |
| g.MP9401 | 5,860 s | 1.20σ | 25.33 | >25.09 |
| r.MP9601 | 25,351 s | 1.72σ | 25.16 | >24.84 |

The last column adds five noise standard deviations to the largest measured positive aperture signal over both regions before correcting to point-source flux. The depth is conditional on the local noise estimate and a stationary, unresolved object. With the observed DR15 0.2–2 keV band-flux sum of 1.6005 × 10⁻¹⁴ erg cm⁻² s⁻¹, the conservative limit yields F_X/(ν f_ν,r) > 7.9 for a representative r wavelength of 6231 Å. Using nominal sensitivity gives 10.6. This convention uses monochromatic ν f_ν. Published ratios defined using an integrated optical passband have a different zero point. Spectral model uncertainty in the X-ray flux is additional.

## Photometry and checks

We downloaded calibrated CADC MegaPipe.211.256 g and r science cutouts and the corresponding weight planes. The g stack includes six exposures and the r stack twelve. The header MJD-OBS records the start of the images, February 2006 for each stack. Total exposure in a stack header need not equal the exposure at a specific pixel. Proper motion after the October 2005 XMM exposure and over a stack's constituent dates remains unconstrained. The public FITS headers specify AB magnitude = PHOTZP − 2.5 log10(summed flux) with PHOTZP = 30. [CADC stacking and calibration](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/proc.html)

Each cutout spans 545 × 545 pixels at 0.18577 arcsec per pixel. We measured the field at (134.077254, +38.119831) deg and at the position after applying the independently estimated (+0.347, +0.949) arcsec east/north shift. Both circles, including every pixel under the 0.743 arcsec radius photometry apertures, have fully positive uniform weights, finite science values and counts below the saturation threshold. The two masks have identical WCS grids to the corresponding science images. The light from the nearby peak did not produce a significant aperture anywhere in either search circle.

The 49-pixel digital aperture sums on a ten-pixel grid measure noise on the same scale as the source extraction, including correlated noise from resampling. A 10–40 arcsec annulus provides 1,368 nonoverlapping apertures. The median aperture flux defines the sky. We used the largest of the four quadrant MAD-based noise estimates, each scaled by 1.4826. Their range is 7.13–11.97 image flux units in g and 8.20–13.61 in r. Positive field sources are retained for this robust estimate. This is a sensitivity estimate, not an injection/recovery completeness measurement.

Two unsaturated local Legacy Survey PSF sources have 15-pixel-radius CFHT magnitudes within 0.10 mag of their Legacy photometry in both bands. The smaller of the two 4-pixel to 15-pixel flux fractions was adopted in each band, 0.809 in g and 0.791 in r. This calibrates the effective area of the compact aperture to a 2.79 arcsec radius. Unknown light beyond the reference aperture and color-dependent differences between the Legacy and MegaCam filters limit the absolute comparison. At the native and shifted X-ray centers the forced photometry has S/N of +0.53 and +0.28 in g, and −0.38 and −0.45 in r. A search over all integer pixel centers in both error regions produced the maxima tabulated above.

The replay requires NumPy, SciPy and Astropy. Its four URLs and both raw-file and pixel-array SHA-256 hashes appear in the JSON. CADC can regenerate header creation timestamps for cutouts. The pixel digest removes those header-only differences. The nearest detected optical peaks, depth estimates, calibrator counts, offsets and weight limits can be checked directly in the receipt.

```sh
python analysis/cfht_photometry.py cfht_cache results/cfht-photometry.json
```
