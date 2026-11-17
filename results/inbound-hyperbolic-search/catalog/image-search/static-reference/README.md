# Historical PS1 reference images at the ZTF aperture positions

The Pan-STARRS r deep stacks show little quiescent light in the six-arcsecond apertures associated with the October 25 and 27 ZTF signals. For the October 25 frame, integrating the PS1 pixels gives a baseline equivalent to about 80 ± 40 ZTF DN compared with the measured 1,138 DN in that frame. The comparison supplies no bright persistent source that accounts for all of the ZTF flux. A transient, variable or blended source, sky subtraction error, and a noise excursion remain possible. The moving-comet interpretation still lacks repeatable astrometry.

| ZTF frame ID | UT date | PS1 r fixed-aperture counts s⁻¹ | Empirical off-sky MAD scale (counts s⁻¹) | Approximate equivalent ZTF DN |
| --- | --- | ---: | ---: | ---: |
| 20251025541736 | Oct 25 | 30.43 | 15.35 | 79.9 ± 40.3 |
| 20251027542164 | Oct 27 | 4.60 | 15.43 | 13.1 ± 43.9 |
| 20251027547164 | Oct 27 | 7.59 | 16.80 | 21.9 ± 48.5 |
| 20251029546134 | Oct 29, known control | 21.86 | 16.08 | 70.6 ± 51.9 |

Four 720 × 720 FITS cutouts were measured at their native 0.25 arcsec scale. The table uses a 6.1 arcsec radius and an 11.2–15.3 arcsec median sky annulus. The 48 shifted apertures per cutout give a local robust noise scale, which includes sky structure. We use the MAST FITS-cutout linear counts, exposure time and published stack calibration, then translate to the science frame's ZTF zeropoint assuming the same r-band spectral response. The uncertainties shown are spatial noise. The conversion omits color terms, PSF matching, variability and differential aperture corrections. Weak positive fluxes in the machine-readable file are noisy measurements and should not be treated as magnitude detections.

The PS1 images predate the event. They test the stationary background at these positions, independently of the October 2025 ZTF exposures. [Four FITS and query receipts](results.json) · [Measurement code](../../ps1_static_image.py) · [PS1 stack calibration](https://outerspace.stsci.edu/spaces/PANSTARRS/pages/298812339/PS1%2BStack%2BImages)
