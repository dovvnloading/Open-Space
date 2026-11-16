# Relative motion of (556068) 2014 JR80 in the 2026 Rubin alert stream

ALeRCE’s moving-object alerts identify five visits in which Rubin resolves the changing orientation of the nearly equal-brightness trans-Neptunian binary 2014 JR80. The *i*-band axis turns **49.46° in 11.963 days**, from June 27 to July 9 UTC. Hubble images an additional 13.7 days earlier give an independent orientation. The binary discovery belongs to Proudfoot et al.; their September 2026 CBET reported its separation and orientation with HST and six earlier CFHT/Blanco epochs spanning separations of 0.5–1.0″. The new measurement uses Rubin’s dated relative positions to constrain the mutual motion. [Discovery circular](http://www.cbat.eps.harvard.edu/iau/cbet/005700/CBET005739.txt) · [ALeRCE public alerts](https://science.alerce.online/alerce-lsst/)

| UTC date in 2026 | Instrument and band | Independent exposures | Separation | Position angle east of north, modulo 180° |
| --- | --- | ---: | ---: | ---: |
| June 13 | HST/WFC3 F350LP | 4 | 0.8461″ | 15.98° |
| June 27 | Rubin LSSTCam *i*, *z* | 2 | 0.633″, 0.676″ | 55.77°, 56.50° |
| July 7 | Rubin LSSTCam *i*, *z* | 2 | 0.541″, 0.597″ | 94.91°, 93.55° |
| July 9 | Rubin LSSTCam *i* | 1 | 0.584″ | 105.23° |

Rubin exposures use the ALeRCE identifiers and ephemerides; the matching cutouts were recovered from Fink’s copy of the same Rubin packets. Three observing nights and four HST exposures within one visit supply four 2026 epochs. Because the components are similar in brightness, every geometric or orbital inference must allow a 180° component swap.

![Resolved Rubin cutouts, single-source residuals and fitted relative positions](output/rubin_astrometry.png)

**Physical constraint.** For a bound pair, the earliest Rubin *i* measurement with the June HST position sets a system-mass floor of **1.71 × 10¹⁹ kg**, conditional on relative-position error envelopes of 0.02″ for HST and 0.06″ for Rubin. Enlarging both envelopes to 0.10″ gives 8.33 × 10¹⁸ kg. This finite-time estimate permits either component assignment and unconstrained depths in the tangent-plane projection. It includes distance and light-time sensitivity. The error radii represent explicit assumptions without a probability calibration. [Derivation, code and every epoch pair](dynamics/README.md)

**What limits the image-only interpretation.** The telescope rotator moved about 56° during the interval, nearly tracking the change in celestial position angle. A fixed distortion in detector coordinates is a serious alternative using Rubin alone. The nominal cutout PSF, alternate fitting radius, broadened PSF and difference images all give comparable turns. The science and difference products share the same photons. The likelihood improvements for a two-source fit also reach large values in some other moving objects. The strongest test uses independent telescope data and same-detector reference stars. [Rubin measurements and controls](METHODS.md) · [Hubble remeasurement](hst/README.md)

**Reproduction.** Start with the frozen ALeRCE queries, source and peer tables, per-exposure cutout hashes and the commands in the methods. The HST measurement and mass constraint each have an offline reconstruction. Archive identifiers, timestamps, pixel conventions, selection and controls accompany the respective products. [Selection and reproduction](METHODS.md) · [Other questions tested in the same database](../rubin-alerce-2026/README.md)
