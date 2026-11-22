# Review of individual Chandra triggers

## ObsID 26248

Two real triggers passed the instrument cuts in this bright Galactic remnant field. The inspected images and initial counterpart receipts leave both events without an established extragalactic origin. Both have a measured 0.5–7 keV aperture peak flux. The first catalog run recorded two service failures. Corrected catalog classifications and independent visits must be read alongside these notes.

| Trigger | J2000 RA, Dec (deg) | MJD UTC | Aperture t95−t05 (s) | 60 s peak flux ± counting sigma (erg cm⁻² s⁻¹) |
|---|---|---:|---:|---:|
| real/7 | 350.994036, +58.817796 | 59605.079441 | 164.82 | (1.55 ± 0.84) × 10⁻¹² |
| real/8 | 350.881778, +58.768756 | 59604.853971 | 47.14 | (2.68 ± 0.94) × 10⁻¹² |

**real/7.** The aperture at Galactic latitude −2.15° lies outside the bright shell in the full-field image. Nine photons occupy nine detector pixels and eight frames during the 180 s search window, in two groups near its ends. Their changing detector coordinates broadly follow the field-derived dither trace. Two Gaia objects with significant stellar astrometry lie 2.11″ and 2.43″ from the centroid. Each has a local chance-projection upper estimate of 0.11, above the 0.01 threshold for an individual foreground identification. A stellar flare or a fluctuation on the local background remains plausible. [Detector image and light curves](evidence/history/44981778618/chandra-evidence-2/photometry/26248/acisf26248N002_evt2.fits.gz/real/7/diagnostics.png)

**real/8.** The aperture intersects structured emission at the shell rim at Galactic latitude −2.18°. Eleven photons occupy eleven detector pixels and ten frames during the 60 s window. The aperture light curve shows comparable excursions elsewhere in the visit, and the long-exposure image shows a brightness gradient through the background annulus. SIMBAD lists molecular cloud MWISP G111.734−02.178 at a nominal separation of 3.86″ with a 5″ major-axis catalog uncertainty. Diffuse Galactic emission and chance clustering remain viable explanations. No secure counterpart follows from the cloud position. [Detector image and light curves](evidence/history/44981778618/chandra-evidence-2/photometry/26248/acisf26248N002_evt2.fits.gz/real/8/diagnostics.png)

The event input is `acisf26248N002_evt2.fits.gz`, 108,058,563 bytes, SHA-256 `bc566c8eca0764a38c62819e4708c966723382ffdc1f50caf3b1675082ce8f6e`. The aperture-flux uncertainties include counting statistics. The original observation ledger records each match, rejected trigger, control and injection. [Event input](https://cxc.cfa.harvard.edu/cdaftp/byobsid/8/26248/primary/acisf26248N002_evt2.fits.gz) · [Original evidence artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44981778618/artifacts/14739004080) · [Archived event bytes](https://github.com/dovvnloading/Open-Space/actions/runs/44981778618/artifacts/14739223358)
