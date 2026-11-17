# ZTF difference-image search on C/2025 V1

The archived difference frames yielded no validated precovery on October 26 or 27. This exploratory run rejects nonfinite pixels but does not ingest the ZTF image-mask products; finite flagged pixels and artifacts can affect both target and controls. The tabulated aperture depths are unvalidated noise scales, not calibrated limiting magnitudes. October 27 had the highest off-arc signal, 3.67 times the local robust noise in a 3.5-pixel aperture. Eight of 168 translated tracks produced at least as large a maximum when searched over the same uncertainty ellipse. Two of three difference frames contributed to the signal. A fixed sky location also retained positive residual flux. We carried no new astrometric position into the orbit fit.

| UT date in 2025 | Loaded / targeted frames | Trial sky offsets | Robust S/N at nominal position | Robust S/N at best offset | Local maxima ≥ target | Rank p | 5σ aperture flux at zero point 25 | Approximate r aperture magnitude at 5σ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Oct 18 | 0 / 2 | — | — | — | — | — | — | — |
| Oct 25 | 0 / 3 | — | — | — | — | — | — | — |
| Oct 26 | 2 / 2 | 69 | 1.23 | 1.62 | 48 / 142 | 0.343 | 350.3 DN | 18.64 |
| Oct 27 | 3 / 3 | 55 | 2.11 | 3.67 | 8 / 168 | 0.0533 | 204.9 DN | 19.22 |
| Oct 29, published control | 2 / 2 | 29 | 6.06 | 6.26 | 9 / 168 | 0.0592 | 328.7 DN | 18.71 |

All five missing requests on October 18 and 25 returned HTTP 404 for the cutout and the full difference quadrant. The [individual receipts](stack-receipts.json) identify both requests per missing frame. The known October 29 comparison was selected as a known-source control; it is not an independently blinded validation. It peaks one arcsecond west of the nominal ephemeris. Most of the signal comes from one of the two exposures. The broad upper tails of the control maxima, even at this previously measured position, limit the diagnostic power of spatial peak rank.

## Measurement

[stack.py](stack.py) opens the IRSA science-minus-reference FITS, uses its WCS and shutter epoch, and retrieves a fresh I41 astrometric ephemeris for the exposure midpoints. The aperture radius is 3.5 pixels, the background annulus spans 7–10 pixels, and the local signal uses a linear pixel interpolation and a mean-annulus subtraction. Frames are rescaled to a common zero point using their own header MAGZP and inverse-variance weighted with empirical off-track aperture noise. We search at one-arcsecond spacing inside the Horizons three-sigma ellipse padded by 1.5 arcseconds per semiaxis. The same east/north offset is applied to all exposures that night as the nominal position moves.

Each comparison uses a nearby disjoint ellipse with the same size, search grid and velocity. We drop a control if any point in it has missing image coverage. The reported rank is (1 plus the number of control maxima at least as large as the target maximum) divided by (1 plus the number of controls). This is an empirical local chance-alignment diagnostic. It does not assign a posterior probability to the comet’s identity. The quoted flux scales are five times 1.4826 × MAD of nominal-offset stack apertures on the translated controls. The magnitudes omit color, encircled-energy and coma corrections. An extended low-surface-brightness source or a source outside the trial ellipse can be missed.

The final October 27 exposure has INFOBITS = 33554449, including ZTF’s bit 25 for an overall problematic quadrant. It contributes the lowest source S/N of the three and the smallest inverse-variance weight. The [per-image results](stack-results.json) retain flags, per-frame motion and fixed-sky S/N, ellipse limits, calibration, and provenance. [ZTF quality advisories](https://irsa.ipac.caltech.edu/data/ZTF/docs/ztf_extended_cautionary_notes.pdf) · [ZTF difference flux scale, §6.6](https://irsa.ipac.caltech.edu/data/ZTF/docs/ztf_explanatory_supplement.pdf) · [JPL Horizons definition of sky-plane ellipses](https://ssd.jpl.nasa.gov/horizons/manual.html)
