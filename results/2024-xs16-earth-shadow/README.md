# A Sicilian telescope measured 2024 XS16 fading as it entered Earth's shadow

On 11 December 2024 the Galhassin Robotic Telescope at Isnello (MPC L34) caught the ingress of a predicted Earth eclipse of the few-metre near-Earth asteroid **2024 XS16**. The asteroid faded from **G = 16.4 ± 0.17 to 17.0 ± 0.21** in 5 min 18 s. Two separately fitted orbits predict 3.3% and 55.9% average coverage of the solar disk during the respective exposures. The same telescope recovered the asteroid at **G = 15.0 ± 0.16** after egress. The orbital geometry predicts about **43 min 38 s of totality**. Only partial ingress has supporting MPC photometry. [All 42 MPC records](mpc.csv) · [JPL model and measurements](measurement.json) · [ESA orbit and independent calculation](esa/README.md)

| L34 midpoint (UTC) | Reported exposure | Measured G | Mean solar disk covered | Expected V without shadow |
| --- | ---: | ---: | ---: | ---: |
| 17:48:59.7 | 660 s | 16.8 ± 0.09 | 0% | 17.12 |
| 18:07:02.9 | 330 s | 16.5 ± 0.12 | 0% | 16.92 |
| 18:27:41.1 | 219 s | 16.4 ± 0.17 | 3.3% | 16.64 |
| 18:32:58.9 | 219 s | 17.0 ± 0.21 | 55.9% | 16.57 |
| 20:26:58.2 | 219 s | 15.0 ± 0.16 | 0% | 15.70 |

The orbit predicts an increase in brightness of 0.076 mag between the two ingress records if the asteroid stays sunlit. Including the shadow predicts a **0.775 mag fade**. L34 measured **0.600 ± 0.270 mag**, using the two quoted errors in quadrature. The observed change differs from the smooth sunlit prediction by 2.50 such errors and from the eclipse prediction by 0.65. These are residuals conditional on the model and quoted errors. Rotation, trailing, magnitude calibration, variable exposure weights and sparse sampling prevent the photometric pair alone from settling the cause. The modeled hourly trend uses V as a proxy for the differential G-band change over five minutes. The absolute magnitudes have a filter and calibration offset.

## A second telescope before ingress

Schiaparelli Observatory (MPC 204) recorded six one-second exposures of CCGRUC2. All six ESA display images show a source at the moving position predicted from the station-specific orbit and the image astrometry. The mean signal on that track exceeded all 296 translated comparison tracks. The last shutter closed at **18:24:28.7 UTC**, 176 s before modeled first contact. These images support the pre-ingress identification and position. The public display stretch limits calibrated photometry. Only L34 supplies a measured fade during ingress. [Six-frame replay, FITS-header astrometry and control tracks](schiaparelli/README.md)

## Independent orbit test

JPL Horizons geocentric geometric vectors and an integration from ESA Aegis elements both place first contact near **18:27 UTC**, second contact near **18:37**, third contact near **19:21**, and last contact near **19:28**. They agree to within one second for this nominal spherical-Earth calculation. The two services can use the same optical astrometry. Their agreement tests the orbit calculation. It is separate from the telescope flux comparison. Aegis covariance draws keep the second exposure more than half in shadow. The JPL estimate ranges from 55.1% to 56.2% across polar and equatorial Earth radii. Shifting the nominal exposure window by ±30 s instead gives 49.5% to 62.2%. [ESA source, propagated covariance and sky-position residuals](esa/README.md) · [L34 observer ephemeris](l34-unshadowed.csv)

The geometric model integrates unocculted solar area over 219 s with uniform time weighting. It assumes a uniformly bright solar disk and an opaque Earth of mean radius 6,371 km. The nominal contact times use geometric Earth and Sun positions with asteroid-to-observer light time. Apparent solar direction, limb darkening, atmospheric extinction/refraction and observatory-specific exposure timing can affect the detailed light curve and contact seconds. The eclipse model makes no use of the five magnitudes to choose its time or depth.

## Search population and observation history

The JPL close-approach query returned 933 nominal Earth approaches within 0.002 au in 2018–2050. We finished the 174 passes closer than **0.0005 au** (173 had occurred by the search date). Each received a ±1 day geometrical scan, with 5-minute Horizons state vectors, cubic Hermite interpolation and a 10-second eclipse grid. Fourteen past encounters passed the nominal 5% solar coverage screen. The MPC exports for those targets supply 398 optical records with magnitudes. One had its reported midpoint inside greater than 5% obscuration, the 18:32:58.9 L34 measurement. Sparse detections and reduced discovery efficiency during eclipse strongly affect this count. The search did not examine every past asteroid eclipse. [Selection, denominators and response hashes for all 174](encounters.csv) · [Archive queries](sources.json)

Earth eclipses of near-Earth asteroids have been predicted and measured before, including 2012 XE54 and 2016 VA. The discovery designation for 2024 XS16 was CCGRUC2. The MPC announcement and Gianluca Masi's image document the close approach. Web and publication-index searches on 2026-11-15 used both identifiers with “eclipse”, “Earth shadow”, “penumbra” and “umbra”. The records located describe the approach and imaging. We found no previous event-specific ingress measurement. A private observing log or unindexed notice could establish prior reporting. [Birtwhistle's 2018 eclipse photometry](https://mpbulletin.org/issues/MPB_45-3.pdf) · [Original MPEC 2024-X303](https://minorplanetcenter.net/mpec/K24/K24XU3.html) · [Virtual Telescope observing account](https://www.virtualtelescope.eu/2024/12/12/near-earth-asteroid-2024-xs16-extremely-close-encounter-an-image-11-dec-2024/)

## Reproduce

Python 3, NumPy 2 and SciPy reproduce the photometric comparison from the compact checked-in state, magnitude and observation tables. The live screen reads the frozen selection, caches archive responses in a directory supplied by the caller and prints one row per processed approach. Omit `--only` to rerun the 174-pass screen. Later orbit updates can change the output. Horizons also timestamps response headers.

```sh
python measure.py
python esa/shadow.py
python replay_screen.py /tmp/xs16-screen --only '2024 XS16'
```

The most discriminating further observation would be the individual calibrated L34 frames and their exposure timestamps around 18:27–18:37 UTC. They could test the predicted ingress curvature against a rotating light curve and show how much weight the combined measurements gave each frame.
