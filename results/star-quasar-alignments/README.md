# Two stars crossing the sightlines to confirmed quasars

*In memory of Scott.*

**2MASS J18123699+4956048 is on a trajectory that could carry it within 23 milliarcseconds of a quasar in June 2109.** DESI took a spectrum of the blue background object in 2021, when it was still well separated from the red foreground star. The quasar has redshift 1.42456. The earlier 2MASS infrared position agrees with the predicted star position to 0.076 arcseconds in 1998. The prediction depends on extending Gaia's measured motion almost a century. Fresh resolved images in the coming years can test that extrapolation while the objects are still easy to separate.

**A second quasar identification resolves a six-month discrepancy between published predictions for 2028.** Klüter and colleagues already identified the pair involving white dwarf WDJ093026.65+503345.70. Their DR2 catalog put the passage in March and their EDR3 catalog put it in September. A Sloan spectrum identifies the background source as a quasar at redshift 2.58127. Holding its position fixed puts the geocentric minimum on **March 20, 2028**, about 71 milliarcseconds away. We found no detection of gravitational deflection in either system. [Earlier DR2 predictions](https://dc.g-vo.org/amlensing/q2/q/info) · [EDR3 predictions](https://dc.g-vo.org/tableinfo/plc3.data)

## The 2109 prediction

| Quantity | Foreground star | Background quasar |
|---|---:|---:|
| Gaia DR3 ID | 2124576942440058752 | 2124576942438079872 |
| Gaia G magnitude | 16.43 | 18.02 |
| Gaia BP − RP, mag | 2.77 | 0.46 |
| Parallax, mas | 6.087 | −0.005 |
| Eastward motion, mas/yr | +41.826 | +0.093 |
| Northward motion, mas/yr | −0.504 | −0.082 |

On June 17, 2021, DESI targeted Gaia **2124576942438079872** explicitly. The ID is in the coadded spectrum's fiber map. Redrock and QuasarNET give redshifts 1.42456 and 1.42425. The latter identifies C IV, C III] and Mg II, with line-based redshifts 1.4200, 1.4219 and 1.4242. Both reductions use the same observation. The independent Gaia quasar motion is consistent with zero within approximately one standard error in each component. [DESI DR1 observation and archive IDs](evidence-2109.json) · [Gaia astrometry](gaia-2109.csv)

The stars' projected separation in the tangent plane was 3.90″ at epoch J2016.0. The geocentric central prediction is **June 22, 2109**, with a nominal minimum of **22.54 mas**. Propagating the lens's full five-parameter Gaia covariance and the quasar's two-position covariance gave a central 95% range of **13.62–32.53 mas** in 20,000 draws. The formal 95% date interval falls on June 18–25. There are annual local minima because of Earth's orbit. Eight draws had their deepest minimum in 2108. The widest search for the global minimum in each realization spans 2106–2113. Those percentiles cover reported catalog errors within the model. Additional acceleration, small astrometric biases, quasar photocenter variation or unresolved companions would change them. A free propagation using the quasar's small measured proper motion and parallax gives **13.75 mas**. Assigning the star radial velocities of ±100 km/s changes the nominal miss by under 0.08 mas. [Global-minimum calculation and sensitivity tests](long_forecast.py)

Conditional on a stellar mass of 0.2–0.5 solar masses, a dark point-lens calculation gives a **0.43–1.05 mas** shift for the brighter quasar image at the nominal encounter. This is an illustrative mass range. A resolved-image or blending model specific to the observing instrument would be needed to infer a mass. An observing program can already test the relative motion directly. The nominal barycentric separations for 2030 and 2040 are 3.317″ and 2.899″.

### What separate telescopes can establish

- **2MASS** measured a high-quality infrared source on June 27, 1998. Relative to the Gaia prediction, including Earth's displacement, its position differs by −62 mas east and −43 mas north. Its catalog error ellipse has 60 mas axes. This checks the past stellar track.
- **Pan-STARRS** measured both sources on 61 common exposures over 30 nights. A fit to the within-image displacement gives +36.35 mas/yr east and −2.42 north. A night bootstrap gives 95% ranges of +20.05 to +54.41 and −12.18 to +10.62 mas/yr. That fit on its own gives a linear 220 mas passage in 2122. Gaia's motion lies within the Pan-STARRS uncertainty, but the ground data cannot validate the 22.5 mas precision. The separate PS1 mean-catalog motions were +35.47/−4.84 for the star and +2.46/+2.79 for the quasar. [Paired-frame analysis](ps1_motion.py)
- **Palomar photographic plates** extend the test to 1950. On the May 20 red and blue plates the star is recovered within 0.18″ and 0.24″ of the Gaia extrapolation, respectively, with a three-pixel fitting radius and 41–48 local registration stars. The blue plate places the quasar within 0.23″ as well. The 1990s plates blend the sources more severely. In the 1994 infrared plate, the quasar fit shifts almost 3″ toward the star. We retain all five plates, both targets and three fitting windows per target. Simple independent centroid fits on these plates cannot test milliarcsecond geometry. [Centroid offsets](plate-offsets.csv) · [Registration and pixel replay](plates.py)

### Prior forecasts

We checked both Gaia IDs, 2MASS J18123699+4956048, DESI target 39633249072319863, the coordinate designation, and combinations of these with “lens”, “microlensing”, and “quasar” on November 14, 2026. We found no prior report of this specific 2109 encounter. Klüter et al.'s EDR3 survey ended in 2066 and selected faster lenses. Bramich & Nielsen's almanac ended in 2100. Wang et al. selected minima below 20 mas within this century and reported events in 2032 and 2073 at other sky positions. This is a candidate additional prediction. The search cannot establish universal publication priority. [Klüter et al.](https://arxiv.org/abs/2112.12152) · [Bramich & Nielsen](https://arxiv.org/abs/1806.10003) · [Wang et al.](https://doi.org/10.3847/1538-4357/aea148)

## Why March matters for WDJ093026.65+503345.70

SDSS plate 7448, MJD 56739, fiber 398 targeted the quasar in 2014. Its catalog redshift is 2.581268 ± 0.000352, with a clear warning mask. The fitted Lyα, C IV and C III] line fluxes exceed their reported errors by factors of 21.7, 8.5 and 8.6. Gaia resolves two distinct sources. Pan-STARRS also resolves a blue, nearly stationary object and a red object moving southwest. [Sloan line fits](sdss-lines.csv) · [Archive products and hashes](sources.json)

| Geocentric prediction | Closest date | Angular separation |
|---|---|---:|
| Klüter DR2 catalog, event 1018659143306715392-125 | March 15, 2028 | 61 ± 128 mas |
| Klüter EDR3 catalog, same event | September 28, 2028 | 145 ± 64 mas |
| Gaia DR3 lens and stationary spectroscopic quasar, this work | March 20, 2028 | 71.34 mas, central 95% 63.71–78.94 |

The EDR3 analysis estimated background-source motion from a position change between Gaia releases. Its proxy was −2.43 mas/yr east and −17.85 north. The spectrum gives a physical reason to set its Galactic proper motion and parallax to zero. Sampling the lens covariance and quasar position errors over 50,000 draws gives **March 14–27, 2028** at 95% within that assumption. There is no population-based significance attached to those dates. [Published values](literature-values.json) · [Model and uncertainty calculation](forecast.py)

White-dwarf atmosphere fits give 0.44 ± 0.12 solar masses for hydrogen and 0.47 ± 0.10 for helium. The nominal resolved shifts of the major quasar image are **0.67 and 0.71 mas**. At the catalog G-band flux ratio the unresolved centroid shifts are only 0.14 mas and the predicted source brightening is about 0.2 millimagnitudes. Each value uses the central mass of the appropriate fit. The older catalog's assumed 0.65 solar masses gives a larger 0.98 mas resolved shift. Deep, spatially resolved astrometry around March can test the prediction. An ordinary brightness light curve would be a poor test. [Gentile Fusillo et al.](https://doi.org/10.1093/mnras/stab2672)

On Pan-STARRS alone the star's fitted motion is −229.5/−345.7 mas/yr compared with Gaia's −224.8/−353.6. Fourteen foreground and 18 background nights yield a nominal 231 mas minimum in April 2028. The night-bootstrap 16th–84th interval is 72–442 mas. These observations identify the two sources and direction of travel. They leave substantial uncertainty in the observing date and miss distance. The March date in Klüter's first prediction is prior work.

## Search and rejected cases

We screened 1,276,314 unique foreground Gaia IDs against 907,144 type Q/A entries in Milliquas v8 and 1,589,836 unique eligible DESI DR1 cumulative quasar targets. The source catalogs overlap. We queried three stellar samples, in order of precedence for duplicate rows. They were Gaia stars with total motion above 100 mas/yr, parallax signal-to-noise above 5 and RUWE below 1.4, the Gaia EDR3 Catalogue of Nearby Stars, and Gaia stars with parallax above 4 mas, parallax signal-to-noise above 10, RUWE below 1.4 and motion from above 25 through 100 mas/yr. Rows with missing requested values were dropped within each input. DESI selection required a clean coadd fiber status, clean redshift warnings and redshift above 0.1. The screen used a linear J2016 tangent-plane track, the interval 1990–2120 and a 4″ corridor. It returned 128 Milliquas pairs and 206 DESI pairs. These exploratory matches are retained, including many weak or unresolved cases. The two headlines were selected for different reasons: the closest future DESI pair with an initially resolved counterpart and a particularly near-term spectroscopically anchored Milliquas pair. [Milliquas](https://quasars.org/milliquas.htm) · [Gaia nearby-star catalog](https://dc.g-vo.org/tableinfo/gcns.main) · [DESI DR1](https://data.desi.lbl.gov/doc/releases/dr1/)

The closest apparent Milliquas pair separated at epoch 2016 by more than 1.8″ had a 2.7 mas linear encounter in 2020. Its DESI fiber map names Gaia foreground star **1064727654013838336** as its own reference and contains its ≈406 mas/yr eastward proper motion at reference epoch 2020.9651. The reported target coordinate followed the star to the fiber. The spectrum has already been published as the unusual white dwarf WD J0951+6454. We reject a two-body alignment. The next apparent historical case, toward SDSS J222444.15−015901.8, had a Milliquas redshift of 1.527 and an SDSS DR17 pipeline galaxy fit at 0.803 with no pipeline line flux above 3σ. Its background-source interpretation remains unverified. [DESI white dwarf spectroscopy, Fig. 6](https://doi.org/10.1093/mnras/stae2205) · [Screen results, Milliquas](screen-mq.csv) · [Screen results, DESI](screen-desi.csv)

## Replay

Use Python with NumPy, SciPy, pandas and Astropy. From this directory:

```sh
python forecast.py
python long_forecast.py
python ps1_motion.py
python plates.py /tmp/quasar-plate-cache
python screen.py /tmp/quasar-catalog-cache
```

The first three commands run offline. The last two recreate the five photographic plate checks and both catalogs' retained rows. The plate replay verifies pixel hashes to allow regenerated FITS headers. The screen downloads about 1 GB and compares every retained pair with the saved evidence. Archive responses, bands, exposure IDs and SHA-256 hashes identify inputs. The Gaia, SDSS, DESI and Pan-STARRS tables in this directory are selected or derived columns. The simulations propagate stated measurement errors. No simulated observations enter the input catalogs. [Source manifest](sources.json)

Data credit goes to ESA/Gaia and DPAC, DESI Collaboration **Data Release 1**, SDSS, the Pan-STARRS1 Surveys and PS1 Science Consortium, NASA/IPAC and the University of Massachusetts for 2MASS, STScI and the National Geographic Society–Palomar Observatory Sky Survey and Second Palomar Observatory Sky Survey for the DSS images, Eric Flesch for Milliquas, and CDS, GAVO and NOIRLab Astro Data Lab for data access. DESI numeric excerpts were selected and reformatted under CC BY 4.0. [DESI release paper and data](https://data.desi.lbl.gov/doc/releases/dr1/) · [DESI acknowledgments](DESI-acknowledgments.txt) · [SDSS acknowledgments](https://www.sdss.org/collaboration/citing-sdss/) · [Pan-STARRS archive](https://catalogs.mast.stsci.edu/panstarrs/)
