# A 31,000 km encounter between 407P/PANSTARRS–Fuls and (75547) 1999 XV257

On 22 February 2026, numbered main-belt asteroid (75547) 1999 XV257 passed within approximately 31,000 km of comet 407P/PANSTARRS–Fuls at 8.06 km/s, about five months before comet perihelion. The Sun–comet distance was 2.427 au. We identified this three-dimensional encounter by screening catalogues and independently evaluated the nominal geometry with two ephemeris services. An image search found the asteroid moving on the predicted track on the nights bracketing the encounter.

| Numerical ephemeris | Closest-approach time | Separation | Relative speed |
| --- | --- | ---: | ---: |
| JPL Horizons, geometric heliocentric vectors | 22 Feb 2026 19:55:35 TDB | 31,373 km | 8.05951 km/s |
| IMCCE Miriade, INPOP19A | 22 Feb 2026 19:54:07 TT | 31,688 km | 8.05947 km/s |

JPL's minimum corresponds to **19:54 UTC**, to the nearest minute. Independent service results differ by 315 km and 88 seconds. They can share underlying astrometry. The JPL separation grows to 34,325 km on either side at ±28.8 minutes and 176,890 km at ±6 hours. The encounters here describe the centers of two known objects. The heliocentric origin, common epochs and geometric coordinates remove the line-of-sight ambiguity of a projected alignment.

## How the pair was found

The exploratory screen included 894,699 numbered JPL asteroids (element epoch MJD >59000, semimajor axis <10 au, eccentricity <0.98) and 758 comet entries in the final search. We propagated osculating elements across an eight-day grid from 2025 through 2028, using heliocentric Horizons samples for a subset of comets. For each grid point we searched a spatial tree within 0.19 au and estimated closest approach within ±4.05 days using the relative velocity. The MPC osculating solution for 407P and the numbered asteroid elements yielded **30,675 km** before refinement. This approximate search can both miss encounters and generate false positives, especially far from an elements epoch. These counts describe the search undertaken and do not establish completeness or an encounter ranking.

For refinement we used each body's heliocentric position and velocity to solve for the zero of relative position dotted with relative velocity, then requested new Horizons states at the refined epoch. We independently repeated the local minimization with Miriade mean J2000 heliocentric rectangular states at the JPL minimum. The bundled evidence includes both sets of state vectors, the JPL time brackets, the exact service URLs and checksums of the screening catalogues. [JPL numbered-body elements](https://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR.gz) · [MPC comet elements](https://minorplanetcenter.net/iau/Ephemerides/Comets/Soft00Cmt.txt) · [Horizons documentation](https://ssd.jpl.nasa.gov/horizons/manual.html) · [Miriade documentation](https://ssp.imcce.fr/webservices/miriade/api/ephemcc/)

## Archive images and limits on the interpretation

We measured all 12 public ZTF science exposures returned for 22–25 February in g, r and i. The predicted asteroid position was within the central usable portion of six cutouts on the first two nights. Four exposures have a matched-filter peak within 0.63 pixel of the predicted moving position, with point-source signal-to-noise ratios of **5.84, 10.51, 10.57 and 13.10**. Two other positive measurements have their peak on the edge of a two-pixel search box and were excluded from the four clean positional matches. The final six images only constrain the comet position. Positions use exposure midpoints and the Palomar topocenter. A median-filter sky estimate, seeing-matched Gaussian and robust local noise estimate are described by the code. These image statistics are diagnostics and do not provide calibrated comet photometry. [IRSA ZTF image data and cutout interface](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html)

At the **nominal** comet position the 12 point-source measurements range from −1.59 to +0.61 times the empirical noise. Diffuse emission and source locations displaced within the comet's orbit uncertainty require different tests. JPL solution 19 uses 734 observations from 2013-05-05 through 2020-06-12. At closest approach the reported formal three-sigma RA, Dec and range errors for the comet are 8.55 arcsec, 0.416 arcsec and 2,631 km. The larger angular error corresponds to approximately 10,350 km taking the full angular value at 1.669 au. The 315 km service difference does not measure the total orbit uncertainty. We have no direct measurement of the comet nucleus or of coma material at the asteroid's location. A dust crossing, impact or resulting physical change remains unestablished.

Paired catalogue-number, provisional-designation and comet-name searches found no published report of this encounter during preparation on 14 November 2026. An independent prior calculation or an unindexed report could exist. The literature already contains asteroid–comet trajectory screenings, including the 3I/ATLAS survey. The proposed new result is this specific passage. [3I/ATLAS encounter survey](https://doi.org/10.3847/1538-3881/ae2ea6)

## Reproduction

Run from this directory. The geometry calculation uses the Python standard library and recorded inputs. The optional image calculation needs NumPy, SciPy and Astropy, and downloads the cutouts listed in the manifest.

```sh
python verify.py
python verify.py --frames /path/to/cutout-cache
```

- `ephemerides.json` contains the independent numerical inputs, source URLs, astrometric uncertainties and exposure-ephemeris queries.
- `ztf_frames.csv` lists exposure midpoints, passbands, both objects' predicted positions and the public FITS cutout URLs.
- `image_measurements.csv` contains the full forced-position results, including image-edge exclusions and crowded-field flags.
