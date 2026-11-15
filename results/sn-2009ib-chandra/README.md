# A Chandra source at the late-time site of SN 2009ib

Chandra recorded a faint X-ray source **0.342 arcsec from the CHASE position of SN 2009ib**, the Type II-P supernova in NGC 1559. The exposure began on 2016 June 9, **2,507 days after the estimated explosion**. Selecting the photons at the published optical position gives 12 events in a 1.5 arcsec aperture, with 3.11 expected from the local background. The source appears in the Chandra Source Catalog as 2CXO J041739.9−624638. The new result is the measured positional association and late-time photometry. One Chandra epoch leaves open whether the supernova produced those photons. [CHASE discovery coordinates](https://www.das.uchile.cl/proyectoCHASE/listado_sne_2009.pdf) · [Explosion epoch and distance](https://doi.org/10.1093/mnras/stv857) · [Chandra observation 16745](https://cda.harvard.edu/chaser/startViewer.do?menuItem=details&obsid=16745)

## Measurement

The ACIS-S exposure targeted another supernova, SN 1986L, and spent 45.37 ks on this field. We projected the CHASE coordinate (04h17m39.92s, −62°46′38.7″) into the original level-2 event list, selected 0.5–7 keV photons, and placed circular apertures at the optical coordinate. At the primary radius, the 3–8 arcsec background annulus contains 76 events. For each aperture, the local tail probability conditions on the total source-plus-annulus counts and assumes a uniform local background. Nearby sources in the annulus can raise the estimate. [Original photon data](https://cxc.cfa.harvard.edu/cdaftp/byobsid/5/16745/primary/acisf16745N002_evt2.fits.gz)

| Source radius | Background annulus | Counts | Expected background | Local binomial tail |
|---:|---:|---:|---:|---:|
| 1.0″ | 3–8″ | 7 | 1.38 | 0.00074 |
| 1.5″ | 3–8″ | 12 | 3.11 | 0.00017 |
| 2.0″ | 3–8″ | 14 | 5.53 | 0.0029 |
| 1.5″ | 2.5–5″ | 12 | 2.52 | 0.000093 |
| 1.5″ | 4–7″ | 12 | 2.66 | 0.000069 |

The catalog centroid is (64.416461°, −62.777342°). The offset lies inside the 95% positional ellipse, which has semiaxes 0.575″ and 0.433″. CSC reports a 0.5–7 keV aperture flux of 9.82 × 10⁻¹⁶ erg cm⁻² s⁻¹, with a 68% interval of 3.27–16.07 × 10⁻¹⁶. At the supernova-based distance of 19.8 Mpc, the inferred observed luminosity is **4.61 × 10³⁷ erg s⁻¹**. The distance uncertainty of ±3.0 Mpc adds to the photometric uncertainty. The catalog and our aperture measurement use the same photons. [CSC v2.2 field definitions](https://heasarc.gsfc.nasa.gov/W3Browse/chandra/csc.html) · [Supernova distance](https://doi.org/10.1093/mnras/stv857)

## Chance association and prior work

We inspected the CHASE position after screening 24,981 distinct transient names from the 2009–2026 Bright Supernova tables against all 432,713 entries in CSC v2.2. The selection includes published transient class strings starting with Ia, Ib, Ic, II, IV, SLSN or TDE. It yielded 68 matches within 1.4″ and 25 within 0.342″. They include known X-ray supernovae and nuclear sources. SN 2009ib was selected after examining several candidates. To measure the effect of searching this many positions, we moved every optical coordinate by 10″, 20″ or 30″ in 360 equally spaced bearings. At least one match as close as SN 2009ib occurred in **87, 53 and 47 of the 360 shifted catalogs**, respectively (24%, 15% and 13%). These shifts preserve sky coverage approximately and sample different locations in hosts. They show how a chance match can enter a survey-wide search. [Bright Supernova tables](https://www.rochesterastronomy.org/snimages/snredshift.html) · [CSC v2.2](https://heasarc.gsfc.nasa.gov/W3Browse/chandra/csc.html)

To estimate how crowded this particular host is, we placed 100,000 random positions uniformly in area within each of three annuli around the galaxy center. The SN is 36.29″ from the center in projection. The fractions falling within 0.342″ of any CSC source were 83/100,000, 120/100,000 and 76/100,000 for annular half-widths of 5″, 10″ and 20″. These fractions describe the geometry of this field. Spiral arms and young stellar populations can concentrate both supernovae and X-ray binaries. Neither these fractions nor the local photon probabilities account for selection across the survey.

Two prior studies analyzed this observation. Nazri & Annuar plotted the supernova positions and tabulated 24 X-ray sources. Ma et al. tabulated 33 X-ray sources. Neither source list contains the point at the position reported here. The deeper CSC detection is public. We found no explicit published attribution of this X-ray source to SN 2009ib in our search. Late HST imaging resolved the earlier optical progenitor candidate into a blend, underlining the crowding at this site. [Nazri & Annuar, Table 1 and Fig. 1](https://www.ukm.my/jkukm/wp-content/uploads/2020/si3/1/02.pdf) · [Ma et al., Table 1 and Fig. 1](https://arxiv.org/abs/2308.06287) · [Late-time HST imaging](https://arxiv.org/abs/2508.00550)

The closest other CSC source is about 8″ away. A coincident host X-ray binary or an unrelated source therefore remains possible at Chandra resolution, and lower-resolution measurements face blending. A second Chandra exposure with relative astrometry and comparable sensitivity could test for fading. The archived observation establishes the 2016 baseline for that test.

## Reproduction

With Python, NumPy, SciPy, and Astropy installed, run the following commands from this directory.

```sh
python verify.py
python verify.py --download
```

The first command recomputes the aperture counts, luminosity and position controls from the committed radial count table and catalog cutout. The second downloads the original Chandra event file, checks its SHA-256, reproduces the counts in all nine disjoint radial bins, and compares the catalog neighborhood against the live HEASARC TAP service. A local copy of the event file can be checked with `--event PATH`. The background annuli, projected control geometry, RNG seed, event release and optical reference coordinate are specified in the script. The machine-readable numerical output is in `measurement.json`. [HEASARC data access](https://heasarc.gsfc.nasa.gov/xamin/vo/tap/)

The optional `python screen.py CACHE_DIRECTORY` command replays the wider exploratory search and offset controls. It also requires lxml and downloads the public annual optical tables and the 545 MiB compressed HEASARC catalog. Releases can change. `survey-inputs.json` records the SHA-256 hashes of the files analyzed here.
