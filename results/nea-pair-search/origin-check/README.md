# Origin constraints on the VE138/VV15 radial acceleration

A hypothesized radiation-pressure acceleration and a high area-to-mass ratio raise the question of space hardware. The four fitted trajectories miss Earth throughout the space age. The closest nominal approach occurred in November 2007 at **1.17 million km from Earth's center** with relative speed **15.1 km/s**. This constrains a direct Earth launch or lunar departure followed by continuous passive propagation under these force models. It leaves open spacecraft maneuvers, separation from an already-heliocentric spacecraft, time-variable forces and errors in the astrometry or force model.

## What was checked

The time interval is 1957 October 4 through 2026 November 15 TDB. An independent heliocentric integration uses the campaign's six-state gravity-only fit and three fits with inverse-square radial/transverse forces. These are alternative nominal histories, **not a covariance or trajectory probability envelope**. In the table, each interval spans the four models. Relative speed is the speed at the encounter.

| Earth encounter | Distance range (au) | Relative speed (km/s) |
|---|---:|---:|
| 1958 Nov 9 | 0.08816–0.08853 | 17.09–17.11 |
| 1988 Nov 14 | 0.02401–0.02404 | 14.61 |
| 2007 Nov 14 | 0.007829–0.007830 | 15.11 |
| 2026 Nov 14 | 0.010556–0.010558 | 14.97–14.98 |

These are all local minima closer than 0.15 au. Under the radial-only fit, the minimum separations across the interval were 0.00932 au from the Moon, 0.20676 au from Venus, 0.15282 au from Mars and 3.31868 au from Jupiter. The dates and distances for every modeled planet and the Moon appear in `encounters.csv`. Published single-apparition ESA results give 0.01056 au and 15.0 km/s for the 2026 encounter. The ESA Rosetta spacecraft really did fly past Earth on November 13, 2007. Its surface altitude of approximately 5,300 km differentiates its trajectory from this object's November 2007 passage. [ESA asteroid record](https://neo.ssa.esa.int/search-for-asteroids?des=2026VV15&sum=1) · [ESA Rosetta operations](https://www.esa.int/Science_Exploration/Space_Science/Rosetta/Rosetta_swing-by_a_success)

The November 15, 2026 snapshots of McDowell's deep-space object and trajectory catalog and Project Pluto's NEOCP artificial-satellite history contain none of the two designations or their temporary IDs, 7V30455 and LVg0046. McDowell lists 458 Sun-centric trajectory phases starting through 2007, including 330 rows giving numerical perihelion, aphelion and inclination for bound orbits. None falls in the explicit broad comparison box q = 0.69–0.89 au, Q = 1.73–2.13 au, i = 14.3–24.3°. This comparison uses rounded catalog elements at different epochs, with incomplete coverage of untracked spacecraft components. The receipt records source URLs, checksums and the exact search. [McDowell's Deep Space Catalog and conventions](https://planet4589.org/space/deepcat/desc.html) · [Project Pluto's NEOCP artificial-satellite log](https://projectpluto.com/sat_id2.htm)

## Physical scale

The reference A1-only fit gives 1.362 × 10⁻⁸ au/day², equivalent to 2.73 × 10⁻⁷ m/s² at 1 au. If this were entirely due to solar radiation pressure, then **Cᵣ A/M ≈ 0.060 m²/kg**, using P = 1361 W/m² divided by the speed of light. Cᵣ includes the optical response. For H = 26 and geometric albedo 0.04–0.40, the usual H–diameter relation gives a diameter of about 42–13 m. At Cᵣ = 1 a sphere of that diameter would need mean density **0.60–1.88 kg/m³**. For comparison, the measured radiation-pressure area-to-mass ratio of small asteroid 2009 BD is 2.97 × 10⁻⁴ m²/kg. Unusual shape, porosity, scattering properties, outgassing and astrometric systematics can invalidate a density inferred under the absorbing-sphere hypothesis. [2009 BD measurement](https://arxiv.org/abs/1106.0564) · [JPL diameter relation](https://cneos.jpl.nasa.gov/tools/ast_size_est.html)

As an order-of-magnitude alternative, a 22.4 m sphere at 1,500 kg/m³ expelling gas at 500 m/s would need about **4.8 g/s** at 1 au if the outflow were perfectly collimated. This follows from momentum balance and shows why the acceleration alone cannot establish an origin. In another NEO a nongravitational acceleration motivated follow-up that found a tail and coma. No corresponding direct activity or spacecraft signature is established here. [Nature Astronomy, 2026](https://doi.org/10.1038/s41550-026-02913-7)

## Reproduce

The script takes the stored gravity fit, the reference fits with acceleration and the JPL DE440s kernel. It verifies the kernel SHA-256, uses the campaign epoch in TDB and takes input positions and velocities in the heliocentric J2000 ecliptic frame. Run with Python 3.11+, NumPy, SciPy and jplephem.

```sh
curl -fsSLo de440s.bsp https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp
python encounters.py de440s.bsp ../orbits/fits/VE138.json ../pressure/fits/VE138.json
```

The point-mass model includes Earth and the Moon separately, planetary-system barycenters, solar first-order relativistic acceleration and the indicated inverse-square components. It omits perturbations from individual asteroids, solar oblateness and planetary relativistic cross-terms. DOP853 integrates the states with relative tolerance 3 × 10⁻¹². Planet positions and velocities are interpolated by cubic Hermite polynomials on a half-day grid. For each model, the code brackets minima on a six-hour grid, then refines in continuous time. Repeating the integration with `--step 0.25` changed the 52 reported values of minimum/encounter distance by at most 0.0000001 au at the stored precision. Screening accuracy and force-model assumptions limit the reported history. [DE440/441 ephemerides](https://ssd.jpl.nasa.gov/doc/de440_de441.html) · [JPL mass parameters](https://ssd.jpl.nasa.gov/astro_par.html)

For the independently maintained hardware catalogs, substitute the local directory containing the three downloaded sources. Changes since the snapshot may change its hashes and counts.

```sh
curl -fsSLO https://planet4589.org/space/deepcat/tab1.tsv
curl -fsSLO https://planet4589.org/space/deepcat/tab2.tsv
curl -fsSLO https://projectpluto.com/sat_id2.htm
python catalog_check.py .
```

## 1988 photographic-plate search


At the 1988 encounter the nominal H–G light curve (H = 26, G = 0.15) reaches approximately V = 20.2. We checked all six November 1988 plates in STScI’s DSS master list against the predicted track at their recorded epochs. The closest plate center was **42.21°** away. The Wide-Field Plate Database returned 696 records in the interval 1988.845–1988.884. Of the 689 records with coordinates, the nearest was **Asiago ASI040 017772** on November 12 with a 600-second exposure. The nominal track at exposure midpoint was **5.39°** from its center, outside the 4.85° half-diagonal implied by its 12 cm square plate and 100 cm focal length. The next closest records are more than 8° away. These catalog checks yielded no third-apparition observation. Seven WFPDB records lack positions and the catalogs are incomplete. The predicted magnitude assumes a standard asteroid phase law; the search does not establish an image-based upper limit. [DSS master-list documentation](https://outerspace.stsci.edu/pages/viewpage.action?pageId=261164277) · [WFPDB record definitions](https://dc.g-vo.org/tableinfo/wfpdb.main) · [Asiago telescope and plates](https://openaccess.inaf.it/bitstream/20.500.12386/27630/1/aa32131-17.pdf)

```sh
python plate_search.py de440s.bsp /path/to/VE138.json ../pressure/fits/VE138.json
```

The script retrieves both catalogs, propagates the A1-only orbit at the individual exposure times, and prints `plate-receipt.json` to standard output. The receipt includes time-filter counts, checksums, the nearest three exposures per catalog, center and target coordinates, and the assumed magnitude. Predictions are geocentric and include first-order light time. Potential topocentric offsets are of order 0.1° at closest approach. Thus even the maximum solar-system parallax and the motion during the nearest exposure do not bridge its center-to-corner shortfall.
