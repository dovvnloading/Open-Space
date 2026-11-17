# L 788−34 field source — colours and background surface densities

The extended source lies about 3.2 arcsec from L 788−34 (GJ 4274). The large-aperture F200W flux exceeds F444W. Its colour is incompatible with a cool point-source companion drawn from the below model set at the distance of the star. Comparable blue galaxies occur in a measured JWST control field. These are independent checks on the source morphology, with no assumed photometric redshift for this source.

**Conditional companion check.** At the host parallax of 138.228 mas (7.23 pc), we selected all unstarred Sonora Bobcat grid rows with F444W 17–45 µJy. This deliberately broad bracket spans alternative source apertures, rather than representing a statistical flux interval. The selected 50 solar-C/O, cloud-free equilibrium models span 200–275 K and [M/H] = −0.5, 0.0, +0.5. Their model radii give F200W ≤0.00134 µJy; F200W−F444W in AB magnitudes ranges from +11.30 to +17.99. The broad-aperture source colour is approximately −1 in the same AB system. All calculations use the tables of log flux densities (mJy at 10 pc), parallax dilution, and −2.5 log10(F200W/F444W), avoiding Vega/AB zero-point mixing. The comparison does not constrain all atmospheric chemistries, or universally rule out a much warmer, unrelated object. The previously catalogued 1.63-day and 69.57-day RV planets have different separations. [Marley et al. model files](https://doi.org/10.5281/zenodo.5063476), [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/overview/GJ%204274)

**Background galaxies.** We counted sources in the JADES DR2 GOODS-S catalogue using F444W-resolution PSF-matched Kron photometry. A 30×30 arcsec cell was retained when each of its 15×15 arcsec quadrants contained at least five catalogue sources observed in both bands. We counted all sources with such coverage inside the same cells used for the area estimate. The effective 61.5 arcmin² denominator is a catalogue-based coverage approximation. The control is a single deep field and cosmic variance is outside the Poisson errors. We used the catalogue's star/galaxy flag equal to 0 for galaxies.

| Lower F444W flux (µJy) | All sources | Catalogue galaxies | Galaxies with −1.5 ≤ F200W−F444W (AB) ≤ −0.5 | Blue galaxy density (arcmin⁻²) |
|---:|---:|---:|---:|---:|
| 20 | 241 | 214 | 45 | 0.732 |
| 25 | 186 | 165 | 38 | 0.618 |
| 35 | 114 | 99 | 29 | 0.472 |
| 40 | 91 | 79 | 26 | 0.423 |

The ≥20 µJy blue-galaxy sample scales to 0.0066 expected objects within 3.21 arcsec of a single random position. This conditional expectation depends on the post hoc radius, flux and colour cuts. The full search covers a larger area and multiple targets. Spectral or astrometric classification of the individual source takes precedence over number counts. [JADES DR2](https://archive.stsci.edu/hlsp/jades), [CDS catalogue and flags](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/V/159)

**Local optical field.** The 5 arcmin Gaia DR3 cone centred at (335.78238°, −17.61240°), excluding the host, contained 28 entries with G<18 and 67 with G<20 (0.36 and 0.85 arcmin⁻²). These are counts of Gaia detections as a foreground-star comparison, with different completeness and bandpass from the infrared galaxy control. The exact ADQL, return values, query date and response hashes are in `gaia_field_counts.json`. [Gaia TAP](https://gea.esac.esa.int/tap-server/tap)

Run `python background.py CACHE OUTPUT` with NumPy and SciPy. The code verifies row identifiers between the CDS tables and records SHA-256 and byte counts for its three source files in `inputs.json`. The full conditional model subset and count thresholds are in the accompanying CSV files. No catalogue rows are interpolated or extrapolated.
