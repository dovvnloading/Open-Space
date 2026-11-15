# What the earlier catalogs identified

* NASA's DR25 and supplemental KOI tables list KOI-7865.01 as a planet candidate with the stellar-eclipse flag set. The centroid and ephemeris-match flags are zero. DR25 reports a 5.637 ± 0.067 arcsecond displacement relative to the target coordinates, toward the contaminating source. The catalog’s separate offset relative to the out-of-transit light centroid has a different reference position. [NASA KOI columns](https://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html)
* Lissauer et al. list KIC 8110757 / 7865.01 as RPPN in Table 1. The authors’ R designation rejects it because the fitted radius is too large. Its fitted radius is 45.118 Earth radii. This is published evidence against a planet. [Table 1 and status definitions](https://doi.org/10.3847/PSJ/ad0e6e)
* Gaia DR3 already classifies source 2079577157039120128 as an eclipsing binary. A 15-arcsecond VizieR query of the AAVSO International Variable Star Index returns the same Gaia identifier, class E, with period 0.518060 days. The original Kepler Input Catalog resolves the source as KIC 8110754, between KIC 8110757 and KIC 8110758. [Gaia binary table](https://doi.org/10.17876/gaia/dr.3/82) · [VSX entry](https://www.aavso.org/vsx/index.php?view=detail.top&oid=7537728)
* We checked the released tables from Abdul-Masih et al. (2016), Kirk et al. (2016) and Thompson et al. (2018, Table 8) for the three KIC identifiers and both KOI identifiers. None occurred in the downloaded tables. The checks are reproducible with `python prior_catalogs.py`. These checks cover named catalogs, without establishing universal publication priority. [Abdul-Masih et al.](https://doi.org/10.3847/0004-6256/151/4/101) · [Kirk et al.](https://doi.org/10.3847/0004-6256/151/3/68) · [Thompson et al.](https://doi.org/10.3847/1538-4365/aab4f9)

## A previous model explicitly rejected the star that eclipses

On 3 December 2016, the Kepler DR25 astrophysical positional-probability report for KOI-7866.01 fitted an aperture depth of 0.0008961. Table 1 tested KIC 8110754 by name. It assigned the binary a modeled fractional eclipse depth of **2.006** and a relative probability of **7.778 × 10⁻⁹⁶**. It assigned the generic background hypothesis a relative probability of **1**. The two Kepler image localizations and the resolved ZTF observations identify the cataloged star the model had disfavored. DR25 had correctly flagged KOI-7866.01 as a false positive. The change concerns the identity of the eclipsing star. [Report archive](https://exoplanetarchive.ipac.caltech.edu/data/KeplerData/docs/appPdfReports_20170420.tar) · [Method documentation](https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19108-001.pdf)

The report is member `kplr008110758_K07866_01_dr25_app.pdf`. Its payload starts at byte 563366912 and has length 35,215 in the cited tar file. To retrieve the individual report and check the table:

```sh
curl -fsSL --range 563366912-563402126 https://exoplanetarchive.ipac.caltech.edu/data/KeplerData/docs/appPdfReports_20170420.tar -o K07866_app.pdf
sha256sum K07866_app.pdf
pdftotext -layout K07866_app.pdf -
```

The expected PDF SHA-256 is `bfce79c44ed480419d814eefeeb6f2e45062b99c0506470b5f0e126f74de4192`. Adjacent members in this numerically ordered archive go directly from KIC 8110733 to 8110758, then 8110767. This specific archive contains no intervening report for KIC 8110757.
