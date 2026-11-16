# Independent counterpart searches

Run these commands from the mission directory after the catalogue and DR14 replays. Python's standard library suffices for the cone requests. Astropy, SciPy, NumPy and Matplotlib support the image, astrometric and figure replays. The cache holds the full response rows, query URLs and failures. Only summaries and scientifically material cutout and observation receipts are kept in Git.

~~~sh
for queue in measured boundary rescue bright difference
do
  python analysis/counterparts.py results/$queue.csv cache
  python overlap/match_prior.py results/$queue.csv overlap/prior-positions.csv overlap/$queue-prior.csv
done
python analysis/ledger.py results overlap cache
python analysis/boresight.py 5XMM_DR15cat_source_v1.0.fits 3030258180100013 cache results/j085618-boresight.json
python analysis/optical_depth.py dr9/north 1339p380 134.077254 38.119831 3.204 1.6005265e-14 cache results/j085618-depth.json
python analysis/ps1_depth.py cache > results/j180019-ps1-depth.json
python analysis/cfht_photometry.py cache results/cfht-photometry.json
python analysis/plot_field.py cache/cfht-R.MP9601.fits products/P0302581801PNS003SRSPEC000E.FTZ results/j085618_fit.json results/j085618-field.png
~~~

- Gaia DR3 was queried to 90 arcsec. The separation is the closest approach over the source's XMM epoch range, with a straight-line proper-motion correction from Gaia epoch 2016.0, treating unavailable proper motion as zero. Missing motions, annual parallax and stars fast enough to enter from beyond 90 arcsec can affect individual fields. The Gaia field registration for J085618 used the actual 2005 epoch. Eight reciprocal Gaia matches were found among 37 usable DR15 sources from the same observation. The eight control translations of 30 arcsec produced counts of 1, 0, 0, 0, 0, 0, 0, 0. The median offset is a registration check, with the complete nominal uncertainty retained.
- AllWISE was queried to 30 arcsec. Its listed nearest object and W1 magnitude are contextual. Blending, contamination flags and the broader infrared PSF matter. A positional neighbor alone supplies no association probability. The code does not infer an infrared limiting magnitude from a missing record.
- Legacy Surveys DR10 Tractor was queried to 30 arcsec using NOIRLab's public endpoint. The northern footprint uses DR9 imaging. The table records the nearest Tractor entry, even where r flux is unavailable, and converts positive r-band nanomaggies to AB. Maskbits, forced flux errors, number of exposures, and depth at the source position matter for a physical limit. J085618 received independent map and deep CFHT checks. Empty Legacy cones in e.g. the Galactic plane can be outside the footprint.
- Pan-STARRS1 DR2 uses a 30 arcsec stack cone north of −31 degrees. The nearest primary detection is retained for audit, including low-significance stack rows. Separate columns give the nearest primary row with positive r photometry, an error below 0.22 mag and at least one detection. Valid responses were received for 218 sources; one additional query at declination −30.33° returned an empty body without a table header and remains unavailable. The other 119 positions fell south of the query boundary. None of these categories yields an optical limit.
- The release-wide SIMBAD and eRASS1 queries use larger review cones than the counterpart localization. The ledger's triage gives priority to in-radius prior-work hits and then lists surveys with a possible neighbor inside the nominal XMM 99% circle. The infrared and optical tests use the XMM radius alone and do not account for every comparison catalogue's centroid errors. They flag fields to review, without asserting astrophysical identity. Source-specific decisions with original evidence supersede automatic labels. There was no high-ratio discovery threshold applied to the entire sky, since local depth measurements are available here only for the individually followed fields.

[Counterpart summary](results/counterparts.csv) · [Review ledger](results/review-ledger.csv) · [Gaia TAP](https://gea.esac.esa.int/tap-server/tap) · [AllWISE](https://irsa.ipac.caltech.edu/Missions/wise.html) · [Legacy Surveys photometry and depth](https://www.legacysurvey.org/dr10/files/) · [MAST Pan-STARRS archive](https://catalogs.mast.stsci.edu/panstarrs/)
