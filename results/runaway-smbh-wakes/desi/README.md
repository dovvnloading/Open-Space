# DESI DR1 coverage and redshifts

Across the 68 host-associated Euclid segments reviewed in calibrated VIS, only one distinct DESI DR1 target falls in the search cones. Two overlapping cones contain that galaxy. Both record the same secure pipeline solution (z = 0.531167, ZWARN = 0, Δχ² = 453.9). No segment has distinct spectra along its length and at its proposed host that could establish a shared redshift. The HST western line lies outside published DESI DR1 HEALPix coverage.

The [69 cones](direct_cones.csv) use a radius of half the measured segment length plus 5 arcsec for Euclid and 10.8 arcsec for HST. Of these, 39 have zero intersecting HEALPix products and 30 have at least one. Every one of the 63 distinct intersecting survey/program/pixel redrock files was searched successfully. The code queries the NSIDE 64 nested pixels intersecting the complete cone, checks them against every survey and program in the public iron tilepix manifest, then measures great-circle separations using the fiber coordinates. Spectra in the same cone can belong to other sources. [Two cone hits](direct_matches.csv) · [Parent FITS digests](direct_inputs.json) · [DESI DR1 archive](https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/)

From this directory run `python direct.py --cache /path/to/desi-cache` with Astropy, astropy-healpix and NumPy installed.
