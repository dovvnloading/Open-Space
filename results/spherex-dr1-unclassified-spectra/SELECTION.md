# Pre-ranking selection

Use all entries in the AllWISE Source Catalog with 7 ≤ w1mpro < 12 in the 0.50-degree cone centered on (RA, Dec) = (270.00°, +66.56°), ICRS. Intersect with the order-5 nested equatorial HEALPix tile containing the cone center. Retain catalogue rows and explicit edge/quality exclusions. Require 16 pixels of clearance from the tile edge for the primary extraction. Coverage and flags determine usable channels. No source-class, color or residual cut defines the parent sample.

The northern ecliptic pole gives repeated coverage in the mission deep field. One tile bounds the archive transfer and allows local blank-sky and peer checks. The geometry, bright/faint limits and edge clearance were selected before extracting and ranking any spectra. The sample supports a statement about those selected bright sources only.

Compare signed features against known spectral families and other sources; retain broad continuum diagnostics separately. Injections into observed extractions condition sensitivity on the measured field and the recovery criterion. A coadd, a second ranking of the same coadd or a match in channel wavelength alone cannot independently validate a source feature.
