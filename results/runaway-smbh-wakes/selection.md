# Observational selection

The coordinate-free blind search samples a central square of the released VIS mosaic at every available Q1 MER tile ID. The four tiles used by the earlier Q1 lens search receive 102.4 arcsec squares. The remaining tile centres receive 51.2 arcsec squares. The source of the list is the public IRSA MER prefix. A tile can be partly empty. Deduplicate overlapping areas where appropriate and account for image margins before reporting searchable area. The wider squares reuse tile identities, not the positions or triggers selected by the ring filter.

Separately query IRSA's full Q1 MER catalogue for VIS-detected entries with reported ellipticity >0.95 and <1 and semimajor_axis >=7.5 native pixels. Ellipticity is 1 minus the axis ratio; the catalogued semimajor axis is not an endpoint-to-endpoint length. Retain detection flags and measure apparent length and width in calibrated images. A catalogue object may omit or split a diffuse wake. No catalogue-area completeness follows from this filter.

RBH-1 is a positive control for linear emission. A known target supplied to a cutout service can test a blind detector within the cutout. Measuring flux along its known position angle is a separate, targeted measurement.
