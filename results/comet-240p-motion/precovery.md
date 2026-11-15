## A check for an earlier image of 240P-B

Two Blanco 4-m DECam exposures from UT 25 September 2023 cover JPL’s backward extrapolations for A and B. The images were taken only 2.5 minutes apart. Both positions fall on CCD N19. The later orbit solutions used in this calculation are K252/47 and 19.

| Calibrated exposure | Time, JD UTC (midpoint) | Band | Exposure | B predicted RA, Dec (ICRF degrees) | Distance to Gaia DR3 4041916898212244608 |
| --- | ---: | --- | ---: | --- | ---: |
| c4d_230925_000357_ooi_r_v1 | 2460212.5043713 | r | 140 s | 266.834282835, −33.672899057 | 0.74″ |
| c4d_230925_000646_ooi_z_v1 | 2460212.5058640 | z | 100 s | 266.834416239, −33.672890014 | 0.34″ |

The tabulated position of this G=21.14 Gaia star falls within one seeing disk of B’s extrapolated track. The DECam headers give FWHM near 3.7–3.8 pixels, about one arcsecond at 0.263 arcsec per pixel. A second G=20.78 Gaia star lies 1.45–1.64 arcsec from the predicted positions. Both stars lack Gaia proper motions. The two bright image peaks closest to A also coincide with catalogued field stars. The field is heavily crowded and the predicted motion between bands is only about 0.40″. We could not identify either comet component securely from these two exposures. An unmodeled change in outgassing would also move B’s extrapolated position. This check supplies neither an earlier detection nor a flux upper limit. [DECam r exposure](https://astroarchive.noirlab.edu/api/retrieve/998c2f06c3a93bab58bb2b245548e4b4) · [DECam z exposure](https://astroarchive.noirlab.edu/api/retrieve/b5e672097f62fe3d16f5a5c9b364704e) · [Gaia archive](https://gea.esac.esa.int/archive/) · [Horizons API parameters](https://ssd-api.jpl.nasa.gov/doc/horizons.html)

The positions can be regenerated with a Horizons observer-table query for targets 90001213 and 90001214, center W84, quantities 1 and 9, extra precision, decimal degrees and the two UT Julian dates in the table. The Gaia DR3 positions can be retrieved by source ID. Station W84 is the Blanco 4-m on Cerro Tololo.
