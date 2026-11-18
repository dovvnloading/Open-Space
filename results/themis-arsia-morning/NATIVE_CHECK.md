# Native radiance diagnostic

This added diagnostic was frozen after the EQR test and before accessing VIS-RDR science arrays. It tests whether a scene-level brightness increase is present in the radiometrically calibrated detector arrays before the EQR projection and framelet filtering.

Use the November and January native RDR labels to select 0.42 and 0.75 µm. Apply each label's signed 16-bit scale and value offset after masking stored special values. Compare median positive radiances over the central half of lines and samples, and repeat over the entire available plane. These are scene medians in detector coordinates. The footprints differ in length and cannot provide an exact same-ground test in this diagnostic.

Also compute an illustrative illumination ratio from the two center incidence angles and heliocentric distances. Assume a flat Lambertian surface plus inverse-square illumination. Report raw and geometrically scaled ratios. This calculation omits local slope, phase dependence, aerosols and differences in footprints. Uncertainty in those terms leaves the origin of a scene-wide change unassigned. [THEMIS calibration](https://doi.org/10.1029/2005JE002568)
