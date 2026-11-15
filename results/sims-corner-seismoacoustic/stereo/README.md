# Two satellites resolve a lower, brighter pulse over Sims Corner

At 07:48:36.48–36.52 UTC on 2026 August 14, the energy-weighted centers of light from GOES-18 and GOES-19 place the main flash at **47.707° N, 119.389° W and 31.85 km above the reference ellipsoid**. Each satellite contributed 20 optical groups to that 40 ms interval. Their sightlines miss by 0.65 km. The NASA Ames published event selection retains all 40 groups in the peak bin. Three earlier nominal stereo bins yield heights near 39–44 km. Ames retains their West groups and omits all 16 East groups, which NOAA assigned to a separate flash cluster. Their association with this bolide is uncertain. Missing simultaneous bins limit any trajectory or speed inference. Each group aggregates pixels within one instrument frame.

The Meteoroid Environment Office's event page describes fragmentation at about 28 miles (45 km) using ground videos and GLM. CNEOS lists 30 km at peak brightness. NASA Ames already gives a 44.66 km beginning altitude. This reconstruction measures the lower, brightest optical pulse and tests its association with an independently recorded regional acoustic wavefront. The CNEOS height and the original detection of this meteorite fall were already public. [Ames selection and group-by-group comparison](amescheck/README.md) The event page also used GOES data, so its comparison with our result shares an input. The offshore November control below illustrates the size of position errors even with closely intersecting sightlines. [NASA event account](https://fireballs.ndc.nasa.gov/skyfalls/events/20260814-074833) · [CNEOS fireball data](https://cneos.jpl.nasa.gov/fireballs/) · [Meteorite and radar history](https://planetarypat.com/meteorite-falls/sims-corner-2026/)

## Reconstruction and controls

We opened the GOES-East and West GLM Level 2 LCFA granules containing the integer CNEOS second. The public 2025-04-07 through 2026-11-14 CNEOS search required location and altitude. Four rows, including Washington, lie within 25–50° N and 75–125° W. We processed all four. For each satellite, take all groups within 1.8° latitude, 2.4° longitude and three seconds of the listed coordinates and time, bounded by the granule. The coarse bounds allow for parallax. They do not use the listed altitude. Each simultaneous 40 ms bin yields an energy-weighted geographic centroid per satellite. The bin maximizing the product of the two energies defines the main flash.

| Event (UTC) | Groups, West / East | GLM height | CNEOS height | Closest separation of rays |
| --- | ---: | ---: | ---: | ---: |
| 2025-06-26 | 623 / 612 | 30.65 km | 32.6 km | 3.09 km |
| 2025-11-11 | 145 / 916 | 35.65 km | 42.0 km | 0.24 km |
| 2026-03-17 | 0 / 1,313 | No stereo measurement | 45.0 km | — |
| 2026-08-14 | 319 / 123 | 31.85 km | 30.0 km | 0.65 km |

Use the operational lightning ellipsoid (equatorial radius 6,392.137 km and polar radius 6,362.755 km) to recover the reported ray intercepts. The satellite nominal subpoints and heights come from each file. Intersect the two sightlines by least squares and convert the midpoint of their closest points to GRS80 geodetic coordinates. This measures the height of the luminous centroid. The detector pixel size, navigation, spatially extended emission, group clustering and differing views can move that centroid. [Bruning et al. 2019](https://doi.org/10.1029/2019JD030874) · [Operational navigation implementation](https://github.com/deeplycloudy/glmtools/blob/master/glmtools/io/lightning_ellipse.py)

In 48 combinations of bin width (20–120 ms), bin origin, inter-camera time offsets (±20 ms and zero), and group quality filtering, the brightest Washington bin gives **31.29–31.90 km**. Perturbing each centroid independently by up to 5 km in each horizontal coordinate gives **26.93–36.77 km**. At 10 km the range is **22.02–41.70 km**. These are deterministic stress tests. The ray gap and those ranges do not determine a statistical confidence interval. A calibrated pixel-level optical reduction or a clocked multi-station video solution would test the absolute height and the relation between the earlier and later emission.

## Reproduce

From this directory, run the reconstruction with Python 3 and NumPy. It uses the checked-in group rows. To independently re-extract those rows, install netCDF4 and run the second command with a local download cache. It verifies each original granule against its SHA-256 digest.

```sh
python triangulate.py
python extract.py --cache /tmp/sims-corner-glm --output /tmp/groups.csv
cmp groups.csv /tmp/groups.csv
```

[NOAA file URLs and hashes](sources.json) · [All selected groups](groups.csv) · [All stereo bins](stereo-bins.csv) · [Peak and robustness tests](reconstruction.json) · [Seismic and pressure measurements](../README.md)
