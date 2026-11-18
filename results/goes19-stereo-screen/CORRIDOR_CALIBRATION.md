# Test of a continuous GOES stereo direction corridor

The west-view centroids of two published 2023 bolides lay 4.82 and 17.50 arcseconds from the nominal west-view arcs predicted using the GOES-16 centroids alone. The arcs spanned 584.38 and 555.91 arcseconds. This tests the geometry on two already identified pairs and leaves the acceptance tolerance for a future independent search undetermined.

The observing unit is one 20 ms bin per known bolide. For each event the script selects the highest-energy GOES-16 bin on a Unix UTC-aligned grid, using the earliest bin on ties. It compares the GOES-18 energy-weighted centroid in the same bin. A missing west-view bin is recorded as unavailable. The source tables are pinned to the four exact [published control extracts](../goes-pacific-double-bolide/README.md). Each selected control bin contains ten groups from each instrument. This bin rule differs from the published joint-energy selection and uses no GOES-18 data to select the bin.

The near-side intersections of a GOES-16 sightline with nested ellipsoids whose axes are raised by 40 to 120 km form a straight segment. As seen by GOES-18, the segment traces the minor great-circle arc connecting the two endpoint directions. The code computes the nearest point on that arc, including its endpoints. This covers the intervening nominal shell family without a discrete height grid. Offsets approximate geodetic altitude. The development and receiver directions use the nominal satellite positions and the navigation ellipsoid with 14 km equatorial and 6 km polar axis offsets.

The two controls share a selection history with the published GOES analysis and provide no estimate of accidental matches. Applying a tolerance to a new observation also requires timing, pixel-footprint and pointing uncertainties, visibility from both satellites, and coverage limits. A future candidate outside a frozen corridor with adequate sensitivity would fail that prediction. A direction inside it would need checks for lightning, glint and unrelated associations. The 2026 GOES-18 group data remain unopened.

The run used `4e3849a2dd95b23c2f6b573f42d7262ad74efef0`, Python 3.14.0 and NumPy 2.5.5. The calibration JSON has SHA-256 `227d177e7794153943f214b46496b68641f0c3e2928880608a081e1566695ae9`. On Windows, use a D: directory containing byte-exact copies of the four control CSVs, and supply a new output path under D:.

```powershell
python results/goes19-stereo-screen/calibrate_corridor.py --controls D:/controls --output D:/run/result.json
```
