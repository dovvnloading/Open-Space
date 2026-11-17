# C/2025 V1 is compatible with zero inbound energy once comet outgassing is fitted

At 250 au on the inbound trajectory, C/2025 V1 (Borisov) lies 5.05 formal standard deviations on the unbound side in a gravity-only fit. Fitting radial and transverse water-ice outgassing reduces that separation to 1.02 standard deviations. All four force-law fits use the same 117 ground positions from 27 observatories on 10 UTC nights spanning 9.62 days. The measured arc does not establish an interstellar origin.

| C/2025 V1 force law | RMS per coordinate (arcsec) | Inbound 1/a at 250 au (au⁻¹, formal 1σ) | Nominal v∞ (km/s) |
|---|---:|---:|---:|
| Gravity | 1.048 | −0.01598 ± 0.00316 | 3.77 |
| Free signed radial 1/r² force | 1.044 | −0.01362 ± 0.00426 | 3.48 |
| Marsden water-ice A1, A2 | 1.045 | −0.02682 ± 0.02618 | 4.88 |
| Marsden A1 ≥ 0, A2 free | 1.045 | −0.02682 ± 0.02618 | 4.88 |

Negative 1/a denotes positive orbital energy. The two Marsden runs converge to the same interior solution. The radial coefficient is +63 ± 105 and the transverse coefficient is −59 ± 111 in units of 10⁻⁸ au/day² at 1 au. Neither acceleration is detected. The loss of significance describes how strongly the energy inference depends on an assumed force model.

| V1 positions withheld from the gravity fit | Withheld count | Held-out RMS per coordinate (arcsec) | Inbound 1/a at 250 au (au⁻¹) |
|---|---:|---:|---:|
| First UTC night, Oct 29 | 2 | 0.76 | −0.01834 ± 0.00437 |
| Last UTC night, Nov 8, station G00 | 2 | 2.73 | −0.00923 ± 0.00365 |
| Largest station by number of observations, S89 | 16 | 0.86 | −0.01609 ± 0.00317 |
| All ZTF observations, I41, Oct 29 to Nov 1 | 7 | 1.59 | −0.02197 ± 0.01286 |

The free Marsden fit omitting the last night failed in backward propagation. The remaining complete Marsden holdouts and their residuals are recorded in the results file.

## Independent checks on catalogue examples

| Object | Retained gravity-fit positions | Result from the fitted trajectory |
|---|---:|---|
| 1I/ʻOumuamua | 189 | Inbound v∞ = 26.34 km/s at 250 au. Adding free radial acceleration gives 26.41 km/s. |
| 2I/Borisov | 2,678 | Inbound v∞ = 32.28 km/s at 250 au. |
| 3I/ATLAS | 8,045 | Inbound v∞ = 57.97 km/s at 250 au. |
| C/2024 L5 | 1,281 | Bound before its Saturn encounter. On 2020 Jan 1 at 12.18 au, barycentric 1/a = +0.03997 ± 0.00011 au⁻¹. |
| C/1980 E1 | 180 | Bound on entry at 250 au before its Jupiter encounter, 1/a = (+26.02 ± 1.20) × 10⁻⁶ au⁻¹. |
| C/1997 P2 | 95 | Gravity gives 1/a = (−38.0 ± 23.1) × 10⁻⁶ au⁻¹ at 250 au. An unconstrained inverse-square fit returns a large inward radial acceleration and never reaches the outer boundary. |
| C/2008 J4 | 18 of 25 | Gravity gives −0.01932 ± 0.00367 au⁻¹ at 250 au after seven rejections. With free radial acceleration the fit keeps 21 observations and never reaches 50 au within the integration interval. The inferred energy depends on force law and data selection. |
| C/1877 G2 | 11 | With a 5″ coordinate error floor the eleven micrometer positions give an excess below two formal standard deviations. |

These classifications are consistent with the published planetary-ejection mechanisms for C/1980 E1 and C/2024 L5. C/1997 P2 also has an independent CODE original-orbit solution of (−13.94 ± 13.77) × 10⁻⁶ au⁻¹. C/1997 P2 and C/2008 J4 were both discussed as possible interstellar objects in 2018. [de la Fuente Marcos et al., 2024](https://doi.org/10.1051/0004-6361/202451920) · [CODE C/1997 P2 original solution](https://pad2.astro.amu.edu.pl/CODE/orbit.php?int=1997p2a5&orb=original) · [2018 MNRAS Letters study](https://academic.oup.com/mnrasl/article/476/1/L1/4840245)

## Data and reproduction

The nine compressed `.obs.gz` files in [data](data) contain the exact MPC 80-column records selected by permanent or provisional designation. [The acquisition receipt](data/acquisition.json) gives the source URLs, retrieval times, HTTP status and SHA-256 digests of the source files and the uncompressed object arcs. The HTML observatory table and JPL planetary kernels are external inputs. The short DE440s kernel is identical to the copy verified in the [earlier Borisov–704021 integration](../../borisov-704021/propagate.py). This analysis reuses that script's planetary and solar gravitational constants.

The runner used NumPy 2.5.5, SciPy 1.18.1, pyerfa 2.0.1.5 and jplephem 2.24. Run these commands from this directory.

```sh
python acquire.py
python - <<'PY'
import gzip
from pathlib import Path
data = Path("data")
for p in sorted(data.glob("*.obs.gz")):
    name = p.name.removesuffix(".obs.gz")
    raw = gzip.decompress(p.read_bytes())
    (data / (name + ".obs")).write_bytes(raw)
PY
python refit.py
```

Running `acquire.py` on its own queries public archives anew. Restoring the compressed observations before fitting selects the frozen arcs. Horizons supplies a six-component starting guess at a fixed epoch. The fit optimizes all six components against the measured sky positions. [Seed queries](data/initialization_queries.json) and [states, covariances, residual statistics, holdouts and propagated boundaries](data/refits.json) record the fit inputs and outputs. Position and velocity are in au and au/day in the solar-system-barycentric ICRF. Extra coefficients are in units of 10⁻⁸ au/day².

## Force, observation and error models

- **Sky positions.** Use supported ground-based optical records with identified fixed observatory positions. Remove exact duplicate site/time/position records and exclude satellite, radar and unhandled observation modes. Transform terrestrial observatory coordinates with ERFA's IAU 2006/2000A Earth rotation, apply the geocentric offset from DE440, and iterate light time from reception to emission. The station calculation uses UT1 ≈ UTC, zero polar motion and TT ≈ TDB. Catalogue debiasing, per-exposure centroid calibration and coma-centre corrections remain unmeasured.
- **Weights.** Adopt a 1″ floor per coordinate for modern astrometry and 5″ for pre-1960 data. Limit each station/night's total weight to four positions. Fit with a soft-L1 first pass and two least-squares passes, updating a four-sigma two-coordinate rejection mask between passes. Report retained counts and unweighted RMS per coordinate. Each formal covariance is the inverse normal matrix, inflated when reduced chi-square exceeds one. Sparse arcs can admit non-Gaussian solutions and coherent station or photocentre offsets.
- **Dynamics.** Propagate solar gravity and the eight planetary systems as Newtonian point masses on interpolated DE440 ephemerides with SciPy DOP853. Propagate the state-transition matrix with the trajectory, including the light-time derivative in the astrometric Jacobian. Compute barycentric reciprocal semimajor axis as `2/r − v²/GM_system`. At crossings of 50 and 250 au, include the change in crossing time in the propagated covariance. Record the incoming radial velocity and whether the integration reached each requested boundary. The end-of-integration diagnostic can instead be inside the planet region or outbound. Planet satellites, separate massive asteroids, general relativity and Galactic tides are omitted from this independent fitter.
- **Active comets.** The inverse-square model fits one signed radial acceleration at 1 au. The Marsden variants fit a radial and a transverse coefficient with the symmetric standard water-ice sublimation law. The outward variant constrains the radial coefficient to be nonnegative. These simple laws do not capture arbitrary outbursts, fragmentation or photocentre drift. A formally preferred unconstrained acceleration may be physically implausible. A covariance normal approximation is particularly weak at a parameter bound.

The fits test the designated objects and published arcs. No occurrence rate or survey detection probability follows from this set.
