# Three independent Uranus visits in January 2024

**None of the 267,000 circular trial orbits per filter in the three inner search intervals survived the glare, known-moon and two-visit tests.** A deliberately unmasked control recovers the source at Perdita's longitude in all three visits in both filters. The brightest anonymous-looking outer-gap feature lies in Miranda's diffraction pattern. The two filters are F182M and F210M, each observed in three physical visits. [Raw and vetted grids](grids/) · [Known-body control and native astrometry](perdita/)

## Search footprint and controls

JWST program 2768 provides fifteen 107.37 s exposures per filter, five in each of visits 02768001001, 02768002001 and 02768003001. Ninety-six SCI/ERR/DQ detector tiles intersect the central cutout. Most dithers place the planet itself in a detector gap. Tiles from a single physical exposure share one registered sky mosaic. Empty pixels carry zero weight. Each visit uses its own fitted planet center offset derived from fully enclosed disk images. The inputs and hashed spacecraft ephemerides record exactly which footprint, times and known-moon paths enter each search. [F182M wide inputs](f182_wide_inputs.json) · [F210M wide inputs](f210_wide_inputs.json) · [Source observables](uranus_observables.csv) · [JPL receipts](uranus_sources.csv)

| Region and grid | Trials per band | F182M supported | F210M supported | F182M / F210M initial positive | Surviving |
|---|---:|---:|---:|---:|---:|
| 52,000–65,000 km, 40 km × 1° | 117,360 | 52,743 | 52,682 | 75 / 542 | 0 / 0 after local ring-light control |
| 72,000–79,000 km, 40 km × 0.75° | 84,480 | 79,512 | 79,488 | 0 / 0 | 0 / 0 |
| 108,000–126,000 km, 100 km × 1° | 65,160 | 59,303 | 59,225 | 15 / 1 | 0 / 0 with Miranda halo exclusion |

Each exposure is sampled 17 times. The reference emission JD TDB is 2460347.42. A trial needs a score of five and a worst leave-one-exposure-out score of two in each of at least two visits; the pooled score must reach seven. These counts concern the specified circular, equatorial grid and matched-filter score. The empty regions, eccentric and inclined orbits and sensitivity as a function of source brightness require their own coverage tests. The project injection measurements from other sequences do not set a magnitude limit here.

For the 617 initial inner-region passes, the same-radius local background is the median of ten phase-displaced, separately measured paths in each exposure (±6°, 8°, 10°, 12° and 15°). No nominal pass remains after subtracting that background: the highest F182M/F210M ranks are 0.205/0.177 of threshold. The wider survey uses the measured bright-moon masks. For Miranda profiles, the reduction measures 72 fixed position angles at every radius from 7 to 100 output pixels. A repeatable five-score feature in two visits persists out to 59 pixels. A three-pixel guard gives an updated exclusion radius of 62 pixels (1.92″). The strongest raw object-like F182M trace stays only 21–23 pixels from Miranda in its two strongest visits. The third visit has a score of 1.21 at a 58-pixel separation. The positive and inverted searches return zero after excluding Miranda's contaminated footprint; the two filters retain 56,023 and 56,089 supported traces respectively. [Miranda radial profile](miranda/radial_profile.csv) · [Calibrated exclusion](miranda/receipt.json) · [Inner F182M background audit](grids/F182M_inner_flanks.csv) · [Inner F210M background audit](grids/F210M_inner_flanks.csv)

## Perdita and the 2025 orbital aliases

Unmasking only Perdita gives a January longitude sweep across the full circle, with the remaining 28 known moons still masked. The strongest F182M path is 1.00° behind nominal with visit scores 7.37, 14.09 and 12.17. In F210M the peak is 0.75° behind nominal, with scores 6.30, 18.46 and 11.88. Eighteen calibrated native exposure tiles supply independent source centroids with a trailed PSF plus local planar background. Formal fit S/N is at least 9.0 in F210M and 11.1 in F182M; the nearest other projected moon is over 15 registered pixels away. The median projected longitude residual relative to predicted Perdita is −0.830° (MAD 0.051°). Ten fits to the comparison moon Belinda give a median of −0.779°. This January sub-degree offset is comparable to the bright-moon registration/ephemeris control. Position-fit errors exclude shared registration systematics. [Jan F182M astrometry](perdita/F182M_native.csv) · [Jan F210M astrometry](perdita/F210M_native.csv) · [Comparison moons](perdita/frame_longitudes.json)

The 2025 result supplies 31 circular aliases trained entirely on February and October 2025 astrometry. Propagating those solutions without refitting into January 2024, the 387-turn solution at 76,418.119 km comes closest to the January source, with an RMS residual of 14.24 pixels (0.44″) and pooled forced scores of 0.37 and −3.45. A naive extrapolation of the independent October-only short-arc circular solution also fails to recover the January source. Adding first-order eccentricity and inclination components and J2/J4 epicyclic precession to the 2025-only fit changes the best training alias to 388 turns: February and October RMS residuals of 0.22 and 0.24 pixels, but a January miss of 322 pixels. The closest eccentric/inclined alternative (387 turns) still misses by 14.05 pixels. These sensitivity fits exclude inter-moon perturbations. The 387-turn radius matches JPL's tabulated 76,418 km mean semimajor axis for Perdita. Perdita's published e≈0.005 and i≈1.6°, interactions with neighboring moons, discrete revolution aliases, and registration uncertainty make the two-epoch solution inadequate for precise multiyear extrapolation. The 2024 data alone cannot identify a separate body. [Frozen 2025 aliases](perdita/locked_2025_aliases.csv) · [Eccentric/inclined sensitivity](perdita/epicyclic_aliases.csv) · [JPL mean elements](https://ssd.jpl.nasa.gov/sats/elem/)

## Reproduction

Inside the parent result directory obtain the program's CAL cutouts with the published product manifest and `jwst_cutouts.py` using `--half 330`. Use the F182M and F210M subsets and the URA184 kernel excerpt for the provisional moon. The three sampled ephemeris windows bracket the native exposure times. The 401-pixel cubes cover the inner bands; the 641-pixel cubes cover the 108,000–126,000 km interval.

~~~bash
python moon_spk.py 2460346.55 2460348.26 cache/ura184-2768.bsp cache/ranges
python ephemeris.py Uranus JWST 60346.229 60346.269 cache/horizons cache/2768-eph.npz --ura184 cache/ura184-2768.bsp --window 60347.037 60347.077 --window 60347.690 60347.732
for band in F182M F210M; do
    python neptune_cube.py cache/2768-eph.npz "$band" "cache/2768-$band.npz" "cache/2768-$band.json" cache/uranus2768/*.fits
    python neptune_cube.py cache/2768-eph.npz "$band" "cache/2768-$band-wide.npz" "cache/2768-$band-wide.json" cache/uranus2768/*.fits --size 641
    python broad_scan.py "cache/2768-$band.npz" "cache/$band-inner" --grid 52000 65000 40 1 --support .8 --epoch-jd-tdb 2460347.42
    python broad_scan.py "cache/2768-$band.npz" "cache/$band-perdita" --grid 72000 79000 40 .75 --support .8 --epoch-jd-tdb 2460347.42
    python broad_scan.py "cache/2768-$band-wide.npz" "cache/$band-outer" --grid 108000 126000 100 1 --support .8 --epoch-jd-tdb 2460347.42
    python jwst2768/flanks.py "cache/2768-$band.npz" "cache/$band-inner/passing.csv" "cache/$band-flanks.csv"
done
python jwst2768/miranda.py cache/2768-F182M-wide.npz cache/2768-F210M-wide.npz cache/F182M-clean.npz cache/F210M-clean.npz cache/miranda
for band in F182M F210M; do
    python broad_scan.py "cache/$band-clean.npz" "cache/$band-outer-masked" --grid 108000 126000 100 1 --support .8 --epoch-jd-tdb 2460347.42
    python neptune_cube.py cache/2768-eph.npz "$band" "cache/$band-perdita.npz" "cache/$band-perdita.json" cache/uranus2768/*.fits --control 725
done
python jwst2768/perdita.py cache/phase cache/F182M-perdita.npz cache/F210M-perdita.npz
~~~
