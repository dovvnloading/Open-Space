# VIS morphology screen and injections

Run with Python, NumPy, SciPy and the system CFITSIO shared library. The summary additionally requires Astropy. Coordinates in the stamp catalogs are zero-based with the first array index along y. Segment orientation is measured from image +x toward +y. The summary translates both endpoints and the midpoint through the released FITS WCS.

## Search

Run from the directory containing this file. `INPUTS` is the cutout cache from the Q1 acquisition script, `TILES` is its frozen tile CSV, and `RUNS` is an empty or matching per-stamp checkpoint directory.

```bash
OPENBLAS_NUM_THREADS=1 python screen_centers.py --input-dir "$INPUTS" --tiles "$TILES" --out "$RUNS" --jobs 4
OPENBLAS_NUM_THREADS=1 python calibrate_centers.py --input-dir "$INPUTS" --tiles "$TILES" --out "$RUNS" --jobs 3
python summarize.py --runs "$RUNS" --images "$INPUTS" --out tables
python test_wake.py
```

1. The Q1 mask accepts finite nonzero VIS, with finite positive RMS below 10⁸. The RMS identifies coverage. The detection score uses the dispersion of the filtered coadd. A mesh of 48-pixel median sky estimates is bilinearly interpolated and subtracted.
2. A Fourier-domain matched filter uses axial top hats of 6 and 10 arcsec, an observed transverse Gaussian FWHM of 0.20 arcsec, and negative flanking templates displaced by 2.5 FWHM. The tested unoriented angles uniformly span 180 degrees, with angular spacing at most 1.25 FWHM/length radians. For each filter the response is centered on its median and divided by the interquartile range / 1.349 on a grid of valid interior pixels. A threshold of six is a ranking cut. The correlated noise and trials over position, angle and scale preclude a Gaussian false-alarm interpretation.
3. Local maxima of the strongest response over all templates seed profile measurements. Bilinear sampling compares the axial trace to lateral backgrounds. At least four of five longitudinal bins must contain positive excess above 15% of the central smoothed median. The contiguous extent ends at 35% of that median. The transverse mean profile gives an observed FWHM by interpolation between half-height crossings, floored at one pixel. Retain measured length at least 3 arcsec, at least 70% of the seed length, and measured length/FWHM strictly above 20. A boundary-censored measurement is flagged. Approximate duplicates at similar angles, lengths and centers within the same stamp are suppressed.
4. An exploratory follow-up gate was added after seeing numerous compact-source-centered detections in the first four stamps. Divide the measured segment into eight bins. At least six must have a center brighter than both local flanks. The brightest bin may contribute at most 40% of the sum of positive bin excesses. All sampled pixels must be valid. This tests the distribution of light along the segment. It cannot identify a host galaxy or distinguish an instrumental spike, trail, edge-on disk or chance alignment.
5. Repeat the same search on the sign-inverted image with identical mask and cuts. This control responds to noise fluctuations and negative structures. Counts from the control cannot be subtracted as a complete estimate of false positives from bright astronomical sources. `segments.csv.gz` retains both signs, including rejected profiles.

The coverage table counts all valid, geometrically disjoint science pixels. For a conservative center-area diagnostic it also counts pixels farther than 100 pixels (10 arcsec) from every mask or exterior pixel. This margin contains the longest profile in any orientation. The actual accepted set depends on direction, length, brightness and neighbors. Neither area is an efficiency-corrected survey area.

## Injection design

The model is an observed, straight, uniform-length filament with a Gaussian transverse profile. Its FWHM is 0.20 arcsec. Mean AB surface brightness refers to the flux density per square arcsecond of the model averaged over its transverse FWHM strip. Numerical integration uses 4 × 4 subpixels. A point-source AB zeropoint of 24.6 and a 0.1 arcsec pixel imply mean native counts per pixel equal to 0.01 × 10^(−0.4 × (μ − 24.6)). Injections do not add photon-counting noise, which limits application to bright-source regimes.

Each trial samples a 192 × 192 real VIS patch and a random center within a fixed boundary margin, with an unoriented angle drawn uniformly. The entire patch undergoes the blind search before and after one insertion. Adjacent magnitude bins use the same sites and angles. `new_match` requires a post-injection detection consistent with the known location, angle and length, and no matching raw detection in the original patch. `new_profile_match` adds the eight-bin gate. The patch has its own background and noise estimate. It approximates the full-stamp search and adds local boundary effects. Trials measure morphology recovery for the inserted family. They do not measure the efficiency of host association, rejection steps or spectroscopic confirmation.

The pilot uses 12 positions per length in the three nonempty large stamps. The fourth large center has no valid VIS pixels. The second sample uses eight additional centers per deep field, selected by sorting the SHA-256 of `wake-centres-1/{tile}`, with two positions per length in each. Empty centers remain eligible in this sample. Sampling at each tile is seeded by its numeric tile ID. `recovery_by_sample.csv` distinguishes all attempted positions from the subset for which the injected axial strip and flanking region (±4 FWHM) have fully valid sampling. Repeated brightness trials, and sites sharing a stamp or patch, are correlated. A zero in a finite trial cell is an observed count.

The rows in the shortlist are photometric segments for further inspection. The data establish neither projected physical length nor a common redshift of a putative host and filament. Continuity, gate failure or success, or an isolated imaging score cannot establish the origin.

## Shape measurements and inspection

`morphology.py` measures the distribution of flux in eight bins on each gated segment, using a separate flanking sky estimate in each bin. Its center-of-light coordinate, central-half flux fraction and reversal asymmetry describe the axial shape. `centerline_scatter_in_fwhm` is the RMS residual of positive-light transverse centroids about a fitted line, in units of the measured width, requiring five usable bins. Source confusion and noise affect all of these descriptors. They are reported without disk, stream or jet classifications.

`montage.py` uses the Euclid host and filament tables to render the host-matched, VIS-continuous positive subset. Each panel pairs a full-stamp image with a local view. Cyan rings mark measured endpoints and an amber ring marks a catalogued host proxy. When an endpoint lies short of a catalogued object on the same axis the marker follows the reported axial gap and cross-track offset. An asinh stretch centers on the local median sky and uses the median absolute deviation for display scaling. The source counts and measurements come from the calibrated FITS arrays.

```bash
python morphology.py --runs "$RUNS" --images "$INPUTS" --out tables/shapes.csv
python montage.py --hosts ../euclid/endpoint_candidates.csv --photometry ../euclid/filament_photometry.csv --images "$INPUTS" --out review
```
