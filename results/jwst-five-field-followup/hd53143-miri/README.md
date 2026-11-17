# HD 53143 at 25.5 µm

Neither of the two November 9, 2025 rolls produced a positive point-source threshold crossing between 3″ and 9″. Their joint mosaic also produced zero crossings. All three products contain the same observing epoch. These are limits on point sources in the specified footprint.

## Inputs and method

The [MAST query and receipt](receipts/query.json) and [response bytes](receipts/caom.json) identify GO 8142 observations 003 and 004, their combined mosaic, and a separate offset background pointing. The [manifest](inputs.json) hashes the three science files and the two F2550W BD+60-1753 calibration mosaics from the merged analysis. Both calibration images were verified against its original hashes. The G29-38 analysis routine was called with the F2550W band, the 0.83″ FWHM, each science mosaic in turn and the same two calibration PSFs. The original G29-38 command retains its original targets and bands.

The filter is normalized to the host's 1.33″ core aperture. The detection criterion uses the one-sided Gaussian 5σ false-alarm tail with a small-sample Student correction in each annulus. An injection must produce a new local maximum within half a FWHM at a supported position. The PSF for injection and the matched kernel come from different standard-star visits. The reported recovery fraction is empirical across these position angles at a discrete tested flux level. The two rolls have different detector orientations but share astrophysical structure, including disk emission. No total-flux correction, planet-mass conversion or combined-epoch motion limit is inferred.

## Recovery

| Radius | Roll 003 | Roll 004 | Combined mosaic | Combined usable annular area |
|---|---:|---:|---:|---:|
| 3.75″ | 0.01911 | 0.01931 | 1.558 (4 angles) | 12.7% |
| 5.25″ | 0.009173 | 0.009196 | 0.01193 (39 angles) | 72.5% |
| 6.75″ | 0.003258 | 0.003281 | 0.003190 (48 angles) | 100% |
| 8.25″ | 0.008190 | 0.008196 | 0.008064 (48 angles) | 100% |

Entries give the smallest tested companion/host core flux ratio reaching at least 95% recovery. Each single-roll annulus passed the pixel support rule over its full area with 48 sampled angles. The combined inner annulus has little coverage under that same weight rule and yields no useful faint-companion sensitivity at 3.75″. The single-roll thresholds are the appropriate inner limits. The output is [per mosaic](measurements/) and includes all 3,138 injection trials and full-image versus linear-injection checks.

## Replay

From this directory, run `python run.py /path/to/fits-cache` with NumPy, SciPy and Astropy. The driver downloads missing files, verifies their SHA-256 and writes the measurements. To reproduce the archive selection independently, run `python freeze_inputs.py /path/to/fits-cache --freeze` in a clean working copy. Later archive releases can change the response.
