# WISE checks of the foreground trajectory

Independent WISE infrared positions also carry the T8 dwarf past both Gaia stars. They change which annual-parallax minimum is deepest for the brighter star. Between the two WISE data products and two fitting choices, its minimum separation is 180–263 mas in October 2033. The fainter star reaches 62–216 mas in late October or November 2039. The projected track is sensitive to the infrared astrometric solution well beyond the formal error on the Spitzer fit. The range across these four fits measures model sensitivity, without assigning probabilities to the alternatives.

| Lens track | Gaia DR3 2060116591502525952 | Gaia DR3 2060116595791115776 |
| --- | ---: | ---: |
| Published five-parameter Spitzer fit | 140 mas, Nov 2034 | 65 mas, Jul 2039 |
| unTimely W2, parallax fitted | 180 mas, Oct 2033 | 62 mas, Oct 2039 |
| unTimely W2, Spitzer parallax and scan-direction offsets | 204 mas, Oct 2033 | 117 mas, Nov 2039 |
| NEOWISE L1b, parallax fitted | 247 mas, Oct 2033 | 199 mas, Nov 2039 |
| NEOWISE L1b, Spitzer parallax and scan-direction offsets | 263 mas, Oct 2033 | 216 mas, Nov 2039 |

## Measurements

The unTimely W2-only check uses 16 image-coadd positions from May 2010 through October 2020 on tile 3003p363. The NEOWISE Reactivation single-exposure product fits source positions jointly across W1 and W2. The cuts retain 343 detections across 21 sky visits from May 2014 through May 2024. Exposure overlap between these two WISE products prevents treating their results as independent corroborations of each other. Spitzer observed the foreground dwarf with a different telescope. [unTimely catalogue and extraction](https://arxiv.org/abs/2209.14327) · [NEOWISE single-exposure documentation](https://irsa.ipac.caltech.edu/data/WISE/docs/release/NEOWISE/expsup/sec2_1d.html)

We select reference stars with Gaia G brighter than 16.5, RUWE below 1.4, a published parallax and no second Gaia source brighter than G=18.5 within 6.5 arcsec. For each date we propagate the local stars with their Gaia motion and parallax, then subtract the median east and north residuals from the dwarf position. Nineteen stars supply the W2 frame, and twelve supply the L1b frame. Remaining faint blends and color-dependent astrometric offsets can affect the dwarf differently from the local stars. WISE visits occur at two nearly opposite parallax phases, which makes parallax particularly sensitive to scan-dependent biases.

For W2 the free fit gives proper motions (+12.9, +377.0) mas/yr east and north, with formal errors of 5.6 mas/yr on both and parallax 177.9 ± 16.7 mas. L1b gives (+13.1, +386.6) ± 4.8 mas/yr and 149.6 ± 14.6 mas. Spitzer gives approximately (+6.0, +372.9) mas/yr and 133.7 mas. The fit uses equal weight per sky visit. Quoted errors assume independent residuals of equal variance. The sensitivity fits fix the lens parallax at 133.4 mas and fit east and north offsets of opposite sign for the two scan directions. All fits propagate the same background Gaia source solutions, adopting zero parallax for the source with a negative measurement.

## Reproduction

With Python, NumPy, pandas and Astropy, run `python compare.py` from this folder to fit the frozen epoch positions, register them using the 551 source-by-visit offsets and print the encounters. Coordinates in the epoch tables measure milliarcseconds east and north of (300.209094°, +36.497796°). The parallax factors are dimensionless. The formal-error printout for the fixed-parallax model lists the two scan offsets after the four motion parameters.

To rebuild the small inputs from the public products, add SciPy and run `python extract.py /tmp/wise-cache /tmp/wise-epochs`. The manifest gives exact archive URLs and hashes for the 16 W2 catalogues. The extraction script includes the full Gaia and IRSA field queries and applies the source quality cuts before computing visit medians. Both routes reproduce the listed fits.
