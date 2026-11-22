# Reviewing the blind-search crossings

The first reduction analyzed 62 narrowband targets. Several angular patterns crossed the array threshold near the same epochs. That result motivated the following checks. They are post-discovery tests of the fitted patterns. The preregistered primary trial statistics, their model-conditional probabilities and all rejection records remain available.

## Coherence under pulsar sign changes

A single random sign per pulsar reverses or preserves its entire residual time series. All epochs, widths and angular models in one draw use the same sign for that pulsar. Its cadence, relative amplitudes, internal correlations and transients are preserved. Comparing each observed statistic to the maximum in each draw tests whether the observed alignment of signs is unusual over the whole bank. This check assumes independent sign symmetry across pulsars. A true common Gaussian background can violate that assumption. The primary calculation includes a specified covariance allowance for that background.

The diagnostic-subset check repeats each fixed template using pulsars that passed the individual residual-noise tests. Per-pulsar contributions and system membership help identify cases in which a few sources or instruments supply most of the fitted amplitude.

## Dispersion and instrument dependence

A change in the distribution of observing frequencies can make dispersive variations project into the achromatic fit. The fixed-epoch control fits an independent inverse-frequency-squared delay within a telescope visit wherever multi-frequency observations identify an achromatic intercept. It propagates measurement covariance and the timing model into that intercept. Inspect how much sensitivity remains and whether the amplitude is compatible with the primary result. Fixed Arecibo-only and Green Bank-only comparisons test instrument dependence when their sky coverage is adequate.

## External evidence and interpretation

Wideband processing tests the same NANOGrav observing sessions. EPTA contributes independent telescope measurements at the frozen epoch, shape and sky direction. Both arrays use DE440, so an EPTA recurrence would still share a solar-system ephemeris. For each internally surviving fit, report coverage, signed amplitude, uncertainty, directional support and multiplicity across the external checks. Record ambiguous geometry and missing support separately. Individual and common fits with different labels can refer to the same timing excursion.

The adjudication looks for agreement across the applicable controls and for previously reported timing behavior. It does not recalibrate discovery thresholds using the selected peaks. The amplitude limits retain the precise coverage and recovery definitions of the individual and coherent search protocols. Passing an internal array threshold alone cannot establish a clock fault, ephemeris error or gravitational-wave burst. [Individual protocol](SINGLE_SEARCH.md) · [Array protocol](ARRAY_SEARCH.md) · [Published comparators](PRIOR_EVENTS.md)
