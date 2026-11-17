# GJ 205 flux scale and published searches

HD 36395 = HIP 25878 = GJ 205 is a nearby M1.5Ve star. The W3-anchored grid gives 1.263 Jy in F1065C and 0.603 Jy in F1550C (grid medians). These host fluxes set the denominator for the contrast ratios. The Ks-anchored ranges are 1.041–1.233 Jy and 0.495–0.593 Jy respectively. For F1065C the W3 grid itself spans 1.256–1.268 Jy. These are alternative assumptions about the stellar flux scale, without a probabilistic interpretation.

AllWISE J053127.91−034059.8 has W3 = 3.576 ± 0.031 Vega and no flagged W3 saturation fraction, but the catalogue flags a halo in the W3 band. 2MASS Ks = 4.039 ± 0.260 for 05312734−0340356. The alternate normalization is consequently uncertain. The colour calculation integrates the same PHOENIX spectra and CALSPEC Vega reference as the NIRCam target normalizations, with temperatures 3500–3900 K and log g 4.5/5.0. The synthetic photometry is based on the SVO transmission curves. [Per-model flux scales](photometry.csv) · [Catalogue queries](catalog_inputs.json) · [Bandpass digests](bandpasses.json)

Boehle et al. previously searched for companions around HD 36395 with high-contrast imaging and radial velocities. GO 9056 is an imaging search targeting accelerating stars. These historical observations and an acceleration do not identify the present field source; its Gaia coordinates, parallax and proper motion do. [Previous target-specific search](https://arxiv.org/abs/1907.04334) · [Gaia relative astrometry](../../vetting/field_identifications.json)

Run `python normalize.py /tmp/photometry-cache` to replay the W3/Ks normalization. The underlying spectra are fetched and checked against the sibling stellar input manifest.
