# SkyMapper g and r photometry of two infrared-red stars

The two red infrared-excess stars show no 0.20-mag faint excursion relative to their own sampled g/r medians. Every retained night has a nominal pointwise 3σ upper limit below that level. For FGK31 the dates span 2014-03-19 through 2021-09-11, with 12 distinct observing nights across both bands. For FGK18 there are eight distinct nights between 2014-11-05 and 2020-02-20.

| Star | Band | Detections before / after cuts | Nights | Faintest residual | Largest pointwise +3σ upper limit |
|---|---|---:|---:|---:|---:|
| FGK31 | g | 13 / 11 | 11 | 0.036 mag | 0.116 mag |
| FGK31 | r | 12 / 10 | 10 | 0.064 mag | 0.167 mag |
| FGK18 | g | 8 / 6 | 6 | 0.027 mag | 0.101 mag |
| FGK18 | r | 7 / 6 | 6 | 0.010 mag | 0.094 mag |

Every accepted exposure had at least three comparison stars on the same image and CCD. All six TAP queries succeeded on the official Australian archive. Both counterparts have 0.10″ SkyMapper-to-Gaia separations. FGK31's second-nearest catalogued Gaia source is 7.78″ away. Its largest g-band gap is 1,041 days; FGK18's is 761 days. [Nightly measurements](nights.csv) · [Full summary](summary.json)

The immutable targets are Gaia DR3 3017430205815280384 (FGK18) and 6359523376149191936 (FGK31). A unique SkyMapper DR4 master source must have that Gaia ID as its closest counterpart within 1.5″. Each exposure must be within 1.2″ of the Gaia position propagated to that exposure, have zero Source Extractor flags and bad pixels, an image quality of 1 or 2, and a finite PSF error ≤0.15 mag. The test retains otherwise-good points clipped by the SkyMapper mean-magnitude pipeline so that clipping alone cannot hide a dip. Gaia and SkyMapper separation/neighbor flags are recorded.

Use g and r separately. For up to 16 field stars within 0.12° matched in brightness and color, form each reference's offset from its own median. An exposure with ≥3 usable same-CCD references is corrected by their median offset. Other quality exposures retain their survey calibration. Nightly magnitudes are medians. Errors include the measurement, comparison offset and a 0.02-mag floor. Nightly uncertainty conservatively uses the larger of the largest exposure error and the robust inter-night scatter, with an allowance for finite-sample baseline error. Report raw detection counts, retained exposures, distinct observing nights, gaps, maximum faint excursions and pointwise 3σ upper bounds relative to a sampled median. This establishes constraints only on the catalogued, sampled detections. It cannot exclude an eclipse during gaps, a disappearance below the detection threshold, an event diluted within an exposure or a constant attenuation. No duty-cycle bound follows.

Run `python measure.py` in this folder with NumPy and Requests. Raw catalog answers and their hashes, full queries and timestamps are in `raw/` and `receipts.json`. The script verifies and reuses saved responses. A summary marked incomplete records a failed request separately from zero detections.

[SkyMapper DR4](https://doi.org/10.25914/5M47-S621) · [DR4 photometry and image metadata](https://skymapper.anu.edu.au/table-browser-standalone/dr4/) · [Instrument quality and flags](https://skymapper.anu.edu.au/data-release/dr4/)
