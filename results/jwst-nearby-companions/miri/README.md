# Two-epoch MIRI white-dwarf controls

No peak in the 4–12″ annulus passed the two-epoch common-motion test for either white dwarf. This is a new reduction of public mosaics of fields with a published companion search. At 5″, 8″ and 11″, the faintest injected F1500W companion/host ratios reaching at least 90% recovery were 0.010, 0.010 and 0.040 around WD 1202−232, and 0.040, 0.020 and 0.020 around WD 2105−82. These are the first qualifying points on a discrete brightness grid, with 40 sampled angles per radius.

| Field | F1500W epochs | Early peaks | Stationary matches | No qualifying association | Stationary offset relative to host (east, north) | Registration RMS |
|---|---|---:|---:|---:|---|---|
| WD 1202−232 | 2023-02-09 and 2025-06-03 | 11 | 4 | 7 | (+0.037, −0.587)″ | 0.039″, 33 ties |
| WD 2105−82 | 2023-04-21 and 2025-06-03 | 5 | 5 | 0 | (−0.559, +0.844)″ | 0.031″, 17 ties |

## Measurement

- Frozen input comprises the two public Level 3 mosaics per star (GO 1911 and 4857). We reuse the repository's repeat-imaging product URL builder. Original CRDS context, pipeline version, native host centers and noise are in the [measurement record](measured/measurement.json). The 2025 measurements were released in June 2026. [Archive query and full-file digests](inputs.json)
- Reproject each image to a tangent grid centered on its measured host core, reject low weights and dilate missing coverage. Search a 1.7-minus-5.5 pixel Gaussian-filtered image for peaks above six local MAD-derived standard deviations in each epoch. Fit subpixel peak locations. The empirical stationary offsets include host proper motion, parallax and translation of the frames. Require a late peak within 0.20″ of the early host-relative position and a residual at least 0.15″ smaller than the closest stationary prediction. Denominators include every early trigger. The seven unmatched peaks have no assigned astrophysical classification. [All early triggers and both hypotheses](measured/sources.csv)
- Scale each epoch's measured stellar core by a trial companion/host ratio at the same host-relative offset. The template tapers from 0.65 to 0.90″. Apply the identical linear filter, peak, match and stationary-confusion cuts with baseline noise fixed. Count a trial only if both fitted positions lie within 0.20″ of the injection and the pair was absent before injection. Each angle recurs at eight trial ratios; neighboring trials can be correlated. The test is at the mosaic stage and does not calibrate original detector processing, orbital motion, colours, or masses. [1,920 individual trials](measured/injections.csv) · [Contrast and recovery](measured/contrast.csv)
- For each peak, apply the 70%-encircled-energy MIRI imaging corrections to sky-subtracted native mosaic apertures in all four first-epoch filters. Measure local scatter with 48 control positions. Fit the signed fluxes to a free nonnegative scale across the solar-metallicity, cloud-free 200–1200 K Sonora Bobcat grid, allowing an extra 10% flux-dependent uncertainty. The tabulated minimum weighted residual ranks model colours. Blending, non-Gaussian noise, imperfect aperture corrections and limited atmospheric chemistry prevent a physical temperature or classification from these fits. [Forced colours and grid comparison](measured/colors.csv) · [Marley et al. model data](https://doi.org/10.5281/zenodo.5063476)
- In the 15–29″ control annuli, 34 and 18 early peaks have late stationary matches within 0.20″. Scaled by area, the counts predict about 7.1 and 3.7 stationary matches in the respective search annuli. This includes compact galaxies and stars; bright-star diffraction and spatially varying completeness affect the counts. Two stationary peaks have epoch-propagated Gaia matches within 1″. The NASA confirmed-planet cones contain no rows. The SIMBAD candidate labels in both fields share coordinates with their hosts. [Per-peak catalogue checks](catalogs/crossmatch.csv) · [Field counts](measured/background_counts.csv)
- Mullally et al. previously measured the background motion of both *inner* resolved candidates. Those published detections lie inside our 4″ inner radius. The central infrared excess around WD 2105−82 also lies outside this off-axis selection. [Initial report](https://arxiv.org/abs/2401.13153) · [Two-epoch follow-up](https://arxiv.org/abs/2512.08191)

## Reproduction

Run in this folder, keeping the repository's sibling repeat-imaging campaign available. NumPy, SciPy, pandas and Astropy are required. Downloaded FITS and calibration inputs are verified against SHA-256. Run catalogue checks offline against the frozen answers or use a fresh copy and pass `--refresh` to requery.

```sh
python fetch_inputs.py /tmp/jwst-wd
python search.py /tmp/jwst-wd /tmp/jwst-wd-results
python vet.py /tmp/jwst-wd /tmp/jwst-wd-results
python catalogs.py /tmp/jwst-wd-results catalogs
```
