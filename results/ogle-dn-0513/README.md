# A second brightening in the OGLE microlensing sample

The light curve used for the 2019 microlensing event BLG515.30.119116 contains a second brightening 399.07 days after its fitted maximum. The source is also listed as dwarf-nova candidate OGLE-BLG-DN-0513 under the **same OGLE-IV database identifier**. The event entered the sample of 5,790 objects used to construct Galactic bulge optical-depth and event-rate maps. This is a candidate contaminant with an observed second brightening inconsistent with treating the published light curve as a strictly single-bump time series. **The physical origin of either brightening remains open**: a chance microlensing event on a variable star, variability of a blended neighbor, or repeated lensing remain possible. The finding does not by itself invalidate the original lens fit or warrant removing the entry from the published sample.

## What the survey measured

The original event peaks near HJD 2455641.772 and has a fitted Einstein time of 3.144 d. Its brightest recorded point is I = 18.779 ± 0.064 at HJD 2455640.89007. The brightest sampled point during the second episode is on the third night. Sparse sampling leaves its maximum uncertain. The median of 174 microlensing measurements 12–80 days from the secondary center is I = 19.307. The weighted mean of the three consecutive nights is **0.312 ± 0.036 mag** brighter than that baseline. The stated uncertainty propagates the three measurement errors. In the separately published dwarf-nova reduction the corresponding result is **0.335 ± 0.027 mag**. Both reductions use the same exposures.

| HJD − 2450000 | 2019 microlensing I | 2015 dwarf-nova I |
|---:|---:|---:|
| 6039.86793 | 19.040 ± 0.065 | 18.981 ± 0.047 |
| 6040.84406 | 18.988 ± 0.060 | 18.987 ± 0.045 |
| 6041.78795 | 18.958 ± 0.064 | 18.935 ± 0.047 |

The three nights individually exceed the local baseline by 4.11, 5.32 and 5.45 times their formal uncertainties in the microlensing reduction. Every other microlensing event from the same CCD that meets the date and baseline requirements provides a contemporaneous comparison. The three comparison signals, expressed as the sum of the three signed residuals divided by √3, were −0.86, +1.42 and +2.21 versus +8.59 for this source.

We also scanned every eligible sequence of three observed nights outside each object's listed event. The target's strongest window was the one above. Across 37,422 overlapping windows in 64 other light curves from field BLG515, the strongest minimum single-night significance in any window was 2.62. We used the 2019 sample's light curves, grouped by integer HJD, with at most 2.25 days between successive weighted nightly means. We masked at least 25 days or five fitted Einstein times on either side of each event, whichever is larger. At least 30 reference observations were required. This is an empirical comparison against local survey systematics. Overlapping windows and survey selection preclude interpreting their count as an independent-trial significance.

## Which catalog entries coincide?

An all-pairs 2 arcsec match between the 2019 low-cadence sample and the 2015 sample of 1,059 Galactic bulge dwarf-nova candidates returns seven pairs. Shifting every dwarf-nova position east or west by 20, 40 or 60 arcsec returns zero pairs in each control. Every unshifted pair has a brightening in both releases near the published event epoch. Only three have identical field and star IDs. The catalog positions can differ slightly for the same ID between releases.

| Dwarf-nova suffix | Microlensing ID | Angular distance | Same ID | Strongest off-event minimum, in single-night errors |
|---|---|---:|---|---:|
| 0022 | BLG613.07.36212 | 1.416″ | No | 1.25 |
| 0067 | BLG609.29.47693 | 0.019″ | Yes | 1.18 |
| 0083 | BLG630.25.49822 | 1.281″ | No | 0.64 |
| 0135 | BLG639.16.7020 | 0.011″ | Yes | 1.38 |
| 0481 | BLG513.06.80568 | 0.247″ | No, overlapping fields | 2.01 |
| 0513 | BLG515.30.119116 | 0.164″ | Yes | 4.11 |
| 0787 | BLG573.04.19912 | 1.323″ | No | 1.00 |

For 0022, 0083 and 0787 the assigned dwarf-nova source is 3.6–5.2 mag fainter at baseline. A bright neighbor's flux can affect its extraction. The database overlap alone leaves open which source brightened. The other two identical-ID pairs show the same original peak in both reductions and lack a comparably persistent second peak in this search. The long event for 0135 limits searchable off-event coverage.

## Independent infrared observations match the timing for three of seven, without proving lensing

We queried the ESO VIRAC2 archive at all seven microlensing positions. Six have a VISTA source within 2 arcsec. For each we tested the optical catalog’s published peak time, Einstein time and impact parameter against its independent Ks-band light curve, fitting only the infrared baseline flux and source flux. Three sources have infrared nights within one Einstein time. All three brighten on the optical event's timescale.

| OGLE dwarf-nova candidate sightline | VISTA nights inside one Einstein time | Quiescent Ks → brightest sampled Ks | Improvement over constant flux (Δχ²) | Time shifts at least as good |
|---|---:|---:|---:|---:|
| OGLE-BLG-DN-0083 | 2 | 13.305 → 13.015 | 46.7 | 0 / 233 |
| OGLE-BLG-DN-0135 | 15 | 16.721 → 15.705 | 317.7 | 0 / 341 |
| OGLE-BLG-DN-0481 | 5 | 16.988 → 16.191 | 37.2 | 0 / 101 |

For 0135, the strongest test, infrared measurements span −84 to +41 days around the optical peak. The brightest sampled infrared point is 1.016 mag above the off-event median; the deepest optical brightening is about 1.17 mag. The coincidence of timing and shape across two different instruments and bandpasses is consistent with the published single-lens model, but does not rule out a similarly shaped intrinsic outburst. Each infrared amplitude is free because the star and blended light can have different colors. An independent fit to the 0135 infrared data illustrates the remaining ambiguity. A two-sided exponential burst with five free parameters gives χ² = 20.4, compared with 22.8 for the optical lens shape with two free infrared fluxes. The infrared time series validates the predicted brightening and timing. It cannot distinguish these light-curve families by residuals alone. For 0083 the two nights cover only the declining wing. In 0513 the nearest main-catalog VIRAC2 source is 2.52 arcsec away, so the optical recurrence lacks an independent infrared counterpart in this test. The other three sightlines have no infrared epochs within their respective Einstein times.

We combined clean Ks observations by observing night, adding the off-event median-absolute-deviation scatter in quadrature to the nightly photometric error. We repeated the fixed-shape fit at 1,000 evenly spaced epochs over each source’s observing span, retaining epochs farther than three Einstein times from the catalog peak with at least two sampled nights within one Einstein time. Adjacent time shifts are correlated. We used MJD + 2,400,000.5 to compare VISTA with OGLE HJD; the heliocentric correction of at most several minutes is negligible for the 14–99 day events in this comparison. [Smith et al. 2025, VIRAC2](https://doi.org/10.1093/mnras/stae2797) · [ESO VIRAC2 source catalog](https://www.eso.org/qi/catalog/show/424) · [ESO time-series archive](https://www.eso.org/qi/catalog/show/429) · [Independent measurements](vista-results.json)

The 2015 paper included candidates selected by inspection of the OGLE Early Warning System, acknowledged possible microlensing contamination of alert streams by dwarf novae, and included stars with only one detected outburst. The 2019 analysis rejected stars with multiple detected bumps using thresholds in the original extraction, including a requirement that at least three points register on difference images. The second brightening here is measured in the published photometry. Whether it met the pipeline's original difference-image and significance cuts is unknown. Searches on the exact 2019 survey ID, dwarf-nova label and both paper titles on 2026 November 15 located the two catalog records. The source reports inspected did not discuss this specific recurrence.

Removing this one entry as a sensitivity test would reduce the full low-cadence count by 1/5,790 (0.017%) and the sum of catalog inverse-detection-efficiency weights by 10.01/30,131.17 (0.033%). In field BLG515 specifically, 80 events contribute a total weight of 316.63. Keeping the monitored sources, exposure and efficiencies fixed, removing this entry would lower the field event rate by 3.16%, from approximately (5.7 ± 0.8) to 5.52 events per million stars per year. The implied shift is about 0.23 of the published statistical uncertainty. Its contribution to the field optical depth is smaller (0.33%) because of its short timescale. The calculation reproduces the published field's mean Einstein time, 30.5 days, from the best-fit times and published inverse-efficiency weights. The calculation is a sensitivity test. A change to the official sample would require a physical classification. [Original field estimates, Table 7](https://ftp.astrouw.edu.pl/ogle/ogle4/microlensing_maps/table7.dat)

## Reproduce the audit

Run the commands below with Python 3, NumPy and SciPy. The script downloads the two OGLE catalogs and the two public I-band light-curve archives into a local cache, then reports input hashes, matched identifiers, both extractions' peaks and independent-date measurements, the sky-shift tests and the same-CCD controls. The rate weight comes from Table 3. The downloads total about 71 MiB.

```sh
python replay.py > replay-output.json
diff -u results.json replay-output.json
python vista_test.py > vista-replay.json
diff -u vista-results.json vista-replay.json
```

The optical photometry and original optical observations are from the Optical Gravitational Lensing Experiment. VISTA/VVV near-infrared photometry comes from the ESO Science Archive, VVVX DR3 VIRAC2. The [VISTA replay](vista_test.py) reads the archived [query rows](vista/) and also supports `--refresh` to fetch the original data. [Mróz et al. 2015, Acta Astronomica 65, 313](https://acta.astrouw.edu.pl/Vol65/n4/a_65_4_2.html) · [OGLE dwarf-nova archive](https://ftp.astrouw.edu.pl/ogle/ogle4/OCVS/CV/) · [Mróz et al. 2019, ApJS 244, 29](https://arxiv.org/abs/1906.02210) · [OGLE microlensing archive](https://ftp.astrouw.edu.pl/ogle/ogle4/microlensing_maps/) · [Measurements and hashes](results.json) · [Analysis](replay.py)
