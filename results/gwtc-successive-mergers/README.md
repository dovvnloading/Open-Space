# A search for successive black-hole mergers within the GWTC-5.0 release

Could the remnant of one observed merger be a component of another? Among 104 events with complete mass and distance summaries and sky maps, 5,356 chronological pairs were eligible. The leading match on the primary sky maps is GW240520_213616 followed 32.13 days later by GW240622_004008. The inferred earlier remnant and later primary have nearly equal source-frame masses (18.7 and 18.6 solar masses). Their luminosity-distance medians are 1,400 and 1,470 Mpc. Independently rotating the sky maps produced an equally strong or stronger best pair in 11,181 of 20,000 catalogs. **Null result:** this search finds no statistically unusual association under the stated sky-rotation null; the near-equal masses do not establish that a remnant black hole merged again.

| Sky maps | Sky grid | Leading pair (month/day) | Maximum score | Global randomization fraction |
|---|---:|---|---:|---:|
| IMRPhenomXPHM_SpinTaylor | 256 × 128 | May 20 / June 22 | 31.88 | 0.559 |
| Same maps, refined grid | 384 × 192 | May 20 / June 22 | 32.22 | 0.552 |
| SEOBNRv5PHM | 256 × 128 | June 30 / July 5 | 34.04 | 0.448 |

All four dates are in 2024. The distinct waveform families use the same detector recordings. Their ranking reversal measures model sensitivity.

## Inputs and coverage

The GWOSC release endpoint lists 161 events contributed by GWTC-5.0. Exactly 104 have the four required median parameters and their two-sided uncertainties. Each has both selected sky maps in the stable version 9 Zenodo archive. The other 57 contribute no tested pairs. This experiment covers pairs internal to this release, including the April 2024 inter-run events. Events belonging to earlier releases fall outside the denominator. Event order uses the GWOSC GPS coalescence time. Sky coordinates are equatorial, distance is luminosity distance in Mpc and all masses are in the source frame. Use the raw SHA-256 values in the summary files to check input identity. The skymap archive also has the provider's MD5 digest 6a4b5d6fc3f62d350be9b73dae1796c8.

[LVK data documentation](https://gwosc.org/GWTC-5.0/) · [GWOSC release table](https://gwosc.org/eventapi/jsonfull/GWTC-5.0/) · [LVK version 9 data and citations](https://doi.org/10.5281/zenodo.20348005)

## Ranking and chance control

The code integrates each multi-order HEALPix density onto equal-solid-angle bins in right ascension and sine declination. It uses 3 × 3 subpixel samples per bin (4 × 4 in the refinement) and normalizes each raster. The sky factor is 4π times the integral of the product of the two probability densities. Two uniform densities have a sky factor of one. The saved pre-normalization integral ranges and the finer-grid rerun quantify discretization.

For each mass and distance comparison, divide the difference of catalog medians by the quadrature sum of the two 90% half-widths facing each other, and multiply by 1.644854. Denote these discrepancies by z. Select the more compatible of the two possible later components. Multiply the sky overlap by exp[−(z_mass² + z_distance²)/2]. This is a ranking statistic. The credible intervals are compressed summaries. The calculation neglects sky–distance, mass–distance and mass–spin correlations, signal selection, and the physical prior on merger delays.

Each control applies an independent uniformly distributed, whole-bin right-ascension shift to every event. It preserves each declination distribution, its pairing with that event's intrinsic measurements, and the dependency among pairs sharing an event. Each shifted catalog is searched over all the same chronological pairs and both component choices. The reported fraction includes a plus-one correction. The null assumes right-ascension exchangeability conditional on these quantities. It does not explicitly simulate time-dependent detector duty cycles. A targeted follow-up would require joint posterior samples, spin consistency, and an exposure-aware detection model. The global control gives no reason here to assign a common remnant.

## Reproduce

Use Python 3.12, NumPy 2.3.5, Astropy 8.0.2 and astropy-healpix 2.0.1. The top-20 tables and complete run summaries are in the three output directories. The analysis and data fetch were run on 2026-11-15.

```sh
curl -fLsS https://gwosc.org/eventapi/jsonfull/GWTC-5.0/ -o gw5.json
curl -fLsS https://zenodo.org/api/records/20348005/files/IGWN-GWTC5p0-29ebe06b7_25-Archived_Skymaps.tar.gz/content -o maps.tgz
export OPENBLAS_NUM_THREADS=1
python search.py gw5.json maps.tgz spintaylor
python search.py gw5.json maps.tgz seob --model SEOBNRv5PHM
python search.py gw5.json maps.tgz refined --longitude-bins 384 --latitude-bins 192 --oversample 4
```

## Prior work

Veske et al. previously searched earlier catalogs for black-hole merger families. This adds an independently reproducible GWTC-5.0 within-release result and an explicit all-pairs chance comparison. On 2026-11-15, searches of indexed literature for “successive consecutive gravitational wave mergers same black hole common sky search GWTC-5 O4b” and “hierarchical triple black hole merger pairs in GWTC-4.0 GWTC-5.0” located the earlier search and LVK's catalog release. Publication priority for the general method is established by the earlier paper. [Veske et al., 2021](https://arxiv.org/abs/2011.06591) · [LVK catalog release](https://ligo.org/detections/o4b-catalog/)
