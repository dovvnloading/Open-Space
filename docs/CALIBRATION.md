# Artificial-line recovery on a measured sky background

How bright must a narrow feature be before the pilot detector recovers it on an actual galaxy spectrum? This paired experiment perturbs the calibrated SDSS coadd in memory. Each trial keeps the background flux structure, pipeline inverse variance, masks and wavelength grid from the same observation. The baseline is a full search of the original, untouched spectrum.

Run the experiment from the repository root.

```bash
python -m open_space.calibrate data/sdss-6497-56329-0139.json --trials 24 --seed 17 > recovery.json
```

The command verifies the input CSV digest and pixel count against the manifest. It writes one calibration JSON document to standard output, including input identifiers and digests, search coverage, exact injection locations, outcomes and residuals. It operates independently of the candidate registry.

## Experiment

Draw center pixels uniformly with replacement from positions the detector can search. Each needs twelve contiguous usable samples on either side. The draw ignores flux values, catalog matches and baseline detections. Jitter the centroid uniformly within ±0.45 local pixels and draw a Gaussian width uniformly from 0.9–1.7 local pixels (one standard deviation). Truncate each profile at six profile standard deviations. These widths are specified in native pixels. Instrument-calibrated line-spread functions would support a more representative future population test.

At each drawn position, run six separate searches of the entire spectrum. Add one line per search, using both signs and three peak amplitudes. One amplitude unit is the quoted single-pixel noise at the selected center, `1/sqrt(ivar)`, in calibrated flux-density units. A level of six describes the injected peak before sampling. The detector's six-unit score threshold operates on the three-pixel core and estimated sideband noise. The two quantities have different definitions.

Each search receives the perturbed spectrum and source redshift. Only the post-search scorer knows the injected centroid and sign. A recovery requires a detection of the same sign within two local pixels. A detection of either sign in the baseline within that tolerance of the injection or its matched detection makes the trial confounded. Report all sampled trials, subtract confounded trials from the assessable denominator, and retain the missed trials. Masks and edges excluded 1,335 of this pilot's 4,624 possible center pixels before sampling. This design measures sensitivity conditional on searchable coverage. Sampling center pixels gives each eligible pixel equal weight.

## Seeded pilot result

The command above evaluated six conditions at each of 24 distinct positions in SDSS 6497-56329-0139. The original spectrum had 12 baseline features and 3,289 searchable centers. All 24 trials in each condition were assessable.

| Injected peak in local noise units | Emission recovered | Absorption recovered |
| ---: | ---: | ---: |
| 3 | 0 / 24 | 1 / 24 |
| 6 | 20 / 24 | 19 / 24 |
| 10 | 24 / 24 | 23 / 24 |

One absorption injection at ten noise units still fell outside the recovery criteria on this observed background. The six outcomes at each position share the same noise realization. These counts describe one coadd and one seeded set of positions. Survey completeness needs independent targets, wavelength and width strata, and instrument response. Blank-sky fibers and independent exposure differences are separate controls for spurious detections. [Full science protocol](SCIENCE.md)
