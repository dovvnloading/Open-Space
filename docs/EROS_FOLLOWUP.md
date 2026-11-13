# Resolved ground-based photometry of the 1.115-day LMC field

Open-Space detected a localized signal at 1.11221 and 1.11694 days in TESS sectors 1 and 2. The EROS-2 all-star archive contains **115 resolved source records and 62,748 paired red/blue exposure rows** within 30 arcsec of the discovery position (77.0550166°, −69.2370982°). The sector 2 centroid lies 6.30 arcsec from that position. No EROS source repeated a strong signal near the TESS period in both later red and later blue observations at the frequency selected from earlier red measurements. The EROS observations in this cone span **July 1996 to January 2002**, more than 16 years before TESS sector 1. [Source-by-source measurements, frequencies and checksums](../results/ffi-ground-followup.json) · [TESS field measurements](../results/ffi-lmc-grid.json) · [EROS data release](https://doi.org/10.1051/0004-6361/202557321)

EROS catalog positions are J2000; the search used the original reference-image astrometry without proper-motion propagation. We verified that the row count for each source exactly equaled its catalog exposure count. Each exposure can carry a missing or unreliable measurement in one or both filters. Of the 115 records, 77 had sufficient red measurements, 89 had sufficient blue measurements, and 76 had adequate paired coverage for the split-epoch test. The nearest record (1.96 arcsec, red reference magnitude 19.72) retained only 77 red measurements. A bright source 2.40 arcsec away (red 14.85) retained 423. Multiple stars can contribute to one 21-arcsec TESS pixel. [EROS source catalog and time series](https://cdsarc.cds.unistra.fr/viz-bin/cat/II/390)

## What the comparison measured

We used published positive magnitude errors up to 0.3 mag, added a 0.01-mag error floor in quadrature, subtracted per-season medians at gaps over 75 days, and clipped outliers beyond six times the larger of the robust scatter or measurement error. Each band needed at least 120 retained measurements over one year. The target frequency interval is the mean of the two TESS frequencies ±1/27 cycles per day (periods 1.07038–1.16256 d). This covers the frequency resolution of the short TESS sectors. EROS times are its heliocentric Julian dates, offset by 2,449,999.5.

The first temporal half of the red time series selected a frequency from that interval. At the locked frequency we measured variance explained, semiamplitude and phase in the later red and simultaneous blue data. The highest value of the **smaller of the two held-out variance fractions** across all 76 records was 3.2%. That source, lm0101l4872, is 14.61 arcsec from the TESS position. Its independently selected full red and blue periods in the window differ (1.09273 and 1.13666 d). The evidence retains each source's wider periodogram maximum and first solar-day and sidereal-day alias windows. These diagnostics are sensitive to sampling and seasonal offsets.

We injected a persistent 1.11457 d sinusoid into the cleaned, observed magnitude series with the same amplitude and phase in both bands. Each source used one deterministic phase, shared across its three trials. Recovery required at least 5% variance explained in the early red and both held-out bands, phase agreement within 0.2 cycles and a frequency within one early-sample resolution element of the injected frequency.

| Injected semiamplitude in each band | Recovered among 76 sources |
| ---: | ---: |
| 0.01 mag | 7 (9%) |
| 0.03 mag | 24 (32%) |
| 0.10 mag | 68 (89%) |

All three injections were recovered for lm0101l5232, the bright source 2.40 arcsec away. The field-wide test has limited sensitivity to small amplitudes and weak/per-band signals. The quoted fractions refer to these 76 observed records, injected phases, cadences and errors. They supply no uniform upper limit for the TESS source. Episodic signals, blending, changing color amplitude and the 39 records outside paired coverage remain plausible explanations of a missing match.

## Two previously cataloged neighbors

We independently fitted the published OGLE I-band measurements for the two nearby variables found by the six-family catalog screen. The source 13.41 arcsec away, OGLE-LMC-RRLYR-05589, has 774 epochs. Its cataloged 0.2476053 d period explains 27.5% of the weighted variance in a one-term sinusoidal fit. The strongest sinusoid in the TESS window explains 1.52%. For OGLE-LMC-ELL-01519, at 22.38 arcsec, the best TESS-window fit explains 0.48% over 5,778 epochs. The two OGLE time series span 2,247 and 5,196 days respectively. The seasonal offset correction and single-sinusoid model suppress part of the ellipsoidal variable's 220-day waveform. Neither light curve supplies an independent 1.115-day identification. [OGLE catalog associations](../results/ffi-ogle.json) · [RR Lyrae photometry](https://ogle.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/rrlyr/phot/I/OGLE-LMC-RRLYR-05589.dat) · [Ellipsoidal variable photometry](https://ogle.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/ecl/phot/I/OGLE-LMC-ELL-01519.dat)

## Reproduce

Run from the repository root with the numerical dependencies installed. The script verifies the source roster and number of exposures for every source before producing the result. The JSON includes hashes of the actual returned EROS tables and OGLE files. VizieR inserts its request date in the table header, so independently retrieved equivalent rows can have different raw response hashes.

```bash
python examples/ffi_eros_followup.py --cache cache/ground \
  --evidence results/ffi-lmc-grid.json --output results/ffi-ground-followup.json
```
