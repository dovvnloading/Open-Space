# Blind pulsar timing transient search

**An individual timing search analyzed 68 of 68 pulsars. The array search tested 68 pulsars. 77 template fits survived the specified checks.**

The individual search produced 5395 primary template crossings. Nearby epochs, widths and physical templates can describe the same excursion. The number of fits is not the number of distinct astrophysical events.

Wideband reductions succeeded for 68 pulsars. EPTA DR2 reductions succeeded for 11 overlapping pulsars. The coverage table identifies untested targets. The individual records retain every primary crossing and its disposition. [Coverage](evidence/coverage.csv) · [Per-pulsar dispositions](evidence/individual-dispositions.csv) · [Individual records](evidence/single) · [Surviving fits](evidence/surviving-fits.json)

## Timing signatures and interpretation

The surviving-fits table includes the fitted epoch, amplitude, significance and checks for every surviving template. Individual survivors require independent EPTA support. Array survivors pass the internal leave-one-pulsar-out check. Each array record separately reports wideband and EPTA support, inconclusive checks and missing epochs. A wideband recurrence tests the reduction of the same telescope observations. An EPTA comparison tests independent telescope data at the discovery epoch. Publication history and any physical attribution require an event-by-event review. [Common-mode results and alternatives](evidence/array) · [Published event epochs](PRIOR_EVENTS.md)

## Sensitivity and scope

The individual search reports 368 pulsar, morphology and width combinations with a conditional amplitude limit. Blank limits mark combinations without a calibrated exclusion. The individual records give the injected amplitudes, units, recovery rule and assumptions. The coherent-search tables give injection counts and recovery fractions by signal family. Use the epoch range and sampling distribution accompanying each limit. [Individual sensitivity](evidence/individual-sensitivity.csv) · [Array sensitivity](evidence/array/array_limits.csv)

### Individual amplitude limits

| Timing signature | Width in days | Pulsars with limits | Median limit in µs | Range in µs |
| --- | ---: | ---: | ---: | ---: |
| spin_step_or_memory | — | 13 | 991.9 | 284.7–1.019e+04 |
| shapiro_log | 60.0 | 13 | 165.8 | 27.85–2201 |
| dm_tail | 60.0 | 0 | Unresolved | Unresolved |

For a spin step the amplitude is the timing drift accumulated in one Julian year. The Shapiro coefficient multiplies the logarithmic profile. The DM tail amplitude is the initial delay at 1,400 MHz. The linked data give every tested width and morphology. [Method](SINGLE_SEARCH.md)

### Coherent amplitude limits

| Spatial pattern | Time dependence | First amplitude reaching the recovery criterion | Largest upper bound over eligible grid cells | Units |
| --- | --- | ---: | ---: | --- |
| clock | step | 512 | 22.47 | us |
| ephemeris | step | 16 | 35.59 | us |
| clock | ramp | 128 | 661.3 | us/yr |
| ephemeris | ramp | 128 | 1061 | us/yr |
| memory | ramp | 3.245e-11 | 9.314e-10 | strain |
| clock | gaussian30 | 8 | 38.2 | us |
| ephemeris | gaussian30 | 8 | 43.5 | us |
| clock | gaussian120 | 512 | 1588 | us |
| ephemeris | gaussian120 | 256 | 3259 | us |

The recovery column uses off-grid injections into the measured residuals. The upper-bound column gives simultaneous 95% coverage under the noise model for an exact grid waveform. [Method](ARRAY_SEARCH.md)

## Inputs and reproduction

| Release | DOI | Downloaded archive SHA-256 |
| --- | --- | --- |
| EPTA DR2 | 10.5281/zenodo.8300645 | `2ae546940419ce12a508d0f24b54db91163b0a0d204b85f8bca626eed847b43c` |
| NANOGrav 15-year | 10.5281/zenodo.16051178 | `91476bf20d4c8baa9f5ad39c6e114b581478a5a7e269876489cafc8f62c5708f` |

Individual timing-file hashes, clock corrections and ephemerides are recorded with the archive receipts. [Input receipts](evidence/inputs.json) · [Output hashes](evidence/SHA256SUMS) · [Run revision](evidence/run.json)

Run the manual workflow at the recorded analysis revision with the same software and inputs. The workflow contains the exact fetch, reduction, search, injection and verification commands. [Workflow](../../.github/workflows/nanograv15-transients.yml) · [Reproduction](RUN.md)
