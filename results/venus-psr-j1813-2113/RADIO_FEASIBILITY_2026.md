# Follow-up V2: Warkworth is inside the Venus path but its documented radio sensitivity is insufficient for a baseline

**Additional independent observability result (14 November 2026)**. This addendum supersedes any implication in the earlier [independent geometry follow-up](INDEPENDENT_FALSIFIABILITY.md) that Warkworth's documented 12-m L-band hardware has confirmed *current* operational status. The ephemeris/chord result is unchanged. An occultation prediction must be distinguished from a detectably pulsing occulted source.

## Hardware provenance: older and newer reports are not identical

[Warkworth's IVS 2021–2022 station report](https://ivscc.gsfc.nasa.gov/publications/br2021%2B2022/NASA-TP-20230014975.pdf) explicitly reports an **L-band feed at the 12-m antenna's prime focus** as well as an S/X feed at its secondary focus. That does not establish that L-band operation, a sensitive L-band receiver, pulsar-appropriate backend and scheduled time remained available in 2027.

The newer, **March 2026** [NASA International VLBI Service 2023–2024 report, *Warkworth Space Center*, pp. 94–97](https://ivscc.gsfc.nasa.gov/publications/br2023%2B2024/NASA-TP-20260002317.pdf) instead lists a 12-m **S/X** operational chain, 30-m **S/C/X**, and:
- Both telescopes **S-band 2200–2300 MHz**, only **100 MHz** wide in that published receiver summary.
- **30-m S-band SEFD 2,100 Jy**, with two circular-polarization receiving channels.
- **12-m S-band SEFD approximately 4,800 Jy** following interim LNA upgrade; 5,800 and 8,900 Jy in earlier intervals.
- Both elevation mount limits (12-m 7.2°–88°; 30-m 6.0°–90.1°) allow the ~29° geometric Venus elevation predicted at Warkworth.
- **No 12-m L-band noise/SEFD figure**, active feed/receiver status or pulsar-folding capability is provided. The absence of L-band from the newer table does **not** prove the prime-focus feed was dismantled. Current L-band feasibility is **unknown**.

## Independent radiometer test: S-band only

From the [2024 MeerKAT timing/flux release](https://doi.org/10.5281/zenodo.8430591), the pulsar's **L-band** median flux is **0.479 mJy** across 21 available measurements (range 0.208–0.570 mJy). There is **no measured S-band flux here**. For illustration, optimistically assume that the *entire* L-band median flux remains at 2.25 GHz, that the telescope is tracking, the whole 100 MHz is usable, both polarizations are independent, the pulse's effective duty cycle is 10%, and the baseline can be integrated for six minutes outside the occultation. The standard periodic radiometer equation is

```
SNR = (flux_Jy / SEFD_Jy) × sqrt(npol × bandwidth_Hz × integration_s)
      × sqrt((1 - pulse_duty_fraction) / pulse_duty_fraction)
```

| Fully specified hypothetical 2.25 GHz scenario | Six-minute ideal folded-pulse SNR |
| --- | ---: |
| 30-m at 2,100 Jy, assumed flux 0.479 mJy and 10% duty | **0.184** |
| 12-m at 4,800 Jy, same flux and duty | **0.080** |
| *Very favorable*: 30-m, assume flux equals the brightest available L-band visit (0.570 mJy), 1% duty, and full 100 MHz | **0.73** |

This is a **conditional non-detectability result** for the published S-band station, under explicit spectral and pulse-width assumptions. Under the reference scenario, a 5σ pulsed baseline in six minutes requires **~13.0 mJy at 2.25 GHz**, about **27 times** the source's historical *L-band* median; alternatively, a telescope with **SEFD ≤77 Jy** at the assumed 0.479 mJy, bandwidth and duty is required. Neither threshold is currently established at the geometrically suitable Warkworth site. A 10σ detection imposes twice the flux or half the SEFD.

SNR here is not an occultation significance prediction: it describes only one off-occultation pulse baseline. A robust disappearance experiment needs reliable pulses before and after the event, so low single-baseline SNR is an observational blocker, not evidence that the planetary crossing is false.

## Exactly what survives and what would falsify it

**Survives**: independent Moshier vs DE441 agreement; 2024 position geometry for a ~6-minute Warkworth occultation; 2023 position predicts a Warkworth miss; demonstrated inability of MeerKAT/ASKAP to see this event (below horizon), Parkes/ATCA to track it (mechanical elevation).

**Does *not* survive as a justified claim**: that Warkworth's **documented S-band receivers** can detect a ~0.5 mJy pulsar in the few-minute window, or that the older L-band feed is presently operable with known pulsar sensitivity. Also not justified: an S-band flux extrapolation as an *observed* upper bound; a pulsar spectral turnover, scintillation or mode change could alter flux either way.

**Falsifiable conditional instrument result**: A measured Warkworth **S-band** pulsed flux substantially >13 mJy, or a materially lower verified SEFD/higher bandwidth, would defeat this 5σ example under its stated duty/time assumptions. The current instrumental status can be falsified by a contemporary engineering report establishing operational L-band pulse folding (and calibrated L-band SEFD) on the 12-m or another geometrically suitable station. A direct positive six-minute 5σ pulsed baseline from the actual telescope supersedes this budget altogether. Repeat it *both* before and after the event.

**Next necessary measurements**: obtain pulsed flux/spectral shape at Warkworth's available band *before* 2027-02-06; measure the working backend's actual L- or S-band SEFD, usable bandwidth, RFI/solar interference and pulse duty. Acquire a new subarcsecond **pulsed** position at event epoch. Only if source position, site visibility *and* pulse detectability pass should an observing allocation be presented as a plausible detection campaign.

## Reproduction and numerical tripwires

From the repo root:

```sh
python results/venus-psr-j1813-2113/radiometer_feasibility.py
```

Only Python standard library; 5σ flux and SEFD conditions and an extreme narrow-duty sensitivity bound are asserted. The computed numbers depend on explicit assumptions, the cited 2026 hardware spec and an older L-band flux distribution; they are **not** instrument measurements made in 2026 or forecasts of how bright this pulsar will be in 2027.
