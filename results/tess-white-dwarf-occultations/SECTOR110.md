# The Sector 110 aperture dip on TIC 55997020 disappears in a wider aperture

The 2,123 quality-zero 200-second exposures of TIC 55997020 on 2026 November 1–6 produced one isolated trigger and no periodic triggers. At BTJD 4348.75604 the small-aperture flux fell by 15.6% for one exposure, while the image moved and broadened. A fixed 3.5-pixel-radius aperture collects the displaced light. These observations add 5.89 elapsed days on a previously searched WD. The partial visit carries no additional star in the occurrence calculation.

## Raw frames and image check

The prior Sector 110 search took the last 432 frames on camera 1 CCD 1 for TIC 43529091. This extension takes the entire then-listed 2,549-frame roster on camera 1 CCD 4 for the other overlapping clean Sector 109 target, TIC 55997020. All headers were retrieved. The target lay off the stampable area in three frames. There were 2,546 measured stamps and 2,123 quality-zero finite measurements. Each 9 × 9 cutout uses the distortion-corrected per-frame WCS, proper-motion-propagated Gaia position, a radius of 1.5 pixels, and a perimeter median sky. The TIC match lies 0.00050 arcsec from the propagated Gaia coordinate. All 2,123 frames pass the six Gaia-neighbour colour and PRF-width models with minimum WD fraction 94.71% and median nominal fraction 96.48%.

The isolated 4-minute box contains one 200-second sample. Registered out-minus-in images and concentric fixed-aperture measurements compare that frame with its quality-zero flanks under the original detector rules. The field-star comparison uses the brightest other Gaia source in the 3-arcmin cone in the same exposures. The fixed WD aperture loses 20.66% at 1.5 pixels, 4.14% at 3.0 pixels and 0.15% at 3.5 pixels. The simultaneous comparison star, Gaia DR3 3185645657979937152, loses 26.26% at 1.5 pixels. The two centroids shift by (−0.127, −0.148) and (−0.139, −0.154) pixel. The shared shift and recovery of the WD flux in the larger aperture identify an image-motion artifact. The instrumental event remains in the archived raw search and injection background.

[Detector inventory, catalogue responses and query hashes](sector110/tic55997020/inputs.json) · [Per-frame fluxes, centroids, dilution and byte-range hashes](sector110/tic55997020/frames.jsonl.gz) · [Search](sector110/tic55997020/search.json) · [Image checks](sector110/tic55997020/image_vet.json)

## Single-event sensitivity

The isolated-dip search was replayed with injections in the measured flux. Draw 200 common random times uniformly over the elapsed interval from the first good exposure's start to the last good exposure's end, including gaps. Inject a single symmetric trapezoid with ingress and egress each occupying one quarter of the total duration. Integrate analytically over the 200-second exposure and attenuate the modeled flux deficit by each exposure's nominal Gaia target fraction. Re-run the same blind detector. Count a threshold-crossing isolated box within half of the injected duration or one exposure of the midpoint, whichever is longer. A preexisting matched trigger receives no recovery credit. These efficiencies describe the detection step for this WD and window. They depend on the specified depths and profiles. They do not specify the recurrence, duty cycle or geometry of a debris population.

| Intrinsic WD flux loss | Total duration | Recovered / 200 | Recovery |
|---|---:|---:|---:|
| 10% | 4 min | 16 | 8.0% |
| 10% | 16 min | 159 | 79.5% |
| 50% | 4 min | 164 | 82.0% |
| 50% | 16 min | 188 | 94.0% |

[Injection times and outcomes](sector110/tic55997020/single_dip_trials.csv.gz) · [Coverage and binomial Monte Carlo intervals](sector110/tic55997020/single_dip_recovery.json)

Run from the search directory. The archive requests require internet access. The injection replay uses only the deposited epoch series.

```bash
python sector110_followup.py s110-cache sector110/tic55997020
python sector110_image_vet.py s110-cache sector110/tic55997020
python sector110_single_dips.py sector110/tic55997020/frames.jsonl.gz sector110/tic55997020
```
