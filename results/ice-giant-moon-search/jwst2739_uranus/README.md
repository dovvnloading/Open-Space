# Uranus in February and September 2023

**No new satellite passed the two-visit search between Perdita and Puck.** For observable circular, equatorial orbits from 79,500 to 80,500 km, all 18 off-grid sources planted at **0.50 μJy (AB 24.65)** were recovered in both visits in each filter. The two observations were 211 days apart. At 0.10 μJy, recovery fell to 1/18 in F140M and 7/18 in F300M. [F140M injection record](f140m_plants.json) · [F300M injection record](f300m_plants.json)

**Perdita was detected in the original F140M images on both dates.** Relative to nominal JPL URA184, the median individual-exposure longitude differences are **−1.44° on February 6** and **−1.03° on September 4**. The medians use four and six native detections, respectively. A comparison moon, Belinda, gave median offsets of −0.31° and −0.01° in the same filter. For Perdita the median nominal sky separations were 4.21 and 3.28 short-wave pixels. The F140M fits have median formal amplitude-to-error ratios of 11.3 and 12.2. The table records each position, nearest other known moon, fitted flux significance and internal positional error. [Native Perdita RA and Dec](f140m_perdita_positions.csv) · [Belinda astrometric control](f140m_belinda_positions.csv)

| 2023 result | F140M | F300M |
|---|---:|---:|
| Feb / Sep native exposures with fitted planetary limb | 4 / 8 | 4 / 10 |
| Tested two-visit orbits | 960,480 | 960,480 |
| Passing the midpoint coverage check in both visits | 953,924 | 828,729 |
| Positive / inverted templates passing the initial cut | 0 / 0 | 0 / 0 |
| Random orbits with usable trailed-PSF support in both visits | 119/120 | 105/120 |
| Recovered at 0.10 / 0.50 / 2.0 μJy, among 18 per flux | 1 / 18 / 18 | 7 / 18 / 18 |
| Belinda at nominal position, February / September score | 81.2 / 57.2 | 14.9 / 32.3 |

These are visits 02739003001 and 02739011001. The long-wave NIRCam detector contained the complete central disk in all 14 F300M exposures. The short-wave analysis used the 12 exposures with fully enclosed, measurable limb slopes. The short-wave detector gap and clipped detector edges reduce some orbital coverage. The exposures in a given visit are dithers, and the two filters were observed simultaneously. The bands and dithers do not increase the number of independent dates. [F140M exposure inputs and hashes](f140m_inputs.json) · [F300M exposure inputs and hashes](f300m_inputs.json) · [JPL state and source checksums](uranus_sources.csv)

## Detection and completeness

The search models equatorial circular motion with Uranus's J2 and J4 terms over the two spacecraft light-time-corrected epochs. The grid uses radius steps of 0.5 km and orbital longitude steps of 0.75°. It masks projected rings, the planet and exposure-integrated trails of all 29 known moons, then subtracts an azimuthal ring profile and a broad scattered-light model. The matched filter uses the STScI encircled-energy table for each band and allows DQ values 0 and 4. The two-stage screen first samples PSF-filtered images at exposure midpoint with a threshold of 80% of the final rank. For a track reaching that screen, the full exposure trails must reach scores of at least 5 on both dates and 7 overall, retain at least 2 after leaving out any one exposure per visit, and have PSF support in two or more exposures per visit. No science track reached the initial threshold in either sign. [F140M grid](f140m_gap.json) · [F300M grid](f300m_gap.json) · [Orbital and score implementation](screen.py)

The injection study draws 120 independent uniform radius/phase combinations per band and selects 18 observable orbits without replacement. It adds individually trailed model sources to the calibrated images before fitting the ring profile and broad glare. Every recovered track passes the same midpoint gate and full trail/visit cuts as a real source, and lies within 1.5 pixels of truth in every exposure. Only the subset of the blind grid able to satisfy that spatial association requires retesting for each injected moon, since the science search emits every passing trial. At 0.50 μJy, the exact two-sided 95% interval for 18 successes out of 18 is **81.5–100%** per band. The study fixes the measured noise, DQ distribution, PSF and orbital plane. It does not measure native-integration completeness, other radii, the full sky, or unseen orbital phases. [F140M individual plants](f140m_plants.csv) · [F300M individual plants](f300m_plants.csv)

At 0.50 μJy, representative diameters at band geometric albedos of 0.04, 0.08 and 0.16 are **27.0, 19.1 and 13.5 km in F140M**, and **37.6, 26.6 and 18.8 km in F300M**, using a unit phase function, the CALSPEC Sun, SVO throughput and the mean geometry of the program. Real icy satellites can have different albedos in the two passbands. [F140M size inputs](f140m_size.json) · [F300M size inputs](f300m_size.json)

## Known-moon checks

Separate reductions held out one known moon at a time. For Belinda the nominal track passed in both filters and none of 36 displaced comparison tracks did. Puck had a supported F140M measurement in February only. At Perdita's nominal coordinates, F140M scores were approximately 0.0 and 10.4 for February and September. A sweep across all 360° in 0.25° increments recovered the strongest F140M tracks at offsets of −1.50° (18.58) and −1.00° (24.59). Every track outside ±10° scored below 3.55 on each date. The weaker F300M September data independently peaked at −0.75° (6.50), while in February the F300M score peaked at 3.09. No F300M Perdita exposure met the individual native-centroid threshold. The native F140M positions above come from fitting a trailed Gaussian PSF and a local sky plane, with the moving frame tied to the fitted planetary limb and Horizons. Their formal errors exclude WCS, PSF shape and limb systematics. [F140M longitude search](f140m_perdita.csv) · [F300M longitude search](f300m_perdita.csv) · [F140M Belinda control](f140m_control714.json) · [F300M Belinda control](f300m_control714.json)

## Reproduction

From the parent result directory, retrieve program 2739 calibrated science windows with the following commands. Source range receipts include each original filename, FITS offsets, lengths and SHA-256 values. The cutout ledgers also identify disjoint chips. [F140M source ranges](f140m_ranges.csv) · [F300M source ranges](f300m_ranges.csv) · [F140M cutouts](f140m_cutouts.csv) · [F300M cutouts](f300m_cutouts.csv)

```sh
S=jwst2739_uranus
python jwst_cutouts.py manifests/products.csv cache/f140 --planet Uranus --proposal 2739 --filters F140M --half 205
python jwst_cutouts.py manifests/products.csv cache/f300 --planet Uranus --proposal 2739 --filters F300M --half 175
python moon_spk.py 2459981.5 2460192.6 cache/ura184-u2739.bsp cache/kernel_ranges
python ephemeris.py Uranus JWST 59981.225 59981.245 cache/horizons cache/ephem.npz --window 60191.918 60192.015 --ura184 cache/ura184-u2739.bsp
python prepare.py cache/ephem.npz cache/f300.npz cache/f300/*.fits --centers "$S/f300m_centers.csv" --size 255 --scale .063 --samples 9 --psf-sigma .68 --mask-pad 5 --allow-corrected-jumps
python "$S/screen.py" cache/f300.npz cache/Encircled_Energy_LW_ETCv2.txt output_f300_grid
python "$S/plants.py" cache/f300.npz cache/Encircled_Energy_LW_ETCv2.txt output_f300_plants --flux 1e-7 5e-7 2e-6 --n 18
```

For F140M, use the 12 rows in `f140m_centers.csv`, size 401, pixel scale 0.031″ and initial PSF sigma 0.78. The preparation stage derives the trial geometry while the matched-filter stage uses the measured radial profile table. To recover the supplied limb centers, use `limb_center.py` at initial radius 28 long-wave pixels or 58 short-wave pixels, with `--use-finite-partly-saturated-limb`. To unmask Perdita, run the same preparation with `--control-body 725`, then apply `perdita.py` and `native_positions.py`. Belinda uses ID 714 and Puck uses 715. The longitude sweep is meaningful for a known body since it retains the catalogued orbital velocity and radius. [STScI short-wave EE table](https://jwst-docs.stsci.edu/files/216457060/216457083/1/1762453590926/Encircled_Energy_SW_ETCv2.txt) · [STScI long-wave EE table](https://jwst-docs.stsci.edu/files/216457060/216457084/1/1762453591836/Encircled_Energy_LW_ETCv2.txt)
