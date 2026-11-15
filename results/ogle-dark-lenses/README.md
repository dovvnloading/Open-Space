# OGLE-2026-BLG-0058 is a candidate isolated black hole

A star in the Galactic bulge has been magnified for most of 2025 and 2026 by an unseen foreground object. The public OGLE photometry gives an Einstein-radius crossing time of about 275 days and a microlens parallax of 0.04–0.06. Parallaxes of 0.2 or more are excluded with Δχ² = 213. A long event with a small parallax needs a heavy lens: the fit corresponds to 2.4 or 1.6 solar masses for every milliarcsecond per year of lens–source proper motion, depending on which side of the source the lens passed. The light that is not being magnified is fainter than I ≈ 19.5, which rules out an ordinary star of that mass.

This event ranked first among 5,641 OGLE-IV alerts from 2022–2026 that we fitted for parallax. As of 14 November 2026 we found no public report identifying it as a black-hole candidate. It is a candidate. The mass depends on a proper motion that photometry cannot measure, the two geometric solutions give different answers, and the detection of a non-zero parallax rests on few points after the 2026 peak. [Event page](https://ogle.astrouw.edu.pl/ogle4/ews/2026/blg-0058.html) · [Deposited photometry](data/ews/2026/blg-0058.dat) · [Posterior](posteriors/2026-BLG-0058.json)

![Eight years of flat baseline at I = 18.13, then a slow rise through 2025 to a peak near I = 17.25 in June 2026, with the two parallax fits and the fit without parallax diverging by about 0.05 mag in early 2027.](light-curve.png)

| OGLE-only fit (1,900 points, 2017–2026) | Lens passes with u0 > 0 | Lens passes with u0 < 0 |
| --- | ---: | ---: |
| χ² (no parallax: 1917.7) | 1893.1 | 1893.1 |
| Crossing time tE, days | 277 (+34/−31) | 273 (+34/−30) |
| Impact parameter u0 | 0.40 (+0.06/−0.05) | −0.43 (+0.05/−0.06) |
| Parallax, north | +0.031 ± 0.007 | −0.049 ± 0.011 |
| Parallax, east | +0.022 ± 0.006 | +0.027 ± 0.007 |
| Parallax magnitude (95% upper) | 0.039 ± 0.006 (0.049) | 0.057 ± 0.009 (0.072) |
| Source brightness, I | 18.50 ± 0.22 | 18.38 ± 0.21 |
| Unmagnified share of baseline light | 0.29 (+0.13/−0.17) | 0.21 (+0.13/−0.17) |
| Lens mass per mas/yr of proper motion | 2.4 M☉ | 1.6 M☉ |

The event is at 17:52:38.09, −31:47:36.1 (J2000), Galactic l = 358.3°, b = −2.8°. The peak was near 21 June 2026 (HJD 2461215.7). Values are posterior medians with 16–84% ranges from 40-walker ensemble chains. Uncertainties were rescaled so the best fit has χ²/dof = 1.

## How heavy is the lens

The lens mass is θE/(κ πE), where the angular Einstein radius θE equals the relative proper motion times tE. Photometry measures tE and πE. The proper motion has to come from a model of the Galaxy until someone measures θE directly. We weight each posterior sample by the event rate for a double-exponential disc and a bar, with Gaussian velocity distributions, and draw the proper motion and source distance. The mass function is a power law without an upper cutoff, so the three slopes below bracket the answer and are not a population model. [Prior and kinematics](galaxy.py)

| Mass-function slope | u0 > 0: mass, M☉ | P(M > 3 M☉) | u0 < 0: mass, M☉ | P(M > 3 M☉) |
| --- | ---: | ---: | ---: | ---: |
| −1 (flat in log M) | 18 (+11/−7) | 1.00 | 5.5 (+6.1/−3.0) | 0.77 |
| −1.75 | 15 (+9/−7) | 0.99 | 3.4 (+4.2/−1.8) | 0.56 |
| −2.35 (Salpeter) | 12 (+9/−6) | 0.96 | 2.4 (+2.9/−1.2) | 0.38 |

The two solutions fit equally well. For u0 > 0 the parallax vector points 36° east of north, within 5° of the direction of Galactic rotation, which is what a disc lens in front of a bulge star produces. That solution puts the lens near 4 kpc with θE ≈ 4 mas. For u0 < 0 the motion is across the plane and the prior prefers a slower, lighter and more distant lens. Averaging the two with equal weight gives P(M > 3 M☉) between 0.67 and 0.88.

A main-sequence lens would outshine the measured blend in 81–100% of the weighted samples across every row and total I-band extinctions of 1, 2 and 3 mag. The lens is therefore dark, or it is a star of about one solar mass or less that happens to move unusually slowly. A black hole is the leading explanation for the first case above three solar masses. A neutron star or white dwarf is possible in the low-mass tail of the u0 < 0 solution.

## What the data do and do not establish

**The upper limit on the parallax is robust.** Fixing the parallax magnitude and refitting everything else, with blend flux required to stay above that of a negative I = 20 star, gives Δχ² = 11 at 0.08, 24 at 0.10, 74 at 0.15, 213 at 0.20 and 309 at 0.30. The crossing time stays above 200 days for every parallax up to 0.2.

**The detection of a non-zero parallax is not.** Adding parallax improves χ² by 24.7. Across the 2,328 events in our screen with crossing times of 3–30 days, where parallax should be unmeasurable, 10% reach that improvement and the median is 3.6 against 1.4 expected. Here 15.2 of the 24.7 comes from seven points taken after a 71-day gap that followed 13 July 2026. Without those points the u0 > 0 fit gives tE = 264 ± 42 days and a parallax of (0.019 ± 0.018, 0.021 ± 0.008), still long and still small. Multiplying every uncertainty by 1.6 to match that excess leaves the medians in place, widens the 95% parallax limits to 0.054 and 0.082, and keeps P(M > 3 M☉) at 0.96–1.00 for u0 > 0 and 0.64–0.98 for u0 < 0. Smaller parallaxes mean heavier lenses, so losing the detection does not weaken the case. [Inflated-uncertainty posterior](posteriors/2026-BLG-0058-inflated-errors.json)

**The method recovers a known candidate.** OGLE-2024-BLG-0062 is the target of JWST program 10476, whose abstract describes it as tE ~ 200 days with πE < 0.08. The same code and public data give tE = 231 (+67/−49) or 288 (+134/−86) days and πE = 0.05–0.07 with 95% limits of 0.11 and 0.10. It fell outside the ranked list only because its crossing time is uncertain by more than 25%. OGLE-2026-BLG-0058 has a longer crossing time, a tighter parallax limit and four times as many points. [JWST program 10476](https://www.stsci.edu/jwst-program-info/download/jwst/pdf/10476/) · [Control posterior](posteriors/2024-BLG-0062.json)

**Competing explanations remain.** A slow relative motion below 1.2–1.9 mas/yr brings the mass under three solar masses. Orbital motion of the source about a companion can imitate or cancel parallax. A blend of unrelated faint stars changes the source brightness and crossing time together. Season medians stay between I = 18.115 and 18.129 across 2017–2024, which argues against an intrinsic variable, but one filter cannot test for colour change.

## Frozen forecast

The models separate once the bulge returns to the night sky. Without parallax the source would be at I = 17.885 on 15 February 2027. The parallax solutions predict 17.934 (95% range 17.917–17.950) for u0 > 0 and 17.945 (17.924–17.966) for u0 < 0. On 1 April the values are 17.950 against 17.985 and 17.995. All three are within 0.02 mag of the I = 18.13 baseline by September 2028. OGLE points have uncertainties near 0.03 mag at this brightness, so a few weeks of data decide between a symmetric event and the low-parallax fits. The forecast cannot choose between the two parallax solutions. [Forecast table](data/forecast.csv)

The measurement that settles the mass is astrometric. For θE of 1–4 mas the apparent position of the source shifts by up to 0.3–1.4 mas, peaking roughly 1.4 tE after the photometric peak, in mid to late 2027. That is the method that confirmed OGLE-2011-BLG-0462 as an isolated black hole. [Sahu et al.](https://doi.org/10.3847/1538-4357/ac739e) · [Lam et al.](https://doi.org/10.3847/2041-8213/ac7442)

## The screen

We downloaded the I-band photometry of every OGLE-IV Early Warning System event from the 2022–2026 seasons on 14 November 2026: 5,642 listed, 5,641 with complete table rows and 5,519 with at least 60 usable points. Each was fitted with a point lens, then with annual parallax for both signs of u0. Outliers beyond five scaled sigma were removed once. The model was checked against MulensModel 3.12 on 20 random events, with a largest fractional difference in magnification of 0.0013.

The ranked list keeps fits with 50 < tE < 1500 days, tE known to 25%, both parallax components known to 0.1, peak magnification above 1.34 and a blend fraction above −0.2. That leaves 125 events. They are ordered by the mass implied at 4 mas/yr after raising the parallax by twice its uncertainty. The amplitude and blend cuts were added after a first pass placed low-amplitude brightenings of bright stars near the top. [All fits](data/screen.csv.gz) · [Ranked list](data/ranked.csv) · [Posteriors for the top 20](data/candidates.csv)

| Rank | Event | tE, days | πE | Blend share | Baseline I | Note |
| ---: | --- | ---: | ---: | ---: | ---: | --- |
| 1 | 2026-BLG-0058 | 277 / 273 | 0.039 / 0.057 | 0.29 / 0.21 | 18.13 | This result |
| 2 | 2026-BLG-0587 | 617 / 279 | 0.12 / 0.08 | 0.87 / 0.61 | 17.51 | Still rising; solutions disagree on tE by a factor of two |
| 3 | 2026-BLG-0072 | 301 / 302 | 0.06 / 0.04 | 0.87 / 0.72 | 14.25 | 0.065 mag rise of a bright star over three years; microlensing not established |
| 4 | 2025-BLG-0397 | 325 / 136 | 0.08 / 0.08 | 0.95 / 0.71 | 13.58 | Bright blend, 349 points, parallax improves χ² by 13 |
| 8 | 2023-BLG-0126 | 239 / 226 | 0.06 / 0.08 | 0.04 / −0.02 | 17.94 | Unblended, but 223 points and the peak fell between seasons |
| 11 | 2025-BLG-0038 | 116 / 118 | 0.060 / 0.072 | 0.01 / 0.01 | 15.36 | Parallax improves χ² by 281; 0.65 or 0.55 M☉ per mas/yr |

Pairs are the u0 > 0 and u0 < 0 solutions. Events with a bright blend can have a luminous giant as the lens, which the main-sequence test does not exclude. OGLE-2025-BLG-0038 is the cleanest secondary case: an unblended I = 15.4 source, 4,565 points and a well-measured parallax that puts the lens at two to four solar masses for typical proper motions. Its light curve ends in October 2025, still 0.1 mag above baseline.

## Data rights

KMTNet lists the same event as KMT-2026-BLG-0521. KMTNet's 2026 photometry is proprietary until 1 July 2027 under its data policy. We examined that light curve once as a private sanity check. No KMTNet measurement is used, reproduced or quoted here, and the forecast above was computed from OGLE data alone. The public release next July allows an independent joint fit. [KMTNet data policy](https://kmtnet.kasi.re.kr/ulens/event/2026/)

OGLE photometry is from the Early Warning System of the Optical Gravitational Lensing Experiment. We credit the OGLE team for the observations, reductions and the public alert service. [EWS](https://ogle.astrouw.edu.pl/ogle4/ews/ews.html) · [Udalski, Szymański & Szymański 2015](https://ui.adsabs.harvard.edu/abs/2015AcA....65....1U)

## Replay

From this directory, with NumPy, SciPy, emcee, PyERFA and Matplotlib installed:

```sh
python candidate.py data/ews 2026-BLG-0058 out --ai 2.0 --steps 12000
python candidate.py data/ews 2026-BLG-0058 out-inflated --ai 2.0 --steps 12000 --inflate 1.6
python candidate.py data/ews 2024-BLG-0062 out --ai 2.0 --steps 8000
python forecast.py data/ews 2026-BLG-0058 posteriors/2026-BLG-0058-chains.npz forecast.csv
python plot.py data/ews 2026-BLG-0058 light-curve.png --forecast 500
```

The first command reproduces every number in `posteriors/2026-BLG-0058.json`: the sampler uses fixed random streams. To repeat the full screen, download the archive first. OGLE updates current-season files daily, so later downloads of 2026 events contain more points.

```sh
python fetch.py ews 2022 2023 2024 2025 2026
python screen.py ews screen.csv
python rank.py screen.csv ranked.csv
```

Times are HJD − 2450000 in the code and tables. Magnitudes are OGLE I. The parallax vector is geocentric with the reference time fixed at the best no-parallax peak, and its sign conventions match MulensModel. File hashes and retrieval details are in [sources.json](sources.json).
