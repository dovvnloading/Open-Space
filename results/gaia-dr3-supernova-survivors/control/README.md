# Control recovery and blend validation

The high-significance, high-transverse-speed search recovers the three Shen objects, 3 of 11 suspected D6 survivors in the reference set. It recovers 0 of 8 LP 40-365-like objects and 0 of 1 runaway helium donor. The parallax S/N gate alone removes eight of the eleven D6 stars. The prototype LP 40-365 has a zero-point-corrected heliocentric transverse speed of 448 km/s, while US 708 moves at just 5.73 mas/yr on the sky. Their statuses as supernova-runaway candidates depend on additional spectroscopy and 3D kinematics.

## Completed comparisons

Twenty unique identifiers and class labels follow the entire Table 2 in El-Badry et al. (2026). This particular convention includes J1637+3631 among suspected D6 stars and J1240+6710 among suspected LP 40-365 stars. The table contains the three Shen objects, the Raddi objects and the hot 2023 discoveries, along with the 2026 additions. The labels describe the published interpretations, some of which are debated. Gaia DR3 has new source ID 688380457508503040 for J0905+2510. [El-Badry et al. 2026](https://doi.org/10.33232/001c.164326) · [Shen et al. 2018](https://arxiv.org/abs/1804.11163) · [Raddi et al. 2019](https://doi.org/10.1093/mnras/stz1618) · [El-Badry et al. 2023](https://arxiv.org/abs/2306.03914) · [Bhat et al. 2026](https://doi.org/10.1051/0004-6361/202558578) · [Hollands et al. 2025](https://arxiv.org/abs/2506.08081)

| Cumulative cut against original 20 | D6 (11) | LP 40-365 (8) | He (1) |
|---|---:|---:|---:|
| Gaia DR3 position, proper motion and parallax, exact ID | 11 | 8 | 1 |
| Proper motion >35 mas/yr and BP−RP <1.2 | 11 | 6 | 0 |
| Plus raw parallax S/N >3 | 3 | 5 | 0 |
| Plus raw heliocentric transverse speed >1,000 km/s | 3 | 0 | 0 |
| Plus uncorrected absolute G >3.5 mag | 3 | 0 | 0 |
| Plus five Gaia quality conditions | 3 | 0 | 0 |

The five conditions require RUWE <1.4, at least 9 visibility periods, no duplicated-source flag, multi-peak fraction ≤2%, and image-parameter determination harmonic amplitude <0.1. All twenty controls pass these five conditions in isolation. A separate sensitivity trial requires a valid individual zero-point calibration, corrected parallax S/N ≥5, corrected transverse speed >1,000 km/s and good astrometric quality. Two D6 objects pass. A separate replay of the exact *published* 2026 selection against these 20 Gaia rows reproduces membership in the authors' Table 4 for all 20, recovering 10 D6 and 3 LP 40-365 objects, with 0 helium donors. The three recovery experiments use correlated Gaia measurements and many controls were discovered with related cuts. Their fractions describe this finite set of known objects. They cannot establish the detection probability, birth rate or spatial density of an unobserved population.

The same literature survey contains 92 classified candidates, 11 explicitly marked blended. The five Gaia quality conditions flag only 4/11 of those blends and 14/81 other entries. The Rybizki fidelity score ≤0.5 flags 3/11 blends and 3/81 other entries. The bright-neighbour statistic norm_dg ≥−3 flags all 11 blends and another 10/81 entries. It flags the genuine control J0927−6335 as well. Thus even a source where the catalogue score rounds to unity can have a close bright neighbour. Missing neighbour values are treated as unflagged. The truth label here records the authors' blend assessment. Unlabelled objects need not all be physically isolated, and these contamination frequencies depend on the literature sample's colour, brightness and astrometric preselection. [Rybizki et al.](https://academic.oup.com/mnras/article/510/2/2597/6460502) · [fidelity and neighbour catalogue](https://dc.g-vo.org/tableinfo/gedr3spur.main) · [prior high-velocity application of the thresholds](https://www.iag.usp.br/sites/default/files/2024-07/LS003_stad1230.pdf)

Only 5/20 controls have Gaia continuous XP spectra (2/11 D6, 3/8 LP 40-365, 0/1 He). The Gaia and Bailer-Jones queries matched all twenty. The transverse-speed diagnostics from the geometric distances multiply the catalogue 16th, 50th and 84th percentiles by measured proper motion. They summarize distance-prior sensitivity with proper motion held fixed. Bailer-Jones distances reflect a Galactic field-star prior. At low S/N this prior can dominate a rare runaway, e.g. US 708 has a geometric median transverse speed of 86 km/s while its published radial velocity is 835 ± 10 km/s. The corrected inverse-parallax point estimates in recovery.csv can be unstable at low S/N and are not velocity posteriors. J0812−5943 has a pseudocolour outside the zero-point calibration domain and its corrected values are empty. Distances and velocities from these methods share Gaia data and must not be multiplied as independent likelihoods. [Bailer-Jones et al.](https://doi.org/10.3847/1538-3881/abd806) · [Lindegren et al.](https://doi.org/10.1051/0004-6361/202039653)

## Reproduction

Use Python, NumPy, Poppler pdftotext and gaiadr3-zeropoint 0.1.0 (wheel SHA-256 39b6ab4c58807a5313674767b494684bd8bea263b3c55a4ea489c8b13b1fe7fb). From this directory:

```sh
python literature.py
python fetch.py
python fetch.py --validation
python recovery.py
```

The literature script verifies the source PDF hash and extracts the 20 Table 2 RVs, RV errors and distance quantiles. These values follow the table used for the authors' kinematic inference, including their revised US 708 RV. It also extracts all 92 Table 4 IDs and verdicts, which can be used as a prior-publication list. The otherwise unmodified query replies have accompanying exact ADQL, source URL, UTC receipt, byte count and SHA-256. The fidelity query contains the 99 distinct IDs in the union of the two tables. Source IDs are strings throughout. Individual distances are in pc in bj.csv and in kpc in published_kinematics.csv. The latter are the authors' runaway-model posterior summaries and include the special spectroscopic-distance constraint for US 708. Gaia parallax is in mas and proper motion in mas/yr.

## Scope

This control calculation does not characterize sources without Gaia colours or astrometry, test discovery-channel independence, or measure selection over luminosity, direction, velocity, extinction, and survey footprint. The control RVs are published inputs for downstream comparisons. Assessments of spectra, escape or trajectories require their own independent data and error propagation.
