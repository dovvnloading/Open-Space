# Companion-search literature and MCC 124 catalogue checks

MCC 124 is a known compact hierarchical triple. Prior ground-based imaging covers the system. No target-specific published JWST companion catalogue was located in the searches below as of 2026-11-17. The 5835 team has reported its survey in progress. A null result for MCC 124 depends on the coverage, candidate vetting and completeness of the present reduction.

## MCC 124

| Check | Result | Primary record |
|---|---|---|
| Identity | MCC 124, HIP 50156, DK Leo, GJ 2079, 2MASS J10141918+2104297, Gaia DR3 628827128175694976, WDS J10143+2104AB, SBC9 3746 | [SIMBAD](https://simbad.idia.ac.za/simbad/sim-id?Ident=MCC+124) |
| Prior wide search | Subaru/SEEDS reported no companion candidate within 7.5 arcsec in its data. Search sensitivity depends on age, band and separation. | [Brandt et al. 2014](https://doi.org/10.1088/0004-637X/786/1/1) |
| Close component B | The first 2011 image resolved B near 90 mas. VLT/NACO recovered it at 77.3 ± 0.8 mas, PA 338.7 ± 0.2 degrees, with L′ contrast 0.172 ± 0.005. These measurements describe historical orbital phases. The VLT study reports no wider visual candidate for this star. Data are from ESO 097.C-0972(A), OB 2. | [Kammerer et al. 2019, Tables 1 and 3](https://doi.org/10.1093/mnras/stz882) |
| Inner components Aa,Ab | Single-lined 73-day orbit inside the outer visual pair with an estimated period around 3 years. The unresolved reference flux includes all three components. | [Sperauskas, Deveikis & Tokovinin 2019](https://www.ctio.noirlab.edu/~atokovin/papers/Sperauskas2019.pdf) |
| Age | Published evaluations differ concerning membership of the beta Pictoris and Columba groups. Treat age-dependent planet mass limits separately from measured flux limits. | [Messina et al. 2017](https://doi.org/10.1051/0004-6361/201629152) |
| Gaia field | The DR3 cone (30.024 arcsec around the approximate J2000 coordinate) returned the host and one G=20.78 star. Both have two-parameter solutions. The latter lies 19.1 arcsec east and 0.6 arcsec north of the host at J2016.0. Faint stars and galaxies below Gaia's selection threshold remain possible. | [ESA Gaia TAP](https://gea.esac.esa.int/tap-server/tap) |
| Host background track | DR2 gives proper motion (east, north) = (−119.552, −191.194) mas/yr, parallax 42.742 ± 0.189 mas, and astrometric excess noise 0.348 mas. Use parallax at actual JWST epochs and allow for orbital photocentre motion when testing a source. DR3 has no host PM or parallax. | [ESA Gaia TAP](https://gea.esac.esa.int/tap-server/tap) |
| Planet registry | NASA Exoplanet Archive `ps` returned zero entries matching the host HIP, Gaia DR3, TIC or four common host names. This does not exclude unconfirmed objects. | [NASA Exoplanet Archive TAP](https://exoplanetarchive.ipac.caltech.edu/TAP/) |

The programme's STScI phase-II file identifies MCC 124 observations 20 and 59. Successful science files, archive access status, and the time separation must be read from MAST. A repeat label on a proposal is insufficient evidence for independent epochs. [GO 5835 Phase II](https://www.stsci.edu/jwst/phase2-public/5835.pdf), [survey team presentation](https://zenodo.org/records/19686485)

### F444W system flux

AllWISE J101419.08+210427.7 gives combined W2 = 5.980 ± 0.034 Vega mag, quality AAAA, `cc_flags=0000`, and saturated fraction 0.068 in W2. 2MASS J10141918+2104297 gives Ks = 6.261 ± 0.023 Vega mag. Near the bright saturation threshold the formal WISE photometric uncertainty omits brightness-dependent biases. [AllWISE column definitions](https://irsa.ipac.caltech.edu/data/WISE/docs/release/AllWISE/expsup/sec2_1a.html), [WISE bright-source behaviour](https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec6_3c.html)

For a catalogue-scaled contrast denominator, W2-normalized F444W is **0.756–0.760 Jy** across 14 solar-metallicity PHOENIX photospheres (3100–4000 K, log g 4.5 and 5.0). The 2MASS Ks normalization to the same model grid yields **0.639–0.873 Jy**. Using 0.758 Jy gives a multiplicative contrast sensitivity range of 0.87–1.19 for the anchor/model choices explored. The catalogue uncertainty alone on the W2 anchor is about 3.1% in flux. The anchor/model range describes a sensitivity test. The true system flux also depends on the three spectra, variability, passband calibration and WISE bright-source bias.

The calculation integrates the model and CALSPEC Vega in each SVO filter profile, using its energy- or photon-response convention. We add the model F444W minus anchor-band Vega colour to the observed anchor magnitude and convert to Jy using the published SVO F444W zero point (184.1022 Jy). This is broad-band system photometry. It carries no coronagraphic spatial throughput correction. Every positive combination of component spectra has a flux ratio bounded by the minimum and maximum of its components' ratios; the scenario only samples the stated temperatures and gravities. The relative weights and spectra of the individual components were not measured here. [PHOENIX](https://ssb.stsci.edu/cdbs/grid/phoenix/), [CALSPEC](https://ssb.stsci.edu/cdbs/calspec/), [SVO Filter Profile Service](https://svo2.cab.inta-csic.es/theory/fps/)

To reproduce, run `python mcc124_queries.py` for the catalogue responses and `python mcc124_sed.py /tmp/mcc124-sed-cache > mcc124_sed.csv` for the model grid (tested with NumPy 2 and Astropy 8). The catalogue receipts include the exact ADQL, URL, UTC retrieval time, byte count and SHA-256. The photometric script verifies the response curves and FITS spectra against `photometry_inputs.json` before computing the colours.

## Other fields with documented prior searches

| JWST program(s) | Target(s) | Search already reported | Primary record |
|---|---|---|---|
| 4534, 8507 | Altair | Two NIRCam F200W/F444W epochs, 2024-10-11 and 2025-05-26. The appendix lists 11 sources. Ten have two-epoch F444W astrometry, object 11 has low S/N, and object 5 is single-epoch with unconstrained status. All eleven were previously reported. | [Meshkat et al. 2026, Tables 3, 6, 7](https://arxiv.org/abs/2611.13132) |
| 4050, 5835, 6005 | Young moving-group survey targets | Survey team reports active NIRCam companion searches. The conference abstract has no individual candidate catalogue. Target-specific publication status remains open unless documented above or in the inventory. | [Biller et al. 2026](https://zenodo.org/records/19686485) |
| 6122 | Wolf 359, EV Lac, AD Leo | 2024-12-09/10. Published sensitivity analysis covers NIRCam F444W coronagraphy and MIRI F2100W imaging for Wolf 359 and EV Lac. AD Leo served as a MIRI reference. Its original NIRCam observations failed acquisition. | [Bowens-Rubin et al. 2025](https://arxiv.org/abs/2505.15995), [MAST dataset](https://doi.org/10.17909/23r9-q596) |
| 1911, 4857 | WD 1202−232, WD 2105−82 | MIRI F1500W two-epoch observations. Both previously proposed spatially resolved planet-like sources remained stationary. They are likely distant galaxies. The first-epoch imaging dates are 2023-02-09 and 2023-04-21 respectively. Both follow-ups were obtained 2025-06-03. The programme number in the follow-up paper text is transposed relative to the archive. | [Identification](https://arxiv.org/abs/2401.13153), [follow-up](https://arxiv.org/abs/2512.08191) |
| 1184, 3840 | 2MASS J09445422−1220544 (G 161−71) | Original F444W source P1 failed recovery in a deeper second epoch. The follow-up dissertation explicitly dismisses P1 and identifies two extended background sources. | [Bogat et al. 2025](https://arxiv.org/abs/2504.11659), [Bogat 2026, Ch. 3](https://drum.lib.umd.edu/items/55e58cc9-7800-4237-854d-b7bf24157eb6) |
| 1184 | HIP 17695, G 7−34 and other nearby low-mass stars | Nine-target F356W/F444W survey. HIP 17695 P1 was demoted after PSF fitting. An initial G 7−34 signal disappeared after updated cosmic-ray correction. | [Bogat et al. 2025, sections IV–V](https://arxiv.org/abs/2504.11659) |
| 3989, 1668 | HD 92945 | NIRCam F200W/F444W on 2024-05-20 and separate MIRI F1140C work; companion searches published. The MIRI search also covered HD 107146 and HD 206893. | [NIRCam](https://doi.org/10.1051/0004-6361/202555602), [MIRI](https://arxiv.org/abs/2601.02156) |
| 3662 | TWA 7 | 2024-06-21 F1140C observation. The MIRI source and two tentative F444W signatures in subsequent NIRCam analysis, including a counterpart at the MIRI position, have already been reported. | [MIRI](https://doi.org/10.1038/s41586-025-09150-4), [NIRCam](https://arxiv.org/abs/2506.19932) |

### Literature search receipt

Primary-source retrieval on 2026-11-17 used the search strings below. Positive matches are linked above. A targeted search without a hit cannot certify that a target has never been searched.

* `"HIP 50156" JWST planet companion NIRCam 2025 2026 "MCC 124"`
* `"HIP 50156" companion imaged binary limits radial velocity imaging Subaru angular separation Biller MCC 124`
* `"HIP 50156" "triple" radial velocity period visual orbit DK Leo Sperauskas`
* `"5835" "JWST" exoplanet direct imaging Spotlight program Carter targets sub-Jupiter first results 2026`
* `"2M J0944" "Bogat" follow-up NIRCam 2026 P1 common proper motion`
* `site:arxiv.org 2025 2026 JWST MIRI WD 1202-232 2105-820 candidate background second epoch`
* `JWST Altair F444W two epoch companion Meshkat`
* `TWA 7 b JWST confirmation 2026 Crotts TWA 7b NIRCam MIRI common proper motion`
