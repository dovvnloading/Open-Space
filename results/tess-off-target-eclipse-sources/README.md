# The source of the 2.99024-day dips in TOI-3143.01

Gaia DR3 5853391820057551616 (TIC 384244835), an eclipsing binary 10.435 arcsec from the catalogued position of TOI-3143.01 (TIC 384244834), is the likely source of its 6,460-ppm TESS signal. The Gaia orbital period agrees within 1.17 seconds. Extrapolating the Gaia primary eclipse to the TESS reference epoch 1,824 days later places its midpoint within 4.26 minutes of the TESS dip. Six TESS sectors separately place the missing light closer to the Gaia binary than to the TOI coordinates. Four of these sectors postdate the three sectors cited in the TOI alert.

This is a source attribution for a previously recorded signal. The Gaia object was already classified as an eclipsing binary. ExoFOP still labels the TOI a planet candidate, lists zero time-series follow-up observations, and gives a nominal radius of 140 Earth radii, itself a reason to question the planetary interpretation. The TESS pipeline had already measured the positional offset. The contribution here is its cross-identification with a specified, independently timed Gaia eclipsing binary. The automated catalog and follow-up history retain their published dispositions.

## Measurements

| Quantity | Measurement |
|---|---|
| TOI and Gaia binary positions, ICRS degrees | (216.568723, −63.050473) and (216.566513787, −63.053193121) |
| TESS period from the TOI archive | 2.9902369 ± 0.000017 d |
| Gaia period from the geometric model | 2.99022333 ± 0.00001133 d |
| Gaia primary eclipse | 0.23393 ± 0.00564 mag in G, or a 19.38% loss of light |
| TESS duration | 3.042 h |
| Primary eclipse phase mismatch | 0.000990 cycles |
| Geometry and flux | The TESS registration changes between observations. The projected TOI–binary distance spans 0.49–0.51 pixels in the six informative sectors. |

The next two columns give the fitted center's distance to each proposed location in the detector's own pixel frame. The control values are integrated difference fluxes at orbital phase −0.25 and +0.25. Each event and control uses the same local baseline rule.

| Sector | Eclipses | To TOI (px) | To Gaia binary (px) | Eclipse flux loss (e⁻/s) | Two controls (e⁻/s) |
|---:|---:|---:|---:|---:|---|
| 11 | 8 | 0.430 | 0.100 | 13.40 | +0.86, +0.17 |
| 12 | 10 | 1.754 | 1.595 | 183.02 | −21.10, +1.92 |
| 38 | 8 | 0.539 | 0.043 | 12.64 | −0.68, +0.53 |
| 65 | 8 | 2.448 | 2.794 | −0.06 | +0.08, −0.04 |
| 99 | 6 | 0.446 | 0.097 | 16.33 | +1.12, −0.62 |
| 100 | 9 | 0.569 | 0.063 | 15.66 | +1.15, +1.51 |
| 101 | 8 | 0.523 | 0.021 | 12.59 | +0.04, +1.25 |
| 102 | 8 | 0.510 | 0.142 | 14.60 | +2.52, −0.13 |

Sector 12 has a large off-source residual and a poorly constrained centroid. Sector 65 has no appreciable eclipse deficit at these coordinates. Both remain in the released measurements. Sectors 99–102 were outside the discovery sectors and were evaluated at the original catalog ephemeris. Epoch bootstraps put the directional offset above zero in each of the six informative sectors. The intervals describe eclipse-to-eclipse scatter and omit absolute astrometric and PSF systematic uncertainties.

## Source attribution and prior work

The field contains 47 Gaia entries within 18 arcsec of the binary, including two G ≈ 21 stars within 2.3 arcsec. TESS pixels alone cannot resolve that entire field. The combination of Gaia's angularly resolved primary eclipse, its period and phase, and the independently localized TESS deficit identifies the binary. Summing the 46 catalogued G fluxes with finite photometry within that radius predicts a roughly 0.96% dip in an idealized G-band aperture. TESS measured 0.646% with a different passband and aperture. This is a dilution scale check, not a model of the TESS aperture.

SIMBAD resolves TOI-3143, its two stellar components and the eclipsing binary as separate records. The binary record cites the Gaia EB survey alone. TIC v8.2 resolves the binary as TIC 384244835 at T=15.300. At the archived target magnitude T=12.048, the two-source, equal-throughput approximation with the Gaia eclipse fraction predicts a dip of about 0.92%, before other neighbors or passband effects. The published 4,584-row TESS-EBs catalog and both live Villanova TIC pages yield no entry for either the binary or the target.

The official SPOC difference-image reports for sectors 99–102 already identified the TOI signal as offset, and the sector 99 finding chart labels TIC 384244835 among the nearby stars. The four pipeline centroids fall respectively 1.950, 0.511, 0.317 and 0.732 arcsec from the Gaia EB, while their reported distances from the TIC position are 9.952, 10.027, 9.202 and 10.240 arcsec. These are independent reductions of the same TESS images and supply no independent exposure count. The PDFs do not identify the Gaia eclipsing period and phase as the origin of the signal. Our accepted primary-eclipse windows span BJD 2461047.308–2461148.976 in these four post-alert sectors. Full official XML product URLs and individual UTC catalogue receipts are in [catalog_checks.json](catalog_checks.json).

The TOI is absent from the published event lists in Panahi et al. and the original full Table 2 in Hadjigeorghiou & Armstrong. ExoFOP lists the faint-star QLP discovery and an active-candidate note. The listed 1.9-arcsec companion is a separate location from the binary measured here. Exact TOI, TIC and Gaia-ID searches revealed no earlier paper assigning this event to this source. The search cannot establish absolute publication priority. Controls against prior art mattered: Panahi et al. already identify TOIs 2998, 3258 and 4876, and Hadjigeorghiou & Armstrong give a nearby star a 99.83% source probability for TOI-2541. 

- [TOI archive and column definitions](https://exoplanetarchive.ipac.caltech.edu/docs/API_TOI_columns.html)
- [ExoFOP target](https://exofop.ipac.caltech.edu/tess/target.php?id=384244834) and [public observing notes](https://exofop.ipac.caltech.edu/tess/edit_obsnotes.php?id=384244834)
- [Gaia DR3 eclipsing-binary catalog](https://doi.org/10.1051/0004-6361/202245330) and [archive schema](https://gaia.aip.de/metadata/gaiadr3/vari_eclipsing_binary/)
- [SIMBAD record for the eclipsing binary](https://simbad.cfa.harvard.edu/simbad/sim-id?Ident=Gaia+DR3+5853391820057551616) and [MAST TESS-EBs survey](https://stdatu.stsci.edu/hlsp/tess-ebs)
- [Panahi et al., Gaia–TESS synergy](https://arxiv.org/abs/2209.05845)
- [Hadjigeorghiou & Armstrong](https://doi.org/10.1093/mnras/stad3286) and [source archive containing Table2.csv](https://arxiv.org/src/2310.15833)
- [MAST TESScut](https://mast.stsci.edu/tesscut/) and [all sector footprints for this position](https://mast.stsci.edu/tesscut/api/v0.1/sector?ra=216.568723&dec=-63.050473)

## Reproduce

Of 8,175 TOIs, 8,066 supplied finite periods and coordinates for a cross-match against 2,184,477 Gaia eclipsing binaries, allowing harmonics and separations up to 120 arcsec. The direct-period subset reports all 11 still-PC pairs separated by 3–30 arcsec with fractional period mismatch below 0.002. This was an exploratory source search. Six successful pixel localizations do not constitute a search-wide false-alarm probability. The Gaia field and matching candidate rows are frozen here. Catalogs were retrieved on 2026-11-15 UTC.

Request CSV from the NASA TAP endpoint with the following query. Download the Gaia file from the VizieR mirror. The hashes below refer to the actual unmodified inputs.

```sql
select toi,tid,tfopwg_disp,ra,dec,pl_orbper,pl_tranmid,pl_trandep,st_tmag from toi
```

| Input | SHA-256 |
|---|---|
| [NASA TOI CSV endpoint](https://exoplanetarchive.ipac.caltech.edu/TAP/sync) | bd9c679df97a61f7e0f6250f54d34678038d551dbe3d0db9b67309a38d2ec9c6 |
| [Gaia DR3 fixed-width VizieR data](https://vizier.cfa.harvard.edu/ftp/cats/I/358/veb.dat.gz) | f8269c63fcedb3a35b87f1c6bb2d347d05d3896903d17647204367ca86572a85 |
| [ExoFOP TOI snapshot](https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv) | f0d1c8abf83022dfc4212f29010639e302b907c757ded5acc85708b0aff638d5 |
| Table2.csv in the Hadjigeorghiou & Armstrong source archive | a3b8e411e7aa92159402d62cd61a738aef0f7fd8721693795d096bb9bce0525f |

```sh
python screen.py tois.csv veb.dat.gz > direct-matches.csv
python remeasure.py ./tess-cache --output measurements-rebuilt.json
```

The pixel analysis uses NumPy, SciPy and Astropy. Each sector is a complete 9×9 TESScut FFI stamp. Retain finite, zero-quality-flag cadences. For each predicted eclipse take the median within ±0.37 times the catalogued full duration, then subtract it from the average of two medians at ±0.7–1.4 durations. Each of the three windows must contain at least two observations. Median-combine independent eclipse images, then fit an elliptical Gaussian with a planar local background in a fixed box around the TOI. Transform both catalog positions with each sector's WCS. Resample complete eclipses 200 times to assess positional scatter. The output contains every cycle used, image centroid, controls, 9×9 deficit image, archive product URL and original FITS SHA-256 for each sector.
