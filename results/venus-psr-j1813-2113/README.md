# A predicted occultation of PSR J1813−2113 by Venus on 6 February 2027

From the South Pole, three published positions of the same radio pulsar put it behind the solid disk of Venus for **3.4 to 6.0 minutes**, near 16:55 UTC. From Warkworth, New Zealand, the 2024 MeerKAT position gives a 6.1-minute crossing. Two other positions put the pulses outside the limb there. All three predict a miss at Green Bank, West Virginia. These are predictions based on measured radio timing positions. The event has yet to occur.

The difference between the two MeerKAT declinations is 9.83 arcseconds. Venus will have an apparent **radius of 8.99 arcseconds**. An independent subarcsecond position measured before February would discriminate between the two predicted shadow paths. The pulse period is 0.426 s. The MeerKAT release records 21 flux measurements with a median of 0.479 mJy and a range of 0.208–0.570 mJy. Detectability near the planet requires suitable collecting area and a pulsed baseline on either side of the window. [Keith et al. radio measurements and digital release](https://doi.org/10.5281/zenodo.8430591)

## What changes with the pulsar position

Geometric ingress and egress on **6 February 2027 UTC**, rounded to the nearest second. Times describe intersections with the nominal solid radius. Venus's atmosphere, refraction, scattering, unknown proper motion, and telescope sensitivity can change the observed pulse behavior. Antarctic coordinates are geographic reference points. Telescope time and pulsed sensitivity at any site have not been established.

| Site | Keith et al. 2024 | Jayacumaar 2023 | Morris et al. 2002 |
| --- | --- | --- | --- |
| South Pole | 16:52:47–16:58:45 | 16:53:37–16:57:23 | 16:53:49–16:57:11 |
| McMurdo | 16:52:05–16:58:06 | 16:53:06–16:56:33 | 16:53:19–16:56:19 |
| Warkworth | 16:50:11–16:56:19 | Miss, closest 11.56″ | Miss, closest 11.96″ |
| Green Bank | Miss, closest 13.86″ | Miss, closest 23.68″ | Miss, closest 24.08″ |
| Earth center | 16:53:26–16:57:58 | Miss, closest 15.97″ | Miss, closest 16.37″ |

All three nominal South Pole intervals include **16:53:49–16:57:11 UTC**. Moving the two modern coordinates separately over a nine-point grid at their quoted marginal errors retains a South Pole chord in every trial. The much less precise Parkes position allows a miss within its quoted declination error. The grid is a sensitivity check. Astrometric covariances and timing-model systematics preclude a calibrated probability from these numbers. Dates and seconds indicate model output precision. [All eight site calculations and sensitivity tests](site-predictions.json)

| Timing solution | Right ascension, J2000 | Declination, J2000 | Quoted declination error |
| --- | --- | --- | ---: |
| Keith et al. 2024, MJD 59400 | 18h 13m 39.9057931s | −21° 13′ 10.22721″ | 2.816″ |
| Jayacumaar 2023, Table 4.1, MJD 58976 | 18h 13m 39.87s | −21° 13′ 00.4″ | 0.3″ |
| Morris et al. 2002, Table 1 | 18h 13m 39.87s | −21° 13′ 00″ | 7″ |

The two modern fits use the MeerKAT timing program and can share observations. The historical Parkes discovery timing is independent. ATNF adopts the 2024 coordinate. We evaluated the solutions separately with zero proper motion because the catalogue gives no measurement of it for this source. The 2024 data come from `tpa_data/J1813-2113_full.par` and `J1813-2113.fluxtable` in the authors' ZIP (SHA-256 `3c4d3993983299337d9653b728eec13c2354b056aae6deda8434b24e8d7092ef`). [MeerKAT paper](https://doi.org/10.1093/mnras/stae937) · [Jayacumaar, Table 4.1](https://research.manchester.ac.uk/files/295567611/FULL_TEXT.PDF) · [Parkes paper, Table 1](https://adsabs.harvard.edu/pdf/2002MNRAS.335..275M)

## An independent radio-image check

We sampled the two modern timing positions in ASKAP's RACS-low (887.5 MHz) and RACS-mid (1367.5 MHz) continuum images, without fitting peak coordinates. The 2024 position is 2.62 and 2.37 local robust scatters above background. The 2023 position gives 0.57 and 0.48. Translating the same coordinate pair across both images gave an equally large or larger combined contrast at 98 of 5,871 comparison points. Adjacent pixels share the synthesized beam and the two surveys image the same sky. Those counts describe the local image structure. Neither band independently identifies a pulsed source, so both timing solutions remain in the prediction. [RACS-low observations and credit](https://doi.org/10.1017/pasa.2020.41) · [RACS-mid archive](https://casda.csiro.au/hips/RACS/mid/I/properties) · [Measured contrasts](radio-check.json)

## Selection and checks

An exploratory seven-planet search used 4,262 equatorial coordinate pairs from the 4,503 entries in ATNF v2.9.0 between November 2026 and January 2077. J1813−2113 was selected for follow-up because of its near-term nominal Venus crossing. The bounded Venus-only rerun from 15 November 2026 through 1 January 2029 returns just this object inside one Venus radius plus one Earth radius from the geocentric sightline. This rerun uses the fixed catalogue coordinates and includes entries with poor astrometry. Objects without explicit RAJ and DECJ are outside the search. The count describes that screen, not an occurrence rate for pulsar occultations. [ATNF Pulsar Catalogue and acknowledgements](https://www.atnf.csiro.au/research/pulsar/psrcat/) · [Bounded screen result](screen-result.json)

The fixed-coordinate screen is a selection step. Its unique hit is conditional on those catalogue coordinates.

The screen uses Skyfield and JPL DE440s with down-leg light time. The separate site calculation uses Horizons' DE441 astrometric ICRF coordinates, one-minute sampling, airless site elevations and apparent angular diameter. Linear interpolation in the pulsar's tangent plane gives the limb crossings. The direct DE440s Warkworth calculation gave 16:50:11–16:56:19 UTC for the Keith position, agreeing to the printed second. The two reductions share JPL planetary dynamics. Venus is about 44.15° from the Sun. Site altitude alone does not establish telescope access. For example the IAR point is inside the nominal 2024 shadow, but its antennas have a published two-hour hour-angle limit on either side of transit and cannot track this event. [Horizons queries including exact coordinates](horizons-queries.json) · [JPL Horizons documentation](https://ssd.jpl.nasa.gov/horizons/manual.html) · [IAR antenna limits](https://doi.org/10.1051/0004-6361/201936525)

Searches of the primary literature via web search on 14 November 2026 used `"J1813-2113" Venus occultation`, the Unicode-minus spelling, February 2027 with pulsar and Venus, and planetary occultations of radio sources. The searches found the pulsar discovery, timing studies, predictions for geodetic radio sources, and work on using background radio sources to probe planetary environments. They did not locate this event or the astrometry-dependent site predictions. This is bounded evidence for novelty. [Prior radio-source prediction work](https://arxiv.org/abs/0907.4026) · [Planetary radio occultation modeling](https://arxiv.org/abs/1702.07075)

## Replay

From this directory, the Python standard library reproduces every contact from the frozen Horizons rows. The query manifest gives the original URLs and hashes of the returned result text. A fresh response can have a different request timestamp. The eight 31-row tracks freeze the numerical inputs.

```sh
python replay.py > /tmp/venus-site-predictions.json
diff -u site-predictions.json /tmp/venus-site-predictions.json
```

To repeat the bounded catalogue search, install NumPy, SciPy and Skyfield (run with Skyfield 1.55), then fetch the pinned pulsar catalogue and the JPL DE440s binary. The script checks both input digests before examining any source.

```sh
curl -fLO https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.v2.9.0.tar.gz
tar -xzf psrcat_pkg.v2.9.0.tar.gz psrcat_tar/psrcat.db
curl -fLO https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp
python screen.py psrcat_tar/psrcat.db de440s.bsp
```


The independent image check requires NumPy, SciPy and Astropy. The two original CSIRO HiPS tiles share the same sky grid. The script verifies their digests, bilinearly samples the two fixed positions, and translates that pair around a 25–90 pixel annulus to measure the local background and contrast distribution.

```sh
for band in low mid; do
  curl -fL "https://casda.csiro.au/hips/RACS/$band/I/Norder8/Dir460000/Npix465424.fits" -o "$band.fits"
done
python radio.py .
```
