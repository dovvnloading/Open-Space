# Archival coverage of the 2016 Phoebe eclipse

Saturn's disk crossed the Sun as seen from Phoebe from 2016 May 22 23:39:41 to May 23 10:59:25 UTC as received at Earth. Cassini was closer to Phoebe and would have seen the same change in illumination about 75 minutes earlier. The search below covers May 20 00:00 to May 26 00:00 UTC. The five ALPO-Japan measurements on May 22 overlap the early modeled ring shadow. Their sightlines also pass through Saturn's bright rings on their way from Phoebe to Earth.

## Archive results

| Archive | Search and outcome |
|---|---|
| Cassini ISS (PDS OPUS) | 65 exposures in the six-day window. The last before the disk eclipse ended May 22 00:54 UTC. The next began May 25 23:23 UTC. Over the full query, the closest ISS footprint bounds were 31.18° from Phoebe (Titan image co-iss-n1842532240). Some WAC metadata collapse to a point; even a full wide-angle field cannot bridge these gaps. |
| Other Cassini instruments (OPUS) | 100 VIMS, 42 CIRS and 28 UVIS products. Only two UVIS hydrogen-deuterium absorption-cell products overlapped the disk event in Cassini reception time, May 22 22:24:39 to May 23 09:44:24 UTC. Their solar-wind pointings were 77.32° away. |
| CADC Solar System Object Image Search | JPL target 609 with a 10″ search allowance yielded 32 products at 12 distinct timestamps inside the six-day window. Ten Las Cumbres exposures of Saturn on May 22 from 06:20 to 06:25 UTC produced 30 archive products. The two other matches were ALMA Titan observations on May 20 and 24. None sampled the disk eclipse. |
| SkyMapper DR4 | Checked all 542 image midpoints in the window against a time-dependent Phoebe ephemeris. The 18 exposures during the disk eclipse were at least 62.77° from Phoebe. Minimum center distance anywhere in the six-day window was 6.75°, beyond the SkyMapper mosaic. |
| ESO raw archive | A 0.5° footprint intersection or a 7°-wide pointing box near the track returned 10 records in the six days and zero during the disk eclipse. The spatial preselection uses the target's slow sky motion across this date range. |
| ALPO-Japan | Three observers imaged Saturn during the predicted early ring shadow. Two observers imaged it the following afternoon after disk egress. The dated stacks below contain annotated Saturn imagery without calibrated faint-source photometry. |

### ALPO-Japan on May 22

The projected location is close to the western ring ansa. Each listed line of sight crossed the ring plane before reaching Phoebe, about 10.98 million km beyond Saturn. The table gives center-to-center angular separation and the Saturn-centric ring-plane intersection radius. The direct ring occultation and the ring-scattered planet light both affect any attempted measurement at these epochs.

| Observer | UTC | Separation | Intercept radius | Foreground region |
|---|---|---:|---:|---|
| T. Kozawa | 13:50:24 | 16.94″ | 111,493 km | B ring |
| T. Kozawa | 14:09:06 | 17.25″ | 113,581 km | B ring |
| K. Endo | 14:16:46 | 17.37″ | 114,439 km | B ring |
| T. Horiuchi | 15:25:58 | 18.52″ | 122,198 km | Inner A-ring edge |
| T. Horiuchi | 15:30:11 | 18.59″ | 122,672 km | Inner A ring |

The 15:25:58 position is close to the inner A-ring edge. Topocentric Horizons checks at the approximate cities of Hachiōji (139.32° E, 35.66° N), Koshigaya (139.79° E, 35.89° N), and Uji (135.80° E, 34.88° N) give 111,523, 113,609, 114,465, 122,217 and 122,690 km respectively. A local optical depth cannot be inferred from the nominal boundary.

Kozawa gives fractional UT minutes. Endo reports 1/13-second subexposures. On May 23, Kumamori's image is stamped 16:02:30 UT, and Horiuchi's are stamped 16:12:15 and 16:29:53 UT. The gallery supplies processed images. Without a raw time series, calibration, comparison-source recovery and a model of foreground obscuration, these frames supply neither a photometric detection nor a calibrated upper limit.

[ALPO-Japan May 22 images](https://alpo-j.sakura.ne.jp/kk16/s160522z.htm) · [May 23 images](https://alpo-j.sakura.ne.jp/kk16/s160523z.htm) · [SkyMapper DR4](https://skymapper.anu.edu.au/data-release/dr4/) · [CADC SSOIS](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/ssois/apidoc.html) · [PDS OPUS](https://pds-rings.seti.org/cassini/iss/access.html) · [ESO archive](https://www.eso.org/cms/eso-data/programmatic-access.html)

## Replay

Run `python audit.py --cache /tmp/phoebe-archive` with NumPy to reproduce `coverage.json`. It contains query URLs, exact response hashes, nearest pointings, counts and the individual ALPO ring-plane intersections. The queries select exposures by time before comparing pointings. Instrument pointing bounds are archive metadata. Cassini pointing directions and the geocentric track use JPL Horizons' observer-specific astrometric J2000 RA/Dec at ten-minute intervals. Cassini event times use observer-to-Phoebe light times to find the same emission instant. A spherical sky-footprint bounding rectangle supplies the ISS distance test. The ALPO calculation intersects the Earth–Phoebe line of sight with the plane perpendicular to Saturn's IAU pole at the light-time-corrected position of Saturn. Geocentric and observatory-specific directions differ by much less than a ring width for this geometry.

The result excludes coverage in the indexed fields and times above. It leaves unindexed or privately held images untested. A contemporaneous amateur discussion identifies the apparent approach of Phoebe to Saturn later that week and asks for a finder chart for May 26–27. It provides no measurement of the May 22–23 eclipse. [May 2016 observer discussion](https://www.cloudynights.com/forums/topic/538033-close-pass-of-phoebe-coming-up/)
