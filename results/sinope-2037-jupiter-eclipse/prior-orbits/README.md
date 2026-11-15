# Sinope’s 2037 eclipse in earlier orbit solutions

JPL’s 2013 JUP300 orbit predicts totality from 15:06:43 UTC on October 17 to 00:34:03 UTC on October 18. JUP343 and JUP344 put those contacts within 76 seconds of the JUP347 calculation. Even the 2013 release, evaluated about 24 years after release, differs by less than six minutes at each contact.

| Sinope orbit | Partial ingress, Oct 17 | Total ingress, Oct 17 | Total egress, Oct 18 | Partial egress, Oct 18 | Position difference at Oct 17 20:00 TDB |
|---|---|---|---|---|---:|
| JUP347, reference | 09:19:59 | 15:01:31 | 00:28:19 | 06:00:39 | — |
| JUP300 | 09:25:22 | 15:06:43 | 00:34:03 | 06:06:13 | 846.7 km |
| JUP343 | 09:20:41 | 15:02:43 | 00:27:52 | 06:00:37 | 190.0 km |
| JUP344 | 09:20:46 | 15:02:47 | 00:27:57 | 06:00:42 | 192.6 km |

Each time is UTC at Earth reception under the main calculation’s UTC convention. The comparison holds the Sun, Earth, barycentric Jupiter, light travel and solar/planet limb models fixed. For each archived SPK, it uses J2000-frame states of both Sinope (509) and Jupiter (599) with respect to the Jupiter-system barycenter (5). Subtracting the latter removes that release’s barycenter-to-planet-center displacement. JUP300 stores type 3 Chebyshev records. JUP343 and JUP344 store type 2 records for Sinope. The comparison uses their position coefficients directly. The arbitrary JPL identifiers do not imply statistically independent fits: their astrometry and dynamical assumptions can overlap. The maximum 344-second shift measures historical orbit-solution sensitivity; it does not establish a 2037 timing confidence interval.

Run from this directory with NumPy and SciPy. Point the Horizons cache at the parent search’s state tables. The SPK reader fetches only the needed DAF directory records and coefficient ranges and saves them for reuse. The result file records their offsets, lengths, and SHA-256 digests.

```sh
python verify.py --horizons-cache /tmp/open-space-horizons --spk-cache /tmp/open-space-sinope-spk
```

[JUP300 release notes](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/jup300.cmt) · [JUP343 release notes](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/jup343.cmt) · [JUP344 release notes](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/jup344.cmt) · [NAIF SPK format](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
