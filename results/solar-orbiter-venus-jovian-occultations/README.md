# Venus blocks Earth and Jupiter in succession during Solar Orbiter’s December 2026 flyby

On 24 December 2026, Venus’s solid limb is predicted to cover the direction of Earth’s center for 139 seconds and Jupiter’s center for 888 seconds as seen from Solar Orbiter. The two intervals overlap for 18.65 seconds. Earth then clears the limb while Jupiter stays behind Venus for another 14 minutes 29 seconds. The Sun’s center stays more than 41° clear of the limb throughout the 40-minute window tested. The distinction matters for radio measurements during the flyby. A sufficiently bright signal that persists when its compact source is geometrically hidden requires another contribution or propagation path.

| Spacecraft reception on 24 December (UTC) | Solid-body geometry |
| --- | --- |
| 22:55:02 to 22:57:03 | Earth’s center is behind Venus. Jupiter’s center is clear. |
| 22:57:03 to 22:57:22 | Both centers are behind Venus. |
| 22:57:22 to 23:11:51 | Jupiter’s center is behind Venus. Earth’s center is clear. |

## Related observing windows

The same calculation identifies these comparison windows. The durations refer to the center of Jupiter and a 6,051.8 km spherical Venus. The selection began with known Venus flybys and was exploratory. The table contains selected encounters only.

| Spacecraft | UTC date | Jupiter ingress–egress (spacecraft UTC) | Duration | Other sightlines |
| --- | --- | --- | ---: | --- |
| Solar Orbiter | 2021-08-09 | 04:27:42–04:43:24 | 942 s | Earth was also behind Venus until 04:38:02, giving 620 s of simultaneous coverage. |
| Parker Solar Probe | 2023-08-21 | 12:48:36–14:29:25 | 6,049 s | Earth and the Sun were clear throughout this window. |
| Solar Orbiter | 2025-02-18 | 20:48:07–20:50:00 | 113 s | Earth's center stayed at least 1.68° clear. Local Venusian plasma emission is a particularly important alternative near this pass. |
| Solar Orbiter | 2026-12-24 | 22:57:03–23:11:51 | 888 s | Earth first, 18.65 s of overlap. |
| Solar Orbiter | 2028-03-18 | 08:15:19–08:26:59 | 700 s | Earth and the Sun remain clear in the tested window. |

ESA’s mission summary gives the Venus flyby dates. NASA published the 2021 optical stellar occultations. Published analysis of the February 2025 Solar Orbiter flyby uses high-cadence RPW spectra to measure local electron properties. Galileo observations provide a precedent for locating Jovian radio emission by occultation behind Ganymede. Searches pairing the spacecraft names and flyby dates with “Jupiter”, “Venus” and “radio occultation” identified these related measurements. A priority claim about the present dates would require a wider literature and mission-planning review. [ESA flyby schedule](https://www.esa.int/Science_Exploration/Space_Science/Solar_Orbiter/Solar_Orbiter_perihelia_and_flybys) · [NASA 2021 imaging](https://www.nasa.gov/image-article/during-close-pass-solar-orbiter-captures-venus-glare/) · [2025 flyby RPW analysis](https://doi.org/10.1051/0004-6361/202557868) · [Galileo radio occultation](https://ntrs.nasa.gov/citations/19970025141)

## Archived radio data at the Jupiter contacts

The source test uses the two solid-limb contact times above, independently of the radio intensity. For each frequency channel we compared the median power immediately outside and inside each contact. A compact, steadily emitting source behind Venus predicts a positive outside-minus-inside contrast at **both** contacts. We then took the median contrast across channels. The 2021 and 2023 tests use two minutes on each side of each contact. The shorter 2025 event uses 20 s. These windows test step-like changes. They cannot identify a source from a naturally variable background.

| Encounter and antenna | Band (MHz) | Ingress (dB) | Egress (dB) |
| --- | --- | ---: | ---: |
| Solar Orbiter, 2021 | 4–8 | +0.33 | −0.09 |
| Parker, 2023, V1–V2 | 4–8 / 8–16 | −0.03 / +0.05 | −0.16 / −0.09 |
| Parker, 2023, V3–V4 | 4–8 / 8–16 | +0.05 / +0.10 | −0.18 / −0.06 |
| Solar Orbiter, 2025 | 4–8 / 8–16 | +0.25 / +0.36 | +0.66 / +0.24 |

At the 2021 and 2023 Jupiter contacts the egress medians have the opposite sign to an emerging steady source. The 2021 upper band fails the coverage requirement. In 2025 both simple contrasts are positive. Only 0/397 and 2/396 eligible time-shifted pairs have at least as strong a smaller-limb score in the lower and upper bands. These are descriptive comparisons within the same observing day. Near closest approach the spacecraft crossed the measured Venusian ionosphere, and a changing local signal can also produce a positive contrast in brackets around a dip.

A local trend test distinguishes a declining-then-rising background from abrupt changes at the two contacts. For each channel we separately fit straight lines to the 20 s on each side, extrapolate to the contact, then take the median outside-minus-inside difference over channels. The fitted ingress/egress jumps are **−0.10/+0.28 dB at 4–8 MHz** and **−0.06/−0.15 dB at 8–16 MHz**. The 2025 data thus do not establish sharp disappearance and recovery of a Jovian source. The lower-band simple contrast also strengthens when both predicted contacts are moved 10 or 20 s earlier. Full results include fits with 5, 10, 20 and 30 s on each side and offsets from −60 to +60 s. [Published Venus-flyby plasma observations](https://doi.org/10.1051/0004-6361/202557868) · [Contact contrasts](radio_contacts.csv) · [Local sensitivity tests](radio_sensitivity.csv)

Solar Orbiter uses positive calibrated PSD_SFU, measured samples (INTERPOL_FLAG=1), publication quality (QUALITY_FLAG=3), a fixed sensor configuration and each channel’s sweep-time offset. Samples below the product background are censored by this selection. The results describe measurable emission above that background. Parker uses positive auto-spectral power, fixed high gain, zero overrange flags and the 40-spectrum accumulation used at the contacts. Each band needs at least three eligible channels, each with 6 (2021), 10 (2023) or 10 (2025) samples in each of the four windows. The same-day controls shift both contacts by equal multiples of four side-window lengths, exclude the real event and buffers, and require the same number of eligible channels. Control epochs can overlap one another. Neither the two bands nor the two dipoles are independent replications. [Radio product manifest](radio_inputs.json) · [All eligible shifted controls](radio_controls.csv)

## Tests and scope

The primary calculation interpolates airless apparent directions from JPL Horizons in Cartesian coordinates and solves for zero angular clearance between the background direction and Venus’s apparent limb. The frozen tables sample Venus every 10 s and the three background directions every minute. Each distance check places Venus in the foreground. For each selected window the replay searches a one-second grid for both contacts, refines each root, and checks the ends of the window. A negative minimum clearance in the CSV means that the direction crosses the modeled disk. Times use spacecraft reception UTC. [JPL Horizons manual](https://ssd.jpl.nasa.gov/horizons/manual.html)

A separate CSPICE calculation uses converged light time, stellar aberration, J2000 vectors, DE421 and successive ESA orbit kernels from 14 October and 4 November 2026. It uses the Jupiter system barycenter as a directional check on the Horizons Jupiter-center solution. At the December flyby the November calculation agrees with the four Horizons contact labels to about 0.002 s. The two ESA kernel versions differ by 0.031 s at most in December and by 0.495 s in March 2028. These are implementation and orbit-revision checks with shared spacecraft navigation data. Neither value measures the absolute timing uncertainty. Future trajectory updates and UTC changes can shift the prediction. [ESA kernel dataset](https://doi.org/10.5270/esa-kt1577e) · [STIX kernel mirror](https://datacenter.stix.i4ds.net/pub/spice/latest/kernels/spk/)

At December’s Jupiter contacts a modeled source anywhere inside a projected disk of radius three Jupiter equatorial radii remains entirely hidden for 887.23 s. Its earliest possible ingress and latest possible egress span 888.16 s. Offsets of −3 and +3 in the CSV denote full and partial coverage of that envelope. This is a sensitivity calculation for compact emitters. Increasing Venus’s blocking radius by 100 or 300 km increases the center occultation to 902.87 or 933.31 s. Actual refraction, diffraction, the ionosphere, source extent, solar or local plasma emission, antenna response and intermittent Jovian emission all affect a measured light curve. The Earth-center calculation does not predict a ground-station-specific loss of lock or model a terrestrial radio source.

The published RPW high-frequency receiver range extends to 16 MHz. To test a Jovian interpretation, an identifiable emission feature would need to be present on both sides of the occultation, with a disappearance and reappearance tied to both predicted limbs after accounting for cadence, frequency and instrumental state. An independent monitor can check whether the source continued radiating. In December the Earth and spacecraft view directions differ by 7.34° at Jupiter, and an identical wavefront would reach Earth 33.91 s later. Narrow beaming and rotation can cause much larger differences in burst timing or visibility. Instrument operations, radio brightness and detection during the future encounter remain open. [RPW instrument and 2025 observations](https://doi.org/10.1051/0004-6361/202557868) · [Cassini and Earth observations of directional emission](https://doi.org/10.1029/2003JA010260)

## Reproduction

From this directory, with Python 3, NumPy and SciPy installed, run:

```sh
python reproduce.py
```

This rebuilds [contacts.csv](contacts.csv) and [inputs.csv](inputs.csv) offline using the 20 compressed Horizons result texts in `raw/`. The manifest gives each exact query, epoch count and SHA-256 of the uncompressed UTF-8 result. `python reproduce.py --refresh` explicitly replaces the frozen tables with current trajectories.

With SpiceyPy also installed, run the orbit-version check using a cache outside the checkout:

```sh
python check_kernels.py --cache /tmp/solo-venus-ephemerides
```

The second command retrieves four pinned public kernels (approximately 54 MB), verifies their SHA-256 hashes and rebuilds [kernel_contacts.csv](kernel_contacts.csv) and [observer_geometry.csv](observer_geometry.csv). Blank contact times mean that the direction stayed clear in the selected interval.

With cdflib installed, the radio comparisons use three SHA-256-pinned mission CDFs (about 47 MB). The following command downloads the missing inputs into the specified directory and rebuilds all three radio CSVs.

```sh
python check_radio.py --cache /tmp/venus-contact-cdf --download
```
