# In situ tests of two STEREO comet-tail geometries

A local-velocity reconstruction brings STEREO-A within 366,000 km of the modeled C/2014 C2 (STEREO) plasma trajectory and STEREO-B within 762,000 km for 169P/NEAT. Independent event lists identify a coronal mass ejection at the first spacecraft and a stream interaction region at the second. Helium pickup-ion counts at STEREO-A sit inside the control distribution. These observations leave cometary material unconfirmed in both intervals.

| Source and spacecraft | Closest receipt UTC with 3 h median | Flight time | Heliocentric ecliptic radial, tangential, normal velocity | Miss with 1 h / 3 h / 6 h median |
|---|---|---:|---|---|
| C/2014 C2, STEREO-A | 2014-02-02 07:51 | 1.399 d | 374.8, +27.5, −7.2 km/s | 339 / 366 / 495 thousand km |
| 169P, STEREO-B | 2009-11-15 02:20 | 2.098 d | 313.7, −9.3, −8.3 km/s | 1,000 / 762 / 985 thousand km |

## Geometry

JPL Horizons vectors are heliocentric and use J2000 ecliptic coordinates. Cubic Hermite interpolation uses supplied positions and velocities. Each ten-minute PLASTIC velocity in the six-day search window is smoothed componentwise by a centered median and transformed into the inertial coordinates at receipt. At each receipt sample we minimize the three-dimensional separation from the comet over launch times in the preceding five days. This assumes the local solar wind velocity was constant along the full path. The baseline uses the three-hour median. Longer smoothing moves the answer appreciably. The 72-hour misses are 1.89 million km for A and 1.58 million km for B. The archive timestamp identifies the center used for the geometry scan to the nearest input sample. We use TDB−UTC offsets of 67.184 and 66.184 s respectively, neglecting the millisecond periodic term.

The geometric proximity depends on an extrapolation across structured plasma. Pickup-ion gyro-motion, diffusion, field-line draping, the unknown transverse flow upstream, the unmeasured tail width, and the comet orbit uncertainty all affect an intersection assessment. The values in geometry.csv describe the ballistic model output.

## In situ measurements

| Six-hour median except He+ count rate | C/2014 C2 event | Median of 14 controls | 169P event | Median of 14 controls |
|---|---:|---:|---:|---:|
| Proton density (cm⁻³) | 8.71 | 2.98 | 3.49 | 3.32 |
| Proton temperature (K) | 35,300 | 57,725 | 40,600 | 34,325 |
| Field magnitude (nT) | 9.66 | 5.48 | 3.24 | 2.61 |
| Detected He+ counts per hour, bins centered between 1 and 2 solar-wind speeds | 206.67 | 186.83 | — | — |

The full six-hour STEREO-A window contains 1,240 selected He+ counts. Six of fourteen control intervals have at least the same count rate. The event's proton density exceeds all fourteen controls. This density enhancement occurred during a cataloged STEREO-A ICME lasting 2014-02-02 00:18 to 2014-02-03 12:15 UTC. The event record describes elevated iron charge state, patchy bidirectional electrons and intervals of low proton temperature and plasma beta. The STEREO-B six-hour window straddles the beginning of a cataloged stream interaction region lasting 2009-11-15 00:00 through 2009-11-16 09:00 UTC. The three-hour-median geometric minimum falls inside that region.

[STEREO ICME event list](https://stereo.epss.ucla.edu/media/ICMES_spase_hpevent.txt) · [STEREO SIR event list](https://stereo.epss.ucla.edu/media/SIRs_spase_hpevent.txt) · [Jian et al. STEREO Level 3 event definitions](https://stereo.epss.ucla.edu/l3_events)

He+ counts are instrument detections. The count-rate comparison does not measure helium abundance or identify a comet as the ion source. CDAWeb includes the STEREO-A level-three He+ spectra. We found no matching STEREO-B level-three product and no calibrated C+, N+, O+ or molecular water-group time series for either epoch in the examined CDAWeb catalog or the PLASTIC level-two and level-three directories. Level-one pulse-height products may permit further species analysis with event classification, detector response and priority-based telemetry corrections. The two published solar-wind identifications and the absence of a diagnostic heavy-ion test prevent a positive cometary plasma attribution.

[STEREO-A hourly He+ product](https://cdaweb.gsfc.nasa.gov/hapi/info?id=STA_L3_PLA_HEPLUS_F_VSW_01HR) · [STEREO PLASTIC archives](https://stereo-ssc.nascom.nasa.gov/pub/ins_data/plastic/) · [Drews, interstellar pickup ions with STEREO PLASTIC](https://macau.uni-kiel.de/servlets/MCRFileNodeServlet/dissertation_derivate_00004818/PickupIonen.pdf)

## Reproduction

Run `python analyze.py` with NumPy and SciPy to rewrite geometry.csv from the supplied Horizons ephemerides and PLASTIC vector snapshots. Run `python reproduce_controls.py` with internet access to regenerate control_blocks.csv directly from CDAWeb HAPI. The source queries and hashes are in sources.json. controls_protocol.md records the comparison fixed before examining the composition and field summaries. The screen that proposed these comets selected on geometry, so these local background ranks do not quantify a discovery false-alarm probability.
