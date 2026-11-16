# Two extreme-softness rejections

| Source and observation | Original imaging events in a 15 arcsec radius | Decision |
|---|---|---|
| 5XMM J150831.7+374025, 0865260201 | Before flag filtering, 3,406 pn events at 0.2–1.0 keV with accepted patterns include 2,667 with FLAG=0x1. Of the 630 events passing FLAG=0 and per-CCD good times, 424 land in a single pn pixel (CCD 3, RAWX 58, RAWY 166) and 549 have PATTERN=2. The three most populated clean pixels share RAWX 58. The neighboring pixel (59,166) contributes 1,299 events before flag filtering. Neither MOS camera has usable exposure at the source. | Reject. The apparent soft source is concentrated at adjacent detector pixels and overwhelmingly in a single split-event orientation. A temperature fitted to this spectrum has no astrophysical interpretation. |
| 5XMM J213004.5−065011, 0780290601 | Of the 143 accepted-pattern MOS1 events at 0.2–1.0 keV, 142 carry FLAG=0x40. A single physical pixel (CCD 5, RAWX 179, RAWY 355) contributes 141 events over 10.14 ks. The strict MOS1 aperture has one clean event versus 0.71 expected background; MOS2 has three versus 0.50; pn has no exposure. | Reject. The detection is dominated by a flagged MOS1 detector pixel. The observation targeted the eclipsing variable EI Aqr at 21h29m20.3s, −06°59′00″, 14.07 arcmin from the catalogue source. |

The event comparisons use standard detector patterns (0–4 for pn, 0–12 for MOS), the imaging exposures only, and the FITS event-list sky WCS. Clean counts use the per-CCD standard good-time intervals and FLAG=0. Unfiltered pixel tallies retain the original flags, with the same energy, radius and pattern cuts. The three measurements for each source record the clean photon distribution, the exposure-scaled local background, and the flagged detector-pixel distribution. The manifest gives the byte-identical input event and exposure products. The clean-photon and background scripts are the same `event_apertures.py` and `exposure_timing.py` used in the pointed-observation audit. The `flagged_pixels.py` command accepts RA (degrees), Dec (degrees) and an input directory containing a single observation's imaging event files.

Target position sources: [ESA observing log](https://xmm-tools.cosmos.esa.int/external/xmm_mission_plan/odf_pps/catstrip_3001_3100.shtml), [Kraków eclipsing-binary ephemerides](https://hpn.as.up.krakow.pl/ephem/AQR.HTM).

Replay the event audit from the mission directory with the combined [archive manifest](analysis/archive_manifest.csv).

~~~sh
python analysis/fetch_manifest.py products1508 --group J150831
python analysis/flagged_pixels.py 227.1321678 37.6736148 products1508
python analysis/fetch_manifest.py products2130 --group J213004
python analysis/flagged_pixels.py 322.5187485 -6.8365591 products2130
~~~

Earlier 4XMM entries for both positions reuse the same exposure. Neither position fell on an active CCD in its respective nearby Chandra pointing. The unmasked Legacy r-band 5σ depth is 23.69 AB at J150831 and 23.90 AB at J213004. The archived spectral fits use the MOS1 5 eV RMF for J213004 and the named pn mode and detector bin with calibration v22.0 for J150831. These instrument-dominated fits are retained as diagnostics. [Chandra footprint receipts](results/boundary-chandra-footprints.json)
