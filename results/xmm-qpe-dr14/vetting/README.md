# Disposition of the six HEASARC archive-screen flags

The frozen 0.2–2 keV EPIC-pn finder raised six preliminary flags in 910 source light curves from 637 of the 670 cohort observations ([HEASARC referee audit](../outputs/HEASARC_REFEREE.md)). All six are explained:

| XMM obsID / SRCID | Position (deg) | Disposition |
| --- | --- | --- |
| 0851180501, 0864560101, 0904640201, 0904640301 | 195.5005, +27.7828 | **RX J1301.9+2747**, a published QPE source, recovered blind in four separate epochs |
| 0144570101 / 101445701010002 | 206.9016, +17.5680 | Known Seyfert 2 AGN; the flag coincides with a particle flare |
| 0841482101 / 108414821010004 | 205.4356, −0.9756 | Particle-background flare (below) |

**The 0841482101 flag.** The source is WISEA J134144.52−005832.5, a galaxy at spectroscopic redshift 0.1476 (NED, queried 2026-11-19), in a field targeted at LBQS 1338−0038. In the PPS 0.2–12 keV source time series binned to 300 s, the only bins above median + 3.5 robust sigma (20.17 h and 21.50 h into the exposure) fall inside the strongest particle-background flares of the observation: the PPS flaring-background rate reaches 112 and 17 counts/s there, against a median of 3.5. Outside the flares the net rate rises slowly from about 0.02 to 0.1 counts/s, as AGN variability does. No recurrent soft eruption survives the particle-flare intervals. [Light curve table](0841482101_src4_lightcurve.txt) · [script](particle_flare_check.py)

Inputs, fetched from the NASA HEASARC XMM-Newton mirror (`FTP/xmm/data/rev0/0841482101/PPS/`):

| File | SHA-256 |
| --- | --- |
| `P0841482101PNS003SRCTSR8004.FTZ` | `76122b7a8d412dbca9fc72376f85f10a4bfcd03af8cbfc37279058b142e835ed` |
| `P0841482101PNS003FBKTSR0000.FTZ` | `574f5e93c558f54f139b27be5a572bc93107e2c4229faf18732ea78f25e9868b` |

Replay: `python results/xmm-qpe-dr14/vetting/particle_flare_check.py P0841482101PNS003SRCTSR8004.FTZ P0841482101PNS003FBKTSR0000.FTZ` (Python 3.12, astropy 8.0.3, numpy).

This check uses the broad-band PPS time series, not the soft-band series the finder used. It shows that the flagged bins coincide with particle flares; it is not an event-level re-extraction.
