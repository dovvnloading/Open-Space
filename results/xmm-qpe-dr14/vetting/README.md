# Disposition of the six HEASARC archive-screen flags

The frozen 0.2–2 keV EPIC-pn finder raised six preliminary flags in 910 source light curves from 637 of the 670 cohort observations ([HEASARC referee audit](../outputs/HEASARC_REFEREE.md)). Five are explained; one remains unconfirmed, low-significance variability:

| XMM obsID / SRCID | Position (deg) | Disposition |
| --- | --- | --- |
| 0851180501, 0864560101, 0904640201, 0904640301 | 195.5005, +27.7828 | **RX J1301.9+2747**, a published QPE source, recovered blind in four separate epochs |
| 0144570101 / 101445701010002 | 206.9016, +17.5680 | Known Seyfert 2 AGN; the flag coincides with a particle flare |
| 0841482101 / 108414821010004 | 205.4356, −0.9756 | **Unconfirmed**: two of three nominal peaks coincide with elevated particle background; the third is about 3σ (below) |

**The 0841482101 flag.** The source is WISEA J134144.52−005832.5 (4XMM J134144.5−005832), a galaxy at SDSS spectroscopic redshift 0.1476, off-axis in a field targeted at the quasar LBQS 1338−0038. The finder reported three soft peaks at MJD(TT) 58880.880, 58881.113 and 58881.498, separated by 5.58 h and 9.25 h ([thread referee note](../outputs/HEASARC_TRIGGER_0841482101_VETTING.md)). An independent check of the merged 0.2–2 keV PPS series, binned to 300 s, with noise taken from the robust scatter of the 279 bins with quiet particle background (median 0.0712, sigma 0.0341 counts/s), finds only two bins above 3.5 sigma, at 20.21 h and 21.54 h. Both fall where the PPS flaring-background rate is elevated (112.3 and 16.8 counts/s against a median of 3.52). The third nominal peak (about 30.7 h, 0.17 counts/s) is about 3 sigma on this noise estimate. So the recurrence rests on one peak during a strong particle flare, one during moderately elevated background, and one near the noise. It is **not** a QPE candidate on this evidence, but it is not excluded either: an event-level soft-band re-extraction with an off-source background region, and the source's eRASS1/Swift history, would settle it. [Soft-band peaks](0841482101_src4_softband_peaks.txt) · [script](soft_band_check.py) · [broad-band table](0841482101_src4_lightcurve.txt) · [script](particle_flare_check.py)

Inputs, fetched from the NASA HEASARC XMM-Newton mirror (`FTP/xmm/data/rev0/0841482101/PPS/`):

| File | SHA-256 |
| --- | --- |
| `P0841482101PNS003SRCTSR8004.FTZ` | `76122b7a8d412dbca9fc72376f85f10a4bfcd03af8cbfc37279058b142e835ed` |
| `P0841482101PNS003FBKTSR0000.FTZ` | `574f5e93c558f54f139b27be5a572bc93107e2c4229faf18732ea78f25e9868b` |
| `P0841482101PNX000SRCTSR8004.FTZ` | `974c99b1b2b42c56b4caafba4905e46dcf2dffa63d95a365f40570228d57bbc4` |

Replay (Python 3.12, astropy 8.0.3, numpy): `python results/xmm-qpe-dr14/vetting/soft_band_check.py P0841482101PNX000SRCTSR8004.FTZ P0841482101PNS003FBKTSR0000.FTZ` and `python results/xmm-qpe-dr14/vetting/particle_flare_check.py P0841482101PNS003SRCTSR8004.FTZ P0841482101PNS003FBKTSR0000.FTZ`.

These checks use PPS time series, not an event-level re-extraction.
