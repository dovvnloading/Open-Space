# Investigation log

## 2026-11-16 · First checkpoint

- Created a dedicated branch from main at a62807165dd86b1c405751042dbd0a98be00a777 and surveyed earlier research and PRs.
- Downloaded the official [5XMM-DR15 source catalogue](https://heasarc.gsfc.nasa.gov/FTP/xmm/data/catalogues/5XMM_DR15cat_source_v1.0.fits.gz). Its MD5 matches the [SSC release](https://xmmssc.irap.omp.eu/Catalogue/5XMM-DR15/5XMM_DR15.html), 3b58ee2c7c8dbe86009b65e23c35a534. Evaluated 818,656 unique source rows.
- Retained 185 measured-hardness sources and eight with HR2 at its lower boundary. Two later sensitivity checks added 12 very soft and 116 bright, absorption-tolerant rows. Each row belongs to one queue. See the cut flow, individual rows, and executable selector.
- Began independent Gaia DR3, AllWISE, Pan-STARRS and Legacy Surveys cone searches. The first multi-survey blank field with original XMM confirmation is J085618.5+380711. The local r-band five-sigma point-source depth across its error region is 23.81 AB. The catalogue soft X-ray flux divided by the corresponding optical nu*fnu gives only a lower bound of 3.07. Candidate status is under investigation.
- Matched its original observation 0302581801 to PPS source 14 by position and retrieved its EPIC-pn spectrum, background, ARF and RMF, and detector event lists. The apparent ultrasoft excess occurs in both pn and MOS2. Earlier releases list the same source. Pending spectral and novelty checks are tracked in this PR.

This checkpoint records the run so far. No new neutron-star identification is established.

This research has made use of data obtained from the 5XMM serendipitous source catalog compiled by the XMM-Newton Survey Science Center, the XMM2ATHENA project and in collaboration with the XMM-Newton SOC.

### Second checkpoint

- Completed four disjoint queues and the recovery test. The selection retained two of eight known thermal controls. All six missed controls passed the primary hardness cuts and triggered the variability veto. Selection and recovery details are in [SELECTION.md](SELECTION.md).
- Completed SIMBAD screens for all 321 rows. Of 195 within the public eRASS1 hemisphere, 121 have a catalogue entry within 30 arcsec. Each is a positional association for review. The other 126 lack public eRASS1 hemisphere coverage.
- Checked 166 prior positions against all queues. Seven selected sky positions require review against previous work. The bright J114051.9−641848 lead is already a published 4XMM INS candidate.
- Fitted J085618's original pn spectrum and verified the sub-keV excess in independent pn and MOS2 events. Its 54 eV blackbody fit and the alternative thermal-continuum fit differ by less than one unit of Poisson deviance. Deeper optical astrometry and photometry are in progress.

### Third checkpoint

- Reached 338 distinct review positions after adding 22 sources more than 15 arcsec from every DR14 slim source. Five already appeared in the first four queues. The independent Gaia DR3, AllWISE and Legacy catalogue queries returned valid responses at all 338 positions.
- Measured deep CFHT g/r photometry at the native X-ray position and after a local Gaia alignment. The conservative full-region limits are g_AB >25.09 and r_AB >24.84. The resulting catalogue-flux ratio exceeds 7.9. See [optical photometry](OPTICAL.md).
- Fitted a solar-abundance line-emitting plasma through the same pn response as the blackbody. At the measured total Galactic atomic column the blackbody and plasma deviances are 25.66 and 28.30. A second pointed epoch and deeper optical imaging are required for a stronger physical interpretation.
- Documented the pointing-level analysis, hash-verified products and eROSITA upper limits in [the X-ray report](XRAY.md). The [review ledger](results/review-ledger.csv) retains ambiguous neighbors and incomplete coverage statuses.

### Final checkpoint

All five queues yield 338 distinct positions. J085618 survives the event, optical, and prior-position checks as an ambiguous follow-up source. The conservative optical bound is 7.9 in the stated convention. An independent replay confirmed the known-control loss and the instrument-level rejection of two apparent supersoft blank sources. No new isolated neutron star is established.

Final replay regenerated all four baseline source lists and the DR14 differential byte-for-byte from the verified native FITS files. The CFHT photometry receipt was identical. The pn temperature and flagged-pixel counts were reproduced from source products.
