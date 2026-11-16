# Mission 03 · A search for thermally emitting isolated neutron stars in 5XMM-DR15

## Investigation log · 2026-11-16

- Created a dedicated branch from main at a62807165dd86b1c405751042dbd0a98be00a777 and surveyed earlier research and PRs.
- Downloaded the official [5XMM-DR15 source catalogue](https://heasarc.gsfc.nasa.gov/FTP/xmm/data/catalogues/5XMM_DR15cat_source_v1.0.fits.gz). Its MD5 matches the [SSC release](https://xmmssc.irap.omp.eu/Catalogue/5XMM-DR15/5XMM_DR15.html), 3b58ee2c7c8dbe86009b65e23c35a534. Evaluated 818,656 unique source rows.
- Retained 185 measured-hardness sources and eight with HR2 at its lower boundary. Two later sensitivity checks added 12 very soft and 116 bright, absorption-tolerant rows. Each row belongs to one queue. See the cut flow, individual rows, and executable selector.
- Began independent Gaia DR3, AllWISE, Pan-STARRS and Legacy Surveys cone searches. The first multi-survey blank field with original XMM confirmation is J085618.5+380711. The local r-band five-sigma point-source depth across its error region is 23.81 AB. The catalogue soft X-ray flux divided by the corresponding optical nu*fnu gives only a lower bound of 3.07. Candidate status is under investigation.
- Matched its original observation 0302581801 to PPS source 14 by position and retrieved its EPIC-pn spectrum, background, ARF and RMF, and detector event lists. The apparent ultrasoft excess occurs in both pn and MOS2. Earlier releases list the same source. Pending spectral and novelty checks are tracked in this PR.

This checkpoint records the run so far. No new neutron-star identification is established.

This research has made use of data obtained from the 5XMM serendipitous source catalog compiled by the XMM-Newton Survey Science Center, the XMM2ATHENA project and in collaboration with the XMM-Newton SOC.
