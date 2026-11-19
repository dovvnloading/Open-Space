# eRASS1 DR1 top-35 follow-up — six known neutron stars, no new XDINS validated

**Result, 19 November 2026 UTC (bounded, incomplete):** Of the **35** top-60 soft eRASS1 sources without a published high-probability Gaia/Legacy/CatWISE association in [PR #204](../erass1-xdins-search/README.md), **six are known neutron stars**: the four Magnificent Seven XDINS **RX J0720.4−3125, RX J0806.4−4123, RX J1308.6+2127 and RX J1856.5−3754**, plus rotation-powered pulsars **PSR J0659+1414 (B0656+14)** and **PSR J1057−5226 (B1055−52)**. The pulsars are independently present in a downloaded, hash-receipted ATNF release and in [SIMBAD B0656+14](https://simbad.harvard.edu/simbad/sim-id?Ident=PSR+B0656%2B14) and [SIMBAD B1055−52](https://simbad.harvard.edu/simbad/sim-id?Ident=PSR+B1055-52). **The remaining 29 are unclassified, not proven AGN, stars, clusters, spurious detections, or newly identified INS. No newly discovered XDINS passes the evidence bar.** This is **not** a calibrated null or a claim that the 29 contain no neutron stars.

## What was checked

| Gate | Measured evidence | Limitation |
|---|---|---|
| Frozen 35-source census | All 35 distinct eROSITA UIDs, plus exact input SHA-256 and per-source classifications | No secure catalogue association is **not** an optical nondetection |
| Astrometry and morphology | eROSITA POS_ERR **0.944–4.447 arcsec**; catalogue EXT=0, EXT_LIKE≤6 by prior fixed selection | Extent and positions are catalogue-model values; no event-image spurious validation |
| Image-level optical photometry | **10/35** have interpretable positional optical g/r image apertures from Legacy DR10 or PS1; **7** have provisional g five-sigma limits | PSF blending, stacks, correlated noise, extinction and band calibration are not fully controlled; many LS cutouts returned HTTP errors; PS1 outside its footprint for 31; DES SIA returned no usable image |
| Native infrared images | **35/35** unWISE coadd cutouts downloaded; W1/W2 fixed-position aperture measurements and background/noise estimates saved | All positions are flagged potentially crowded; IR limits alone do not give fX/fopt |
| Optical flux-ratio test | **0/7** provisional, *observed-band* f_X/(ν_g f_ν,g) limits exceed **10⁴**; highest **1480** at 1eRASS J062619.8-694546 (g≈22.893) | eROSITA flux assumes a catalogue spectral model; these are provisional aperture bounds, **not independently certified physical f_X/fopt limits** |
| Archival X-ray counterparts | 4XMM-DR14/HEASARC **20**, Chandra CSC **8**, Swift 2SXPS **7** source positions with at least one nearby catalogue record | Positional matches are not confirmed same-source spectra, flux-history measurements or long-term constancy; absent match is not unobserved |
| Galactic extinction and distance | IRSA Schlafly–Finkbeiner line-of-sight E(B−V) retrieved for **35/35**, range **0.0241–3.4908 mag**; N_H proxy recorded; observed-flux L_X at hypothetical **0.3, 1 and 50 kpc** for each | Not fitted N_H, distance, unabsorbed luminosity or cluster/LMC membership |
| Neutron-star spectrum and time checks | Softness/P1–P3 counts from the frozen eRASS1 screen; known-pulsar crossmatch | **No** response-folded kT=40–120 eV blackbody vs AGN power law/coronal model likelihoods; **no** eROSITA/archival event-level variability or pulsation searches |
| Prior work and novelty | ATNF catalogue downloaded (4,262 parseable pulsars); named known objects checked; [Kurpas et al. 2024](https://doi.org/10.1051/0004-6361/202450248) reported 33 existing candidates | Kurpas full positional table and 2025–26 follow-up novelty audit **not retrieved/complete**; cannot verify that any unresolved source is novel |

**Best unassociated-source optical ratio is far below the 10⁴ bar**, and even that ratio is provisional. A source without a direct optical image has **no optical upper limit assigned**. A source with a weak association or an XMM catalogue match is **not** automatically an AGN/star/cluster. The catalogue eROSITA softness proxy is not a fitted thermal spectrum; the dust-derived N_H is only a line-of-sight plausibility aid. The 50-kpc example is a **hypothetical luminosity check**, not a distance measurement.

## Reproducibility and status

The [single Actions job](../../.github/workflows/erass1-xdins-35-followup.yml) executes [run.py](run.py) on pushes to its research branch **only when scripts or requirements change**, downloads survey images and public astronomical catalogues, computes the full 35-source evidence, and batch-commits `output/**` without recursion. It has no pull-request-wide trigger. Input hashes and exact source names remain from PR #204; original upstream failures and successful receipts are preserved. Intermediate source classifications and missing evidence are in:

- [All 35 classifications, positional errors, optical ratios, archive statuses](output/classifications.csv) and [full per-source evidence including luminosity and absorption proxies](output/source_evidence.json).
- [Optical/IR FITS and dust receipts](output/image_receipts.json), [XMM/Chandra/Swift positional matches](output/archive_receipts.json), [ATNF pulsar matches](output/atnf_receipts.json), [failed Kurpas table receipt](output/literature_receipts.json) and [machine summary](output/summary.json).
- [Independent evidence-level referee](REFEREE.md) and [machine-generated limit summary](output/REFEREE.md).

**Referee outcome:** acceptable as an incomplete, auditable follow-up of catalogue-selected objects, **not** a confirmed new XDINS or calibrated absence. Under [docs/REMARKABLE_FINDINGS.md](../../docs/REMARKABLE_FINDINGS.md), no **[REMARKABLE]** label or log entry is warranted.
