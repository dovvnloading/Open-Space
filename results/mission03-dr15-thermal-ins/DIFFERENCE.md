# Objects absent from the 4XMM-DR14 slim source catalogue

This supplemental search changes the 5XMM intensity and variability criteria and explicitly tests positional novelty against the previous single-observation source catalogue. The DR14 input contains 692,109 unique source positions built from 13,864 observations. DR15 contains 818,656 source positions and includes 14,616 observations. A missing DR14 neighbour can arise from stacking, different source fits, recent sky coverage, motion, diffuse structure or instrumental effects.

Calculate a second complete cut flow on the native DR15 source table. Require SUM_FLAG ∈ {0,1}, EXTENT = 0, DIST_NN ≥20 arcsec, STACK_DET_ML ≥30, EP_CTS ≥70, 0 < RADEC_ERR ≤3 arcsec and the sum of finite nonnegative 0.2–2 keV band fluxes ≥3×10⁻¹⁴ erg cm⁻² s⁻¹. Allow any HR1. Require interior HR2 with positive uncertainty and HR2 + σ₂ < –0.5, or HR2 = –1 as a censored case. Require at least two EPIC cameras each with ≥3σ in either band 1 or band 2. In both hard bands, require finite nonnegative flux and a finite positive error with flux <3σ. Reject VAR_FLAG = true. For VAR_PROB <10⁻⁵ retain only sources with a measured 0 < FRATIO < 2, allowing precise measurements of modest secular changes. Untested inter-epoch variability passes the screen.

Map all DR14 SC_RA/SC_DEC positions to unit vectors and use a Euclidean nearest-neighbor search in three dimensions. Convert chord lengths to angular separations exactly. Require nearest DR14 source >15 arcsec, with no cuts on DR14 quality. The result includes the nearest DR14 SRCID and separation. Rank the survivors by the same multiwavelength probability and observed soft-flux rule as the other queues. The 15 arcsec cutoff can join unrelated sources in dense regions or miss displaced repeat detections, including high proper-motion objects. The test does not address a prior appearance in 4XMM-DR14s, another mission, a published candidate study, or a source known outside X-rays. Check the 4XMM stacked catalogue as well for positionally novel queue members.

```sh
curl -fL -o 4XMM_DR14cat_slim_v1.0.fits.gz https://heasarc.gsfc.nasa.gov/FTP/xmm/data/catalogues/4XMM_DR14cat_slim_v1.0.fits.gz
echo '1f3c9faf54ebf790139be62f1494c1f0  4XMM_DR14cat_slim_v1.0.fits.gz' | md5sum --check
gzip -dc 4XMM_DR14cat_slim_v1.0.fits.gz > 4XMM_DR14cat_slim_v1.0.fits
python analysis/dr14_difference.py 5XMM_DR15cat_source_v1.0.fits 4XMM_DR14cat_slim_v1.0.fits results
python analysis/simbad_screen.py results/difference.csv results/difference_simbad.json
python analysis/erass1_screen.py results/difference_erass1.json results/difference.csv
```

The replay requires NumPy, SciPy and Astropy. See the main selection document for DR15 source provenance, units, survey footprint, and astrometric cautions. [4XMM-DR14 release](https://xmmssc.irap.omp.eu/Catalogue/4XMM-DR14/4XMM_DR14.html)

The J085618.5+380711 targeted eROSITA check uses a 30 arcsec DR2_Main cone and local aperture limits in 0.2–2.3 keV. The standard upper-limit conversion assumes an absorbed power law with Γ = 2 and N_H = 3×10²⁰ cm⁻². Run `python analysis/erass_followup.py 134.077254 38.119831 results/j085618-erass.json`. [eRODat API](https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/apis.html)

The HEASARC file's MD5 matched the release value. Of 1,176 objects passing the soft-flux and HR2 requirements, 1,041 had two-camera support, 824 passed the hard-band screen, and 765 passed the temporal criterion. Twenty-two had no DR14 slim neighbour within 15 arcsec. All 22 also lacked a DR15-reported probability ≥0.8 for an association with Gaia, WISE or OM. This differential sample overlaps the original four branches and must not be added to their row count without deduplication.

SIMBAD returned at least one indexed object inside 15 arcsec for 21 of the 22 sources. Of the 16 sources in the eRASS1-DE footprint, 13 have a catalogue neighbour inside 30 arcsec. The sole SIMBAD-empty position, J070328.5–371046, has an independent Gaia DR3 source (5566252434228157184) at about 3 arcsec with G = 15.16 and BP−RP = 3.14, so it does not have a clean optical blank field. The very brightest DR14-absent source, J081952.3–131932, is a previously studied eROSITA neutron-star candidate.
