# Published identifications and FGK comparison sample

Three of the original seven Hephaistos WISE excesses have strong source identifications. JWST resolved the galaxies responsible for D and E. Radio interferometry resolved a background AGN beside G, and the longer-wavelength WISE centroids coincide with it. The remaining sightlines have differing evidence for contamination. [Zackrisson et al. 2026, Hephaistos IV](https://arxiv.org/abs/2607.09460) · [Ren et al. 2025](https://doi.org/10.1093/mnrasl/slaf006)

## Hephaistos sightlines

Suazo et al. selected seven M dwarfs from about five million Gaia stars within 300 pc. The reported excesses were strongest in W3/W4 at 12/22 µm. The 100–275 K temperatures in Table 5 describe fits to unresolved photometry. Later optical spectroscopy finds ordinary main-sequence properties and little evidence of youth in the spectroscopic subset. A rotation-period search supplies no sensitivity limit for isolated optical dips. [Suazo et al. 2024, Tables 4–5](https://doi.org/10.1093/mnras/stae1186) · [Korn et al. 2026, Hephaistos III](https://arxiv.org/abs/2607.25701)

| Target | Published evidence and interpretation |
|---|---|
| A | Nearby radio source J1245−2652 has a poorly constrained spectral slope. Our calculation below yields subarcsecond W1/W2 separations. Contamination is plausible. The dominant MIR emitter remains unidentified. |
| B | Long-wavelength WISE centroids approach a radio source with spectral slope 0.63 ± 0.11 under the convention Sν ∝ ν^(−α). Strong evidence for contamination. The contaminant lacks a confirmed extragalactic redshift. |
| C | VISTA resolves a second source 3.75 arcsec away along the displaced MIR emission. Strong evidence for contamination. The nature of the nearby object and any physical association remain open. |
| D | MIRI resolves a point-dominated galaxy at z ≈ 0.9 within about 1 arcsec. Its spectrum resembles a hot dust-obscured galaxy. The galaxy explains the WISE excess. |
| E | MIRI resolves an extended, clumpy galaxy at z ≈ 0.4 about 1 arcsec away. Its spectrum resembles a dusty starburst. The galaxy explains the WISE excess. |
| F | Weak WISE astrometric evidence and an elongated W3 image. Deep optical/near-IR imaging has yet to identify the object producing the MIR light. |
| G | VLASS J233532.86−000424.9 lies >5 arcsec from the epoch-corrected star. Compact emission with a brightness temperature >10⁸ K identifies an AGN. W3/W4 centroids favor the AGN. |
| H, I, J | Hephaistos III added these lower-significance targets. All have W3/W4 catalog S/N below 3.5. H and J show possible interlopers or displacements in weak imaging. I lacks an obvious archival contaminant. |

The original A/B/G radio associations, the resolved G radio observations, the subsequent epoch-corrected optical/WISE investigation and the D/E JWST observations are distinct evidence. The archival null for D illustrates that centroids alone can miss a background galaxy. [Ren et al. 2024](https://doi.org/10.3847/2515-5172/ad5017) · [Ren et al. 2025](https://doi.org/10.1093/mnrasl/slaf006) · [Ren et al. 2026](https://arxiv.org/abs/2607.03619v6) · [Zackrisson et al. 2026](https://arxiv.org/abs/2607.09460)

Zackrisson et al. previously extracted **QR2** SPHEREx spectra of D and E. Stellar light dominates to approximately 5 µm and accounts for roughly 80% in both F560W images. The stellar shares at 10 µm fall to about 10% for D and 20% for E. A rest-frame 3.3 µm band would be observed near 6.3 and 4.6 µm at these galaxy redshifts. [Zackrisson et al. 2026, §§2–4](https://arxiv.org/abs/2607.09460)

## Reproduced astrometric check

An independent coordinate check changes the evidence for A. Using the printed WISE centroids and mean epochs (Tables 4–6) and direct Gaia DR3 positions, we get 0.181 and 0.111 arcsec in W1/W2 for A, compared with 1.60 and 1.48 in Table 4. For A–F, 3/24 separations are within 0.05 arcsec of the printed values at the stated epoch. Evaluating the star 16 years earlier matches 24/24 within 0.0052 arcsec. For H/I/J, 12/12 agree within 0.012 arcsec at the stated epoch. At the corrected epoch, B and C retain W3 offsets of 2.68 and 3.77 arcsec. The radio and MIRI imaging results do not depend on this coordinate check. No remeasurement of image centroids or statistical re-fit of Ren's controls was attempted. The earlier public audit identified the same pattern. Run `python centroid_check.py` on `ren_centroids.csv` to reproduce `ren_offsets_recomputed.csv`. [Ren et al., v6](https://arxiv.org/abs/2607.03619v6) · [Pemberton & Lilitari, reproducibility audit](https://gitlab.com/codalilitari/seti-reproducibility-papers/-/releases/v1.0)

In a separate November 6 erratum to the **2024** paper, Ren et al. corrected the hot dust-obscured galaxy density to about 2.5 × 10⁻⁹ per square arcsecond. In their simple alignment calculation across five million stars, radii of 6 and 12 arcsec give expected counts of about 1.4 and 5.6. The TGSS separation for G was corrected to 5.773 arcsec. Neither the expected counts nor this low-resolution TGSS association can identify each WISE source. [Ren et al. 2026, erratum](https://doi.org/10.3847/2515-5172/aebe03)

In `candidates.csv`, `ra`, `dec`, `pmra`, `pmdec` are Gaia DR3 at J2016.0 in ICRS degrees, with μα* and μδ in mas/yr. Gaia identifiers, literature W3/W4 S/N and labels come from Hephaistos II Table 5 and III Table 3. G magnitudes, parallaxes and `ren_table1` coordinates retain the precision of Ren et al. Table 1. Direct Gaia astrometry underlies the epoch propagation. Read source IDs as strings or exact integers. [Gaia DR3](https://doi.org/10.1051/0004-6361/202243940) · [Ren et al., Table 1](https://arxiv.org/abs/2607.03619v6)

## Other published searches

The targeted FGK sample contains all 53 optical/near-IR selected main-sequence candidates flagged by Contardo & Hogg for excess unWISE W1/W2 flux among 4,898,812 stars. Sixteen published fitted temperatures fall within 100–600 K. These temperatures and luminosity fractions in `contardo_hogg.csv` are hypotheses about unresolved emission. The first three entries in Appendix A Table 2 were specified as the SPHEREx subset before spectrum inspection. All three lie beyond 300 pc. The two stars within 300 pc were added as a distance-complete supplement. A separate extension subsequently measured all other 48 targets and their matched comparison stars. [Complete 53-star spectral analysis](../fgk_spectra/README.md) [Contardo & Hogg 2024](https://doi.org/10.3847/1538-3881/ad6b90) · [author catalog and code](https://doi.org/10.5281/zenodo.13171690)

Both FGK sources with Gaia parallax > 3.333… mas have published natural explanations or disqualifying evidence.

| Gaia DR3 ID | Distance by reciprocal parallax | Follow-up |
|---|---:|---|
| 3183786280737192704 | 255.5 pc | TIC 43488669, a runaway with optical dips and a debris-disk interpretation. Later TESS/ATLAS work recovered a 4.14-day periodic signal alongside the irregular dips. [Tajiri et al. 2020](https://doi.org/10.3847/1538-4365/abbc17) · [Greenberg 2026](https://doi.org/10.3847/2515-5172/ae7322) |
| 2942418533272205312 | 264.3 pc | TYC 5940-1510-1, an extreme debris-disk host. Narrow-band mid-infrared photometry across the 10 µm feature is consistent with silicate-dominated emission. [Moór et al. 2021](https://doi.org/10.3847/1538-4357/abdc26) · [Moór et al. 2024](https://doi.org/10.1093/mnras/stae155) |

Their temperatures in Contardo & Hogg Table 2 are 815 and 835 K. Earlier IRAS/LRS surveys selected bright 100–600 K emitters, including evolved stars and planetary-nebula candidates. A newer machine-learning study ranked photometric outliers using earlier selections in training. This analysis supplies no new completeness estimate for either survey. [Carrigan 2009](https://arxiv.org/abs/0811.2376) · [Mignone & De Carolis 2026](https://doi.org/10.1016/j.eswa.2026.131232)

## FGK measurements

**The W1/W2 photometric excess repeats for all 53 FGK stars in CatWISE.** All 53 have one qualifying CatWISE counterpart within 2 arcsec. Each is >0.10 mag brighter in both bands than the published stellar prediction. Median predicted minus measured magnitudes are 0.208 in W1 and 0.346 in W2. Restriction to the 51 stars with one CatWISE detection within 6 arcsec gives 0.213 and 0.356 mag. Across all 53, median corrected unWISE minus CatWISE magnitudes are +0.002 and −0.016. The multiple reductions use overlapping WISE exposures. A repeated brightness excess alone cannot identify its origin.

The authors' AllWISE matches give upper limits for 15 stars in W3 and 32 in W4. `fgk_baseline.csv` preserves quality flags, photometric uncertainties, source IDs and the published stellar predictions. The data receipt pins the author's full machine-readable catalog to its Git commit, blob and SHA-256. `fgk_astrometry.csv` stores the Gaia coordinates, motions and parallaxes. [Author catalog](https://github.com/contardog/NotATechnosignatureSearch/tree/fc257c9c3fd54f4edb7e94d61cc5b9a20b27bc2b) · [AllWISE column definitions](https://irsa.ipac.caltech.edu/data/WISE/docs/release/AllWISE/expsup/sec2_1a.html)

Year-11 NEOWISE returned 1,121 neighborhood extractions. At least one qualifying visit was available for 52 stars in W1 and 46 in W2. Only four W1 and three W2 stars had a second qualifying visit. None of those first-to-last differences exceeds 5σ under the uncertainty model below. This short, incomplete baseline cannot support a general claim of constancy or a variability-rate upper bound. [NEOWISE final release](https://irsa.ipac.caltech.edu/data/WISE/docs/release/NEOWISE/)

### Association and uncertainty

CatWISE2020 and NEOWISE were queried separately in 12-arcsec cones centered on epoch-propagated Gaia positions. For CatWISE we require one motion-fit position inside 2 arcsec of the Gaia position propagated to MJD 57170, clear `ab_flags` and `cc_flags` in the band, and S/N ≥ 5. The comparison separately counts additional catalog sources within 6 arcsec. Stellar predictions and native unWISE photometry adopt −0.004/−0.032 mag W1/W2 corrections for the AllWISE scale. The 0.10 mag offset is a descriptive threshold. [CatWISE definitions](https://irsa.ipac.caltech.edu/data/WISE/CatWISE/gator_docs/catwise_colDescriptions.html) · [Contardo & Hogg, §III.1](https://arxiv.org/abs/2403.18941)

For NEOWISE the query spans December 13, 2023 through August 1, 2024 (60291 ≤ MJD < 60524). Source rows must lie within 2 arcsec of the per-exposure propagated position. We reject duplicate associations within an exposure, require frame quality >5, image-quality factor >0.5, SAA separation >5 degrees and a single fit component. Per band we require clear moon/artifact flags, zero saturation, quality A/B, S/N ≥5 and profile χ² ≤3. At least four observations with all consecutive gaps ≤30 days make a visit. The reported magnitude is the visit median. Its error combines a 0.03 mag systematic in quadrature with the larger of the formal median error and robust scatter-based median error. `fgk_neo_coverage.csv` reports first-to-last changes only when two visits qualify. Selection on individual detections can bias faint visit medians bright.

### Reproduction and limits

```sh
python fgk_baseline.py
python fgk_fetch.py --kind catwise --cache /tmp/fgk-cache
python fgk_fetch.py --kind neo --cache /tmp/fgk-cache
python fgk_measure.py
python centroid_check.py
```

Use Python 3.10+ in this directory. `fgk_*_receipts.json` records request URLs, epochs of retrieval, row counts, byte counts and digests. Each compressed response bundle contains the exact per-target csv or VOTable payload. `fgk_measure.py` verifies per-response hashes before reduction. Live services can change response metadata or row order.

Image-level deblending, companion exclusion and optical monitoring for every distant FGK member remain incomplete. The exposure-detection approach also leaves low-S/N forced NEOWISE photometry open. These 53 targeted stars were already selected on WISE and add no blind-survey exposure.
