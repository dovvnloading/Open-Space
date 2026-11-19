# Frozen source register and scientific caveats

Access date: 2026-11-18 America/New_York (2026-11-19 UTC). This file records *public source metadata*, not evidence that all links were successfully downloaded.

| Product | Authoritative URL | Frozen reference / caveat |
| --- | --- | --- |
| COSMOS-Web DR1 public-release description | https://cosmos.astro.caltech.edu/page/cosmosweb-dr | 0.54 deg² 4-filter NIRCam; MIRI F770W only ~0.2 deg² |
| COSMOS2025 paper | https://doi.org/10.1051/0004-6361/202555799 | Shuntov et al., A&A 704 A339, published online 2026-01-06; 784,016 matched sources |
| Original catalogue homepage | https://cosmos2025.iap.fr/ | Official **v1.1 dated 2026-05-04**, revised photometry and photo-z for MIRI tiles B5, B9, B10 |
| Original catalogue description | https://cosmos2025.iap.fr/catalog.html | FITS `COSMOSWeb_mastercatalog_v1.1.fits`: PHOTOMETRY, LEPHARE, SE++APER, CIGALE, ML-MORPHO, B+D, GALIGHT-MORPHO; download access requires website form |
| Original usage tutorial | https://cosmos2025.iap.fr/tutorials/catalog_usage_tutorial.html | Separate large LePHARE PDFz pickle with redshift grid, Bayesian PDF and likelihood per source; separate LePHARE/CIGALE SED files. Those are **not joint p(M*,z)** samples |
| CDS/VizieR archive | https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A%2BA/704/A339?format=html&tex=true | Files `lephare.dat` = 388-character fixed width × **784,016** rows; `phot.dat` = 5,563-character fixed width × 784,016. VizieR inserted 2026-01-06, last revised 2026-03-02. **Not established identical to May 2026 official v1.1**. `screen.py` hashes the actual archive bytes, if transfer succeeds |
| Published COSMOS2025 quality / high-z SED | https://doi.org/10.1051/0004-6361/202555799 | LePHARE medians/68% intervals are **marginal posteriors**; star/AGN χ², Chandra flag, blending/artefact flags and multi-band fluxes help triage but do not validate source masses |
| Boylan-Kolchin 2023 | https://doi.org/10.1038/s41550-023-01937-7 | Use rare-halo cumulative mass function, max baryon fraction `f_b`, survey comoving volume, and stellar mass *density* comparison. The existence of one object above an *expected* max is **not a deterministic inconsistency** |
| ST99 / Colossus | https://bdiemer.bitbucket.io/colossus/lss_mass_function.html | Sheth99 only defined for FOF in Colossus; Tinker08 200m. Mass definitions and high-z extrapolation are substantive systematics |
| STScI Dawn JWST Archive (DJA) | https://dawn-cph.github.io/dja/ | Search actual released NIRSpec redshifts matched to coordinates; a catalogue photo-z does not constitute spectroscopy |
| MAST JWST holdings | https://mast.stsci.edu/ | Independent NIRSpec observation metadata, reductions and line identifications must be checked per candidate |
| COSMOS2025 image/SED browser | https://cosmos2025.iap.fr/fitsmap/ | Per-object inspection PDFs and tiles for blend, artifact, star, AGN, and optical low-z diagnostics |

## Frozen analysis choices
- Survey geometry is **nominal** 0.54 deg² NIRCam solid angle for each shell, not a validated depth/selection/mask. MIRI 7.7µm coverage is only ~0.2 deg² and must **not** be treated as an F770W nondetection elsewhere.
- Redshift shells (lower exclusive, upper inclusive): (6,7], (7,8], (8,9], (9,10], (10,12], (12,15]; do not double-count boundary objects. Any z>15 is to be reported separately rather than silently discarded.
- Flat ΛCDM: H0=67.32 km/s/Mpc, Ωm=0.3158, Ωb=0.0493, ns=0.96605, σ8=0.8120; imposed maximum baryon allotment **f_b=0.16** (requested convention). Stellar IMF, mass return and stellar-population model must be harmonized before interpreting masses.
- HMF: Sheth–Tormen 1999 (FOF halo mass) with Tinker 2008 (M200m) as sensitivity; Colossus 1.3.9. These *different* halo definitions must not be equated. Compute `lambda(>M)` by integrating each model over the full shell and define `P(Mmax ≤ M)=exp(-lambda)` (Poisson baseline). Halo cosmic variance and calibration uncertainty are **not** included.
- Compare both `lambda=1` (expected-one threshold) and upper CDF quantiles (95%, 99%, 99.9%) of the distribution. For a galaxy with posterior `p(M*,z)`, an **extreme-halo check** is only a screening statistic. Conclude a discovery only from forward-modeled population under the identical selection and a viable complete posterior, not from 68% intervals or a single threshold crossing.
- Catalog screening in `screen.py` uses Type=0 and `(zfinal>6 OR zPDF-u68>6)` as a **broad triage**, not a full posterior selection. A low-z-dominant source with a secondary high-z peak outside its central 68% range can be missed. Any result that hinges on this screen is explicitly **not a complete high-z null**.

## Pre-specified source-vetting stop conditions
For **every candidate whose independently inferred joint redshift–stellar-mass posterior could put it above even the model-sensitive 95%-maximum-halo baryon bound**, require:
1. Record exact catalogue release, FITS row ID, J2000 position, matched image tiles, foreground extinction, aperture/SE++ photometry and calibrated filter curves, original source-data SHA256, depth/upper limits and masks.
2. Independently run EAZY and a physical-posterior sampler (Prospector or Bagpipes) with documented sampler convergence, stellar+nebular line and AGN prescriptions, evolving dust and SFHs, uncertainties/covariances; retain per-candidate joint `p(M*,z)`, compare prior and IMF sensitivity, test catastrophic photo-z alternatives.
3. Fit separate dusty z≈2–4 galaxy and brown-dwarf stellar libraries, flexible [O III]/Hβ/Hα line strengths, reddened AGN continua/torus, and blended-source forward models. Inspect NIRCam detections and nondetections and image-cutout segmentation.
4. Map whether MIRI F770W coverage exists and evaluate actual measured flux/SNR or upper limit; inspect Chandra and wider X-ray/radio association and flags (not just LePHARE `Type`).
5. Query DJA/MAST and primary spectroscopy literature by source coordinates, exposure/quality and actual line ID (NIRSpec or independent instrument). A failed search is **no spectroscopy located**, not spectroscopic exclusion.
6. Propagate **joint** mass-redshift posteriors into a pre-registered forward halo+selection model; account for Eddington bias, look-elsewhere effects (number of tested galaxies / bins), survey masks, cosmic variance, mass systematics, lensing. Competing halo-mass definition and HMF prescriptions must not disagree qualitatively.
7. Independent referee reruns the arithmetic, checks actual raw-data hashes and source classifications, and explicitly signs off on the exact claimed level. **No [REMARKABLE]** absent independent discovery-grade validation; any unmet item remains in the pending ledger.

Status at document creation: only primary-release documentation and frozen scientific criteria are evidenced. No observational refit or posterior-probability measurement has been performed.
