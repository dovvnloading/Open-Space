# Gaia FPR search for previously unreported lensed quasars

**No new lensed-quasar candidate survived individual review.** We searched 174,737 Gaia quasar and quasar-candidate fields. Two overlapping apertures yielded 675 distinct image configurations whose colours matched and whose Gaia parallaxes and proper motions were consistent with zero. Seven passed both the two-filter lens-galaxy light test and the forward lens model. Five also had an independently catalogued quasar and no stellar or redshift conflict. Four of those five matched earlier lists. The fifth, PGC 2231560, remains an inconclusive quasar–galaxy system. Imaging was inconclusive for 242 configurations and unavailable for 80. The search recovered 2 of 235 eligible confirmed lenses. The null applies only to this selection. [Individual review decisions](individual_review.csv) · [Decisions for all 675 configurations](evidence/secondary/unique_science_configurations.csv) · [Combined results](evidence/secondary/combined_summary.json)

## Search and sensitivity

The October 2023 FPR tables contain 4,760,920 components in 3,760,032 environments. We selected fields with 2–8 components at absolute Galactic latitude ≥20°. Accepted components needed at least five Gaia detections, onboard G ≤21, and unique DR3 counterparts. Image separations ranged from 0.7″ to 5.2″. Both apertures used the same colour, astrometry, survey-image and mass-model tests. Legacy Surveys and Pan-STARRS fits tested for extended light between the point images in two filters. The mass model was a singular isothermal ellipsoid with external shear. [Selection thresholds, geometry and model](METHODS.md)

| Measurement | Primary aperture | Wider aperture |
| --- | ---: | ---: |
| Quasar fields examined | 174,737 | 174,737 |
| Fields passing Gaia cuts | 148 | 642 |
| Fields with a conclusive two-filter light test | 57 | 344 |
| Fields passing both the light and lens-model tests | 3 | 9 |
| Fields also having a quasar identification and no stellar or redshift conflict | 3 | 7 |
| Eligible benchmark lenses recovered through the quasar check | 1 / 235 | 2 / 235 |
| Random comparison apertures | 174,199 | 173,686 |
| Random apertures passing Gaia cuts | 1 | 5 |
| Random apertures passing both the light and lens-model tests | 0 | 0 |

The aperture columns count many of the same systems. Across either aperture, 31 of 235 eligible confirmed lenses passed the Gaia cuts, and 2 passed through imaging, modelling and quasar identification. The full benchmark contained 364 confirmed lenses. The recovered systems were SDSSJ134859.31+122653.2 and SDSSJ1537+3014. The requirement that multiple images have usable DR3 astrometry and colours already excludes many real lenses. The benchmark is an observed recovery fraction in this sample. [Known-lens ledger](evidence/secondary/known_recovery_either_aperture.csv) · [Primary audit and 95% Wilson intervals](evidence/audit.json) · [Wider audit and intervals](evidence/wide/audit.json)

The random apertures were offset at seeded random angles inside the observed FPR fields. We excluded apertures containing a catalogued FPR input. All five distinct configurations passing the random-aperture Gaia screen had inconclusive two-filter imaging. No control triggered a positive lens-and-galaxy classification. Their entry rates measure chance alignments under these local observing conditions. Uncatalogued quasars and shared Gaia measurements can occur in this control. [Matched masks and exclusions](METHODS.md) · [All primary fields](evidence/gaia_fields.csv.gz) · [All wider-aperture fields](evidence/wide/gaia_fields.csv.gz)

## Individual systems

Two image configurations passed the light and lens-model tests without a prior-list match. We retained both for individual review, including the field without a quasar identification. For each, the fitted colour difference and lens parameters are given below. The survey pixels favoured a freely positioned extended component over a single host fixed to either point image. Coordinates are ICRS degrees and the magnitudes are Gaia DR3 G. For these doubles the fitted Einstein radius, mass axis ratio and shear describe one regularized solution. A double supplies four position measurements for nine model and source parameters. Angles run east toward north. [Individual double measurements](evidence/secondary/reviewed_double_measurements.csv) · [Component astrometry, per-band fits and spectra](evidence/review_packets.jsonl)

| FPR field | Image coordinates | G | Separation | Example mass model |
| --- | --- | --- | ---: | --- |
| FW2917878464529177856 | (90.042375, −21.259147), (90.041951, −21.259799) | 18.473, 20.226 | 2.745″ | Centre (90.042155, −21.259709), θE 1.501″, q 0.618, PA 6.1°, shear 0.0546 at 102.1° |
| FW916413362313609472 | (125.796934, +43.888455), (125.798010, +43.887904) | 17.005, 18.453 | 3.425″ | Centre (125.797422, +43.888460), θE 1.930″, q 0.658, PA 17.0°, shear 0.0513 at 101.2° |

- **FW2917878464529177856.** The Legacy r and z fits give a relative image-colour difference of 0.034 mag. A bright point source sits beside a fuzzier companion, with residuals around both. Milliquas, DESI DR1 and SDSS DR19 supply no quasar identification. Gaia’s combined classifier assigns stellar probabilities of 0.869 and 0.99993. The field entered Gaia FPR through the Assef R90 infrared selection. Legacy DR10 and Pan-STARRS place their extended detection on the fainter Gaia source. A stellar blend or projection remains viable. [Data, fit and residuals](evidence/imaging/FW2917878464529177856-c8b7d9fb.png)
- **FW916413362313609472.** Pan-STARRS i and z give a relative image-colour difference of 0.215 mag. The two knots share a luminous envelope. Milliquas identifies PGC 2231560 as a radio-associated active galaxy at z ≈ 0.182. The lens fit predicts the fainter Gaia component would be more magnified. DESI target 39633148258026924 places a quasar at z = 0.182596 in the brighter component. Its measured fibre centre is 0.108″ from that image and 3.484″ from the fainter component. Gaia classifies the fainter source as a galaxy. Its flagged QSOC estimate is 0.178. Legacy’s extended galaxy is centred within 0.001″ of that fainter component, 2.515″ from the fitted lens centre. Pan-STARRS places its extended detection within 0.086″ of the bright component. The two Legacy catalogue releases use the same northern measurement here, a de Vaucouleurs galaxy with fitted half-light radius 1.05″. [Data, fit and residuals](evidence/imaging/FW916413362313609472-c1726fd2.png)

**Why PGC 2231560 remains open.** In Pan-STARRS i and z, a single extended profile with a freely fitted centre outscored a separate symmetric host on each nucleus. The paired-host minus central scores were +66.7 and +71.8 for Gaussian profiles, and +69.9 and +70.0 for exponentials. These are descriptive penalized-fit differences. An asymmetric or interacting host remains viable. The legacy SDSS spectrum 6380-56340-0196 also measured the brighter nucleus, at z = 0.18238. Both coadds show broad Hα blended with [N II], and DESI also resolves a broad Hγ feature. Neither samples the fainter object independently. The NED record for the fainter source identifies a galaxy in SDSS DR6. Smith et al. had already classified the bright object as a type 1.5 active nucleus at z = 0.182. The NED and SIMBAD references reviewed here identify the nucleus and its host. None identified either Gaia component as a repeated quasar image. Forty-three paired Gaia epochs show a zero-lag rank correlation of 0.398, with scan-dependent blending. MAST returned no covering Hubble image in the completed query. [Host-model residuals and scores](evidence/secondary/host-models/host_comparisons.json) · [Fibre positions, spectra and line fits](evidence/secondary/resolved_spectra/FW916413362313609472/spectral_evidence.json) · [Catalogue history and paper checks](evidence/secondary/object_history/object_history.json) · [Epoch measurements](evidence/secondary/focused/epochs/variability.json) · [Hubble archive check](evidence/secondary/hubble/coverage.json)

The low quasar redshift also constrains a lens interpretation. With the fitted 1.930″ Einstein radius, a spherical isothermal deflector would require a velocity scale of at least 259 km/s. At hypothetical deflector redshifts of 0.05 and 0.10 the scales rise to 305 and 389 km/s, using flat ΛCDM with H₀ = 70 km/s/Mpc and Ωₘ = 0.3. The deflector redshift and velocity dispersion are unmeasured. A resolved spectrum of the fainter object would distinguish a second image of the broad-line source from light dominated by a companion galaxy. [Conditional velocity scales](evidence/secondary/isothermal_velocity_scales.csv)

A separate test followed 4 distinct configurations with a positive single-filter residual, compatible two-filter colours and a quasar identification. None supplied the required combination of independently deblended central light and separated foreground and background spectra. A single quasar spectrum can occur in a star–quasar projection. Two consistent quasar redshifts are also possible in a physical pair. The field files retain isolated stellar spectra, discrepant redshifts, prior pair identifications, galaxy deblends, fibre positions and the uncertainties that prevent a classification. [One-band results](evidence/secondary/secondary_summary.json) · [Every follow-up decision](evidence/secondary/field_adjudication.csv) · [Focused field checks](evidence/secondary/focused/focused_summary.json)

## Published systems and data provenance

The 21 primary lists cover GLQ doubles, quads, binaries and nearly identical quasars, CASTLES, GraL III/VI/VIII/IX/X, Gaia’s 381 published FPR candidates, Lemon et al., DESI Dawes, He, Hsu and Jing, the DESI strong-lens value-added catalogue, DULAG, Lenscat, Euclid Q1 and ECLIPSE. A prior candidate or published quasar pair vetoed a new claim. We matched within 10″ of each image, with 6″ for the supplemental close-pair and SuGOHI lists. The Chen close-pair catalogue contributed 4,112 entries, the fixed SuGOHI VI release 2,000, and the authors’ full SuGOHI compilation 3,961. For the two unlisted fields we also checked all 1,867 entries in the Wu et al. golden sample, all 5,286 in their superset and all 3,140 Makarov–Secrest candidate dual active nuclei. None fell within 8″ of either image. [Archived tables and independent matches](evidence/secondary/object_history/object_history.json)

The published binary and nearly identical quasar comparison samples produced no two-band lens-model trigger in either aperture. [Known-pair comparison](evidence/secondary/wider_published_comparisons.csv)

The field FW3793908336710193920, for example, coincides with the published DESI quasar pair J1142−0110. Object-specific reference searches cover the unlisted fields that warranted review. Catalogues and object-linked bibliographies have finite sky and literature coverage. [Catalogue status](evidence/catalog_status.csv) · [Publication matches](evidence/literature_matches.csv) · [Supplemental receipts](evidence/supplemental_receipts.json) · [Individual literature queries](evidence/secondary/literature/individual_literature.json)

Milliquas v8 contributed 355 component–catalogue matches, DESI DR1 contributed 291 and SDSS DR19 contributed 168. The DESI comparison read all 28,425,963 rows of the native zall-pix-iron FITS file and verified the publisher’s SHA-256 digest, `2d95ad99361039b556c402b49e0e7c84df5f00106dc5731d44476a58b128b49b`. These counts include repeated fields and spectra. The complete publication snapshot and input catalogue receipts are pinned to a commit. SHA-256 values for Gaia queries, catalogue downloads, survey FITS stamps and targeted follow-up requests are recorded alongside the corresponding measurements. [Frozen publication snapshot](evidence/published_snapshot.json) · [Catalogue bytes and checksums](evidence/catalog_input_checksums.csv) · [DR3 query receipts](evidence/gaia-dr3-queries.json) · [Wider DR3 queries](evidence/wide/gaia-dr3-queries.json) · [Image inputs and fits](evidence/imaging/imaging.jsonl)

| Gaia FPR file | Downloaded SHA-256 |
| --- | --- |
| lenscand.dat.gz, detected components | `bc1b277a3e099ab4b5c8b477bb09bb479a539adddfad55dcb2e761a8d96895eb` |
| lenscatn.dat.gz, input catalogue associations | `8b92e5b3fb41393d406bc7cb2e99ccf46e1a25103d2145b6685e39a1081e05c3` |
| ReadMe, release dictionary | `3479cc127c08d34cc7b3f5a7c44a97066f39b565947b96f5201abc9073204a40` |

## Reproduce the search

The Gaia extraction, corrected survey-pixel fits, full DESI crossmatch, object-specific follow-up and final catalogue checks were deposited by successive runs of the same branch workflow. [Gaia and corrected image analysis](https://github.com/dovvnloading/Open-Space/actions/runs/45007382313) · [DESI release scan](https://github.com/dovvnloading/Open-Space/actions/runs/45012253872) · [Host models and spectra](https://github.com/dovvnloading/Open-Space/actions/runs/45014076367) · [Final verification](https://github.com/dovvnloading/Open-Space/actions/runs/45015023850)

The branch workflow records the acquisition and analysis commands. A workflow-only push on the research branch can run the entire calculation, including new image fits. The existing complete literature snapshot is replayed for the prior-publication test. [Workflow](../../.github/workflows/gaia-fpr-quasar-lenses.yml) · [Gaia release receipts](evidence/fpr-downloads.json) · [Output checksums](OUTPUTS.sha256)

```bash
git switch research/gaia-fpr-quasar-lenses
workflow=.github/workflows/gaia-fpr-quasar-lenses.yml
printf '\n# Replay the complete Gaia FPR study\n' >> "$workflow"
git add "$workflow"
git commit -m '[fresh-fpr] Replay the Gaia FPR study'
git push origin HEAD
```

To check the committed evidence and repository tests from the repository root, use

```bash
study=results/gaia-fpr-quasar-lenses
(cd "$study" && sha256sum -c OUTPUTS.sha256)
python "$study/audit.py" "$study/evidence"
python "$study/audit.py" "$study/evidence/wide"
python -m pip install '.[archive]' numpy pandas scipy matplotlib
python -m open_space.verify
python -m unittest discover -s tests
```
