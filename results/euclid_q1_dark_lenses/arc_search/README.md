# Independent VIS ring search

This pilot selects sky positions without a catalogue target or central-galaxy requirement. The frozen fields use zero-based entries 45, 120, 190 and 300 from the lexically sorted public Q1 MER tile directory, spanning EDF-S, EDF-F and EDF-N. Each tile contributes eight 384-row strips at the positions in fields.json. Partial and uncovered regions count only when the detector can use them. No tile was chosen using a candidate location. The scope is bright, nearly closed rings. Open arcs, doubles and sub-resolution rings have little or no sensitivity here.

Pre-run detection design: 2×2 binning of the 0.1 arcsec VIS mosaic; radii 0.8, 1.2, 1.8 and 2.4 arcsec; coherent azimuthal emission in at least six of eight sectors; a faint central aperture; local background and noise calibration. Visually inspect any passing configurations and use independent NISP photometry before assigning an evidence level. Artificial rings on real backgrounds test the response of this geometry only.

## Observed sky and seven configurations tested

**Six independent sky configurations and one catalogue-cutout configuration passed both VIS filters. None establishes a shared background source behind an undetected deflector.** The filter tested 27,013,228 eligible 0.2″-square binned centres, or 0.083374 deg² after excluding margins and unusable neighbourhoods, across 32 non-overlapping strips in four widely separated Q1 MER tiles. Of the 1,258 initial maxima, six have an azimuthally extended ridge that is brighter than its inside and outside. This second test was added after inspecting the first-stage trigger distribution. The result is exploratory.

| ICRS RA, Dec (deg) | Trial radius | Native H core S/N | Evidence at the trial position |
|---|---:|---:|---|
| 63.834995, −49.362202 | 1.8″ | −0.43 | The brightest diametrically opposed knots have incompatible four-band colours. A free-knot model gives lower pooled held-out VIS error. [Four bands](validation/102020065_14336_54.png) |
| 63.900793, −49.705064 | 1.2″ | 3.75 | Ambiguous in VIS. A common-source fit wins on pooled error but loses four of five individual folds. Moving the allowed centre changes the recovered radius from 0.75″ to 1.75″ across the wider-bound fits. Two main knots have mutually disjoint MER photo-z 95% intervals, 0.90–1.06 and 0.09–0.64. A MER object lies 0.37″ from the trial centre. No secure shared source or mass centre. [Four bands](validation/102020065_2048_5.png) |
| 63.948638, −49.587271 | 2.4″ | −0.30 | Independent knots predict held-out pixels better in every fold. Multiple resolved galaxies have different colours and photo-z solutions. [Four bands](validation/102020065_6144_6.png) |
| 66.370888, −46.862138 | 2.4″ | 3.45 | Six VIS peaks surround an empty-looking space between bright galaxies. Independent knots win all five folds. Several pairings have compatible colours, but the model fails to reproduce the configuration. [Four bands](validation/102022484_14336_8.png) |
| 66.483179, −47.089540 | 1.2″ | 0.25 | Three of four peaks have MER counterparts with broad, overlapping photo-z intervals. A common colour fits several pairs. Independent knots have lower pooled prediction error, 2.42 versus 6.45; widening the lens-centre prior lowers the lens error to 2.62. Lensing remains unverified. [Four bands](validation/102022484_6144_2.png) |
| 65.920077, −47.090762 | 0.8″ | 5.27 | J and H both detect central flux. A MER entry lies 0.18″ from the trial centre. Independent knots win all five folds. [Four bands](validation/102022484_6144_44.png) |
| 270.465429, +64.153227 | 1.2″ | 6.41 | This is the separate catalogue-cutout trigger. The trial centre is about 1.15″ from the parent MER object and contains J and H light. Several neighbouring MER detections and chromatic structure complicate the putative images. The pooled VIS fit alone is inconclusive. [Four bands](validation/102157956_target_102157956_2704655017641535458_0.png) |

The first and second-stage filters use the native VIS mosaic at 0.1″ per pixel, binned 2×2. Each first-stage test compares three samples per angular sector at the trial radius with samples at 1.65 times that radius. At least six of eight sectors must exceed their empirical 2.2σ contrast and their average must exceed 3.4σ. The core must lie below 2.5σ. The second stage samples 36 angles, requiring positive radial contrast in at least 19 and angular coverage in at least three quadrants. The code records shape concentration and the local core constraint. Near-zero RMS, exact-zero fill, chip edges and inadequate sky estimates limit individual steps. The wide-area eligibility definition uses finite nonzero VIS science pixels, not individual-exposure masks. Native RMS and flags are retrieved for each shortlist location.

The native VIS/Y/J/H cutouts use their own zero points and WCS. Central apertures have radius 0.3″; the second aperture uses the trial Einstein radius. Photometric errors take the larger of the local blank-aperture MAD and the propagated diagonal RMS. Mosaic flags are ORed across contributing exposures; we record the flag union and nonzero fraction in each aperture. Some flagged measurements can retain valid contributing pixels. Colour comparisons PSF-smooth to an approximate 0.5″ FWHM, use two apertures, and fit a common four-band spectrum with one relative amplitude and a 5% per-component error term. Colour agreement alone occurs among unrelated galaxies. The likelihoods omit deblending and PSF-model uncertainty. Neither broad photo-z overlap nor absence in DESI establishes a shared source. [Photometry](validation/photometry.csv) · [Pairwise colour fits](validation/image_pairs.csv) · [MER associations](validation/component_associations.csv) · [Image-model validation](model/cross_validation.csv)

A common-source SIS with external shear and an elliptical Gaussian source was compared on the same VIS pixels with separate Gaussian emitting knots. Fits use five spatial folds with a guard band. Other source morphologies, intrinsic spiral arms, tidal structures, collisional rings, asymmetric source colours and more complex lens potentials remain possible. We have no spectroscopy classifying those individual morphologies. The failed shared-source, centre or colour tests are sufficient to withhold a discovery label.

## Sensitivity and limits

The 1,491 ray-traced injections have no central galaxy. The injected source has elliptical Gaussian surface brightness, a 0.06 external shear and a 0.17″ VIS PSF. Each lands on actual Q1 pixels. Source peak brightness is specified before PSF convolution in units of local binned-pixel noise. Recovery requires the correct discrete radius and a centre within 0.4″. These fixed positions and radii are repeated across scenarios, so the counts are sensitivity measurements, not independent sky trials.

| Source offset / Einstein radius | Recovery at peak 3σ | Recovery at peak 8σ |
|---:|---:|---:|
| 0.02 | 127 / 248 | 213 / 249 |
| 0.18 | 28 / 248 | 150 / 249 |
| 0.45 | 0 / 248 | 0 / 249 |

A separate geometry-only test recovered 24/82 and 72/85 thin circular rings with peaks of 2σ and 5σ. The drop to zero for 497 substantially offset sources precludes an abundance limit for lenses producing open arcs or doubles. Results refer to nominal plateau coverage and the reported noise, masks, scales and models. [Ray-traced injection measurements](lens_injections.csv) · [Circular controls](ring_injections.csv)

For all seven positions the H luminosity bound and hypothetical projected M/L floor change with lens redshift and light profile. Each entry in the scenario grid treats the trial radius as an Einstein radius and assumes a source at z=2. The direct bound assigns all light in that radius to a possible foreground galaxy and uses max(0, flux) + 5σ. Profile-conditioned limits use the central aperture with circular Sérsic n=1 or 4, half-light radius 0.2″ or 0.8″, and approximate Gaussian PSF. Those limits apply to attenuated observed light in the redshift-shifted band; a more extended, dusty or higher-redshift deflector can evade the compact-source test. For (66.483179°, −47.089540°), H central flux is 0.25σ. Its all-light H apparent bound is 22.62 AB; at assumed z_l=0.6, the corresponding M/L floor is 93.9. The common-source lens interpretation failed validation, so this is a sensitivity illustration. [All four-band, redshift-dependent limits](validation/conditional_limits.csv)

## Reproduction

From this directory, Python 3.10+, NumPy, SciPy and Pillow suffice for the search and photometry. For the model comparison see its script and the modelling notes. The cuts and source keys are frozen in fields.json. The ZIP contains the per-strip measurements and control results as small JSON files. Each strip receipt verifies the exact raw FITS header and science byte range. No survey-sized download is required for offline inspection. [Strip receipts](strips.csv) · [All VIS triggers](vis_triggers.csv) · [Native FITS inputs](validation/native_inputs.csv) · [Frozen stage records](frozen_stages.zip)

```sh
# Replay all native strips and both filters. Alternatively, unzip frozen_stages.zip to inspect or reuse the recorded stages.
python search.py
python followup.py
for tile in 102020065 102022484 102044820 102159484; do
    python lens_injections.py --tile "$tile" --out raw_ray_trials
done
python summarize.py --ray-tests raw_ray_trials

# Recheck the 188 catalogue parents without a usable positive model radius.
python targeted.py --cache raw_zip_indexes
python fallback.py /path/to/q1_discovery_engine_lens_catalog.csv
python photometry.py --cache raw_native
python associations.py
python limits.py
```

Archive data and instrument details are documented by [ESA](https://euclid.esac.esa.int/dr/q1/dpdd/merdpd/dpcards/mer_bksmosaic.html) and [IRSA](https://doi.org/10.26131/IRSA601).
