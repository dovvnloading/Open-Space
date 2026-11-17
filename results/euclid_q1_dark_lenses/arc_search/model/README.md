# Pixel-level model comparison of the seven dark-centre triggers

**Five of seven triggers are predicted better by independently emitting knots than by a single lensed Gaussian source.** One blind trigger and the catalogue-adjacent trigger have lower pooled residuals for the lens model, with strong dependence on which region is held out. The adjacent trigger's images also fail the independent colour comparison in the [band-matched photometry](../validation/image_pairs.csv). No fit identifies a lightless mass concentration.

| Trigger | Proposal radius (arcsec) | Knots | Held-out mean squared residual, lens | Independent knots | Folds favoring knots |
|---|---:|---:|---:|---:|---:|
| [102020065_14336_54](102020065_14336_54.png) | 1.8 | 4 | 5.69 | 4.76 | 3/5 |
| [102020065_2048_5](102020065_2048_5.png) | 1.2 | 3 | 3.48 | 5.35 | 4/5 |
| [102020065_6144_6](102020065_6144_6.png) | 2.4 | 4 | 119.71 | 2.82 | 5/5 |
| [102022484_14336_8](102022484_14336_8.png) | 2.4 | 6 | 7.46 | 1.99 | 5/5 |
| [102022484_6144_2](102022484_6144_2.png) | 1.2 | 4 | 6.45 | 2.42 | 4/5 |
| [102022484_6144_44](102022484_6144_44.png) | 0.8 | 3 | 248.85 | 4.55 | 5/5 |
| [Near catalogue system 102157956_2704655017641535458](102157956_target_102157956_2704655017641535458_0.png) | 1.2 | 3 | 10.87 | 13.56 | 3/5 |

Each linked panel shows the observed VIS stamp followed by each fitted prediction and its standardized residual. Panels are representative fold-0 fits. The scores pool all five folds with held-out pixel counts as weights. Coordinates, fitted scales and boundary behavior are in [cross_validation.csv](cross_validation.csv).

**102020065_2048_5**, (63.90079259°, −49.70506435°), is the most consequential geometric ambiguity. Four folds favor knots. The fifth puts a bright area into the held-out region, with mean squared residuals 5.61 for the lens and 22.91 for knots. All five baseline lens fits reach the allowed ±0.35″ centre range. Enlarging this to ±0.75″ lowers the pooled lens error to 2.25, but the fitted isothermal scale then ranges from 0.75″ to 1.75″ across folds, and four individual folds still favor free knots. The scale and location are too unstable for a mass measurement. For **102022484_6144_2**, the wider centre changes 6.45 to 2.62, against 2.42 for knots. Only two of the five wider-centre folds favor knots. For the **catalogue-adjacent trigger**, every baseline lens fit reaches the centre bound. One held-out region gives very large error for both models. Simple model preference averaged over pixels cannot securely classify those objects.

## Fit and predictive check

* The lens maps each input angular position through a singular isothermal sphere plus constant external shear into one elliptical Gaussian emitter. We optimize SIS scale `b`, lens x/y offset, shear components, source x/y, intrinsic major-axis sigma, axis ratio and orientation. The fitted scale spans 0.55–1.55 times the proposed ring radius, with additional absolute limits 0.17–5″. Centre offsets range over ±0.35″ per coordinate for baseline fits, and each shear component over ±0.32. Four starting configurations reduce local-minimum sensitivity. The comparison uses one unconstrained-in-brightness elliptical Gaussian for each detected VIS knot, allowed to shift ±0.45″ per coordinate. Both descriptions have positive light amplitudes, two optional centrally fixed circular foreground Gaussians of sigma 0.14″ and 0.45″, and a free constant sky. The independent-knot model can approximate blends, compact groups and star-forming regions without identifying which interpretation applies.
* Both predictions integrate a 2×2 subpixel grid then convolve with the same circular Gaussian VIS PSF (assumed FWHM 0.17″). The fit uses the calibrated background-subtracted Q1 MER mosaic, flags and RMS. A circular aperture extends to the larger of 1.65 times the search radius or the search radius plus 0.55″, keeping the central pixels and room for predicted counterimages. Pixels with finite intensity and positive finite RMS below the sentinel threshold enter. Nonzero coadd flags can include usable exposures. Their fractions are retained in each JSON.
* Divide each stamp into 6×6 pixel blocks on its native 0.1″ grid. One fifth of the blocks are held out per fold. Exclude an additional one-pixel guard around test pixels from training. Every valid pixel is evaluated once across five independent parameter fits. Squared errors use diagonal native RMS. Such errors cannot be translated into likelihood odds using an independent-pixel interpretation. The sky and residuals are correlated, the PSF is approximate, and the candidate and initial knot positions were selected on the original image. This is conditional validation after candidate selection. A compact source hidden entirely in a test block can drive one fold's error.
* A fitted isothermal scale is a candidate model parameter. These comparison fits neither establish multiple imaging nor assign a scientific Einstein mass. They leave complex/multicomponent sources, additional deflectors, more general mass profiles and explicit spiral/ring morphologies untested. The known-lens [four-image control](../../catalog/refit_control.json) separately verifies a source-plane positional fit in a luminous system.

[Per-fold table](cross_validation_folds.csv) records held-out and training scores and all fitted parameters. Lens parameter order follows the description above, using arcseconds and radians. Knot parameter order is x, y, major-axis sigma, axis ratio and angle, repeated in brightness-ranked detection order. Values are in the image coordinate frame. Nonlinear shapes are fitted on training pixels. Their nonnegative amplitudes and the common background are solved by weighted linear least squares. In fold JSON, amplitudes are in native image counts, with the two shared central foreground components last. The fitted constants use the same units.

## Reproduce

From this directory, supply writable input and fold-cache directories. The program verifies the 21 raw VIS image, RMS and flag payloads against [the archived URLs and SHA-256 digests](vis_inputs.csv); a missing file is downloaded. The [initial knots](initial_knots.csv) are a reproducible subset of the band-matched photometry table. Source centres and proposed radii come from the parent search records. Dependencies are NumPy, SciPy and Matplotlib. The plotting and telescope cutout code operate in the frame of the native FITS WCS.

```bash
python fit_triggers.py /tmp/euclid-native --output /tmp/euclid-folds --jobs 2
python fit_triggers.py /tmp/euclid-native --output /tmp/euclid-folds --jobs 2 --center-bound .75 --ids 102022484_6144_2 102020065_2048_5
python summarize.py /tmp/euclid-folds
```
