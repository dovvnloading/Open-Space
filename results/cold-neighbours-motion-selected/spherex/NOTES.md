# SPHEREx DR1 spectral measurements

WISE 0855−0714 is a positive control at a supplied trajectory. The in-field T9 WISE J041358.14−475039.3 is a separate control selected by the blind CatWISE screen. Seven additional CatWISE tracks received exploratory follow-up. One, J181125, has published astrometry and an early-Y estimate in the Kirkpatrick et al. 20 pc census.

## Inputs and scope

The spectral images are SPHEREx DR1, DOI [10.26131/IRSA678](https://www.ipac.caltech.edu/doi/10.26131/IRSA678). Each per-frame row records the IRSA cutout URL, the fitted position and the SHA-256 of the downloaded multi-extension FITS. The manifests freeze the frame selection.

For WISE 0855 the exposure and pointing list comes from the [Cold Movers snapshot](https://github.com/WilmerGaspar/cold-movers/blob/58ad127da016c2e6d12704b368036f3b25356418/data/precomputed/wise0855/photometry.csv), SHA-256 8bd7a250c152163fcabfcbab46c013e8663ba6d0f1e645ecf702467821081b68. Its 254 rows identify 248 distinct DR1 image products. The five specified wavelength windows selected 62 input rows independently of upstream flux, flags or S/N. Two repeated the same exposure and byte-identical measurement for this target. Counting each target–image pair once gives 60 DR1 measurements, of which 52 passed quality checks. Extraction, spectral summaries, model comparisons and plotting reject repeated target–image pairs. Separate targets can use the same native exposure. We measured the flux anew from the IRSA pixels at the snapshot's epoch-specific positions.

For the seven exploratory tracks, the archive search covered the DR1 wide and deep collections at a position advanced 10.5 yr from epoch 2015.405. The SIA constraint was 3.10–5.0 µm with a 0.001 degree cone. The T9 lookup used a 0.0005 degree cone centered at (63.4930, −47.8430). Detector coverage bounding boxes give a superset of the wavelengths at the source. After deduplicating parent URLs, each detector 4–6 was sampled using up to 48 evenly spaced indices from its time-sorted frame list. The actual wavelengths and widths come from the FITS WCS-WAVE at the target position. This sampling uses part of the available deep-field DR1 data. The CatWISE RA motion is multiplied by sec(dec) when converting to RA in degrees. The candidate trial tracks assume zero parallax.

An IRSA TAP cone around each CatWISE 2015.405 position (0.035 degrees) supplied AllWISE neighbors from allwise_p3as_psd using designation, ra, dec, w1mpro and w2mpro. WISE 0855 used (133.76, −7.2423), radius 0.02 degrees and ordering by W1. Positions were held fixed. AllWISE J041358.14−475039.4 is the known T9 itself and was removed from that control's stationary background list. The seven exploratory trials retained every returned neighbor as a competing stationary source. Closely spaced objects can be unresolved, and an uncataloged background source can bias the fitted amplitudes.

## Reduction

The fit estimates one signed target flux, a local planar background and an independent amplitude for each AllWISE source within five detector pixels of the target. It evaluates the nearest local, oversampled ePSF at native pixel phases with cubic interpolation and unit integrated response. The WCS Jacobian converts MJy/sr to µJy. The fit uses the VARIANCE plane, keeps SOURCE and FULLSAMPLE flags, removes every other flagged pixel and rejects frames with central quality flags or insufficient usable pixels. Frame errors include inflation by the square root of reduced χ² when it exceeds unity. The local plane absorbs the background. Absolute fluxes remain subject to point-source flux calibration, PSF systematics and unidentified neighboring sources. The [IRSA ePSF tutorial](https://caltech-ipac.github.io/irsa-tutorials/spherex-psf/) describes the geometry and sampling.

The five bands retain signed measurements and inverse-variance weighting. Visit jackknife errors group frames by gaps greater than 20 days, then omit each group in turn. A missing visit jackknife means there was no gap that long. For continuously sampled polar fields the additional jackknife uses fixed 60-day blocks anchored to MJD 60000. The number of occupied blocks is reported for each window. Bands average different spectral channels. A rise from 3.3 to 4.5 µm alone cannot identify a specific molecule, classify a new dwarf or determine its effective temperature.

The shape comparison fits a single nonnegative normalization per template over the same qualified measurements from 3.1 to 4.97 µm. The T7/T8/T9 field templates and individual Y0/Y1 spectra are from the [Gagné et al. SPIFF QR2 archive](https://zenodo.org/records/19051216). The two Y spectra are WISE J1738+2732 and WISE J0350−5658. The comparison includes all 25 [Polletta SWIRE galaxy and AGN templates](https://www.iasf-milano.inaf.it/~polletta/templates/swire_templates.html) over redshifts 0 through 4 in 0.1 increments. Templates are converted from F_lambda to F_nu shape and integrated over a rectangular approximation of each frame bandwidth. Template uncertainties, frame covariance, extinction and emission-line diversity beyond the grid are unmodeled. Reported minimum χ² values compare the specified grids descriptively.

## Reproduction

From this directory, with numpy, scipy and astropy installed, run the following. The first pass retrieves the cutouts. Local copies can then be reused.

~~~bash
python fit_dr1.py wise0855_frames.csv wise0855_fits.csv --cache fits_cache
python fit_dr1.py shortlist_frames.csv shortlist_fits.csv --cache fits_cache
python fit_dr1.py wide_frames.csv wide_fits.csv --cache fits_cache
python fit_dr1.py additional_frames.csv additional_fits.csv --cache fits_cache
python summarize.py additional_fits.csv additional_bands.csv
python template_check.py additional_fits.csv additional_templates.csv --templates templates
python summarize.py wide_fits.csv wide_bands.csv
python template_check.py wide_fits.csv wide_templates.csv --templates templates
python summarize.py wise0855_fits.csv wise0855_bands.csv
python summarize.py shortlist_fits.csv shortlist_bands.csv
python template_check.py wise0855_fits.csv wise0855_templates.csv --templates templates
python template_check.py shortlist_fits.csv shortlist_templates.csv --templates templates
~~~

The six-panel [DR1 spectral figure](dr1_spectra.png) uses 0.12 µm display bins. Each marker is at the inverse-variance weighted wavelength, with the formal error of the weighted mean. The vertical scale is independent in each panel. Reproduce it with the following command.

~~~bash
python plot_dr1.py wise0855_fits.csv shortlist_fits.csv wide_fits.csv dr1_spectra.png
python plot_dr1.py additional_fits.csv additional_spectra.png
~~~

J181125 denotes J181125.34+665806.4, J182532 denotes J182532.09+671906.2, and J041358 denotes J041358.20−475038.1. J085014 denotes J085014.24+541123.2 and J092244 denotes J092244.78+640537.8. The three additional follow-ups use full CatWISE designations to distinguish nearby objects. Their [spectral figure](additional_spectra.png) uses the same binning. The field search and WISE image checks set the evidentiary status of these specific tracks.

## Background sensitivity in the three additional fields

The additional targets were remeasured with quadratic spatial terms (x², xy, y²) added to the sky plane. Source, neighbors, pixels and weighting remain the same. At 4.40–4.62 µm, the planar and curved fits give +73±36 and +34±39 µJy for J010638.19+122908.3, −46±43 and +55±45 µJy for J213348.06+391040.9, and −300±89 and +83±95 µJy for J192049.57+111222.6. For J192049 even the curved fit has median image-level reduced χ² of 16.3, and a galaxy template fits at least as well as a known Y dwarf. The negative planar amplitude is evidence of background-model misspecification. J010638 separates a cold template from the best galaxy by less than 1.1 χ² with either background. J213348 favors a galaxy template with the curved background. These results give no robust spectral classification for the three tracks. The planar 3.1–4.97 µm figure can be compared to the corresponding curved-background band and model tables.

~~~bash
python fit_dr1.py additional_frames.csv additional_quadratic_fits.csv --cache fits_cache --quadratic-sky
python summarize.py additional_quadratic_fits.csv additional_quadratic_bands.csv
python template_check.py additional_quadratic_fits.csv additional_quadratic_templates.csv --templates templates
~~~
