# Extended-source response and non-exclusive spectral alternatives

The 100 K NASA PSG CO2 shape, selected from the earlier four-temperature December comparison, remains consistent with the small-aperture January 2, 2026 measurements under several spatial weighting approximations. This is evidence of compatibility, not a unique chemical identification or measured temperature.

## Spatially integrated response sensitivity
We evaluate true detector-coordinate band-center maps, channel-specific fiducial curves and photon-weighted model responses at each image pixel out to 12 native pixels. Profiles include the embedded ePSF, and PSF-convolved softened rho^-1 and rho^-1.5 coma distributions. The response is averaged with the light weights inside radii 2.5 and 3.5. A second variant includes signed annular weights corresponding to mean sky subtraction at 8–12 pixels.

This assumes a spatially invariant spectrum and a separable observed light profile. It is not a complete response-before-optical-convolution model. The actual extraction uses nonlinear median sky, so the mean-sky variant is only a sensitivity check. Truncated PSF wings, spatially varying gas/dust ratios, neighbors and full radiometric uncertainty remain unresolved.

January heldout chi-square (two residual degrees of freedom, amplitude and constant refitted):
- Radius 2.5: **4.86–4.90**, versus center-only 4.879.
- Radius 3.5: **8.34–8.42**, versus center-only 8.354.

Thus this spatial response approximation neither removes larger-aperture tension nor materially changes the small-aperture compatibility. Full output: `spatial-response-test.json`; code: `examples/swan_spatial_response.py`.

## Alternative broad empirical shapes
We also fit a symmetric Lorentzian and triangular emission profile through the same spectral responses. Each has a free center and width selected with December 22/25 measurements; separate amplitude and constant per epoch make six free parameters for six training points. The fits therefore interpolate the training data and do not constitute evidence. Centers were bounded to 4.23–4.29 microns; widths to 0.0005–0.08 microns. January amplitude and constant were refitted with the shape frozen.

| Shape | Radius 2.5 January chi-square / 2 DOF | Radius 3.5 January chi-square / 2 DOF |
|---|---:|---:|
| Lorentzian |338.28|396.32|
| Triangle |19.13|48.37|
| Earlier Gaussian |44.80|91.61|
| Physical CO2 100 K |4.879|8.354|

The molecular template extrapolates better than these particular symmetric shapes. This does **not** exclude flexible asymmetric features, mixtures, contaminants or systematics. Models have different training flexibility and physical priors; this table is not a Bayesian model probability or exhaustive chemical-identification test. The triangle fit requires width about 0.0561 microns, within the final search interval. An earlier narrow width bound was enlarged before reporting results.

See `ALTERNATIVES_AND_INDEPENDENT_DATA.md` for documented fixed CO and H2O alternatives and the bounded external-data search. Neither failed alternatives nor absence of indexed papers proves novelty. No independent contemporaneous infrared confirmation has been retrieved.
