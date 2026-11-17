# HST comparisons at the JWST source positions

The HST pixel checks did not establish an additional co-moving counterpart. The 2021 HD 53143 STIS images have substantial occulting-wedge and PSF residuals. The independently published Gemini/NICI measurement is the useful earlier astrometric test for that field's bright source. Several HST rolls have excursions at other predicted coordinates, but a single high pixel or local maximum carries little weight in these correlated residuals.

## Coverage

The [ESA Hubble Science Archive queries](query_receipts.json) intersect a 20.16″ circle around each host with archived image footprints. The [responses](receipts/) retain all returned rows and the checksums refer to their uncompressed response bytes. They include different product levels and exposures of overlapping fields. The [footprint list](hst_footprints.csv) yielded zero imaging rows for LAWD 66, 897 for HD 53143, 61 for HD 38206, 36 for HD 21997 and 122 for HD 202628. The zero applies to this archive and query.

| Host | HST pixels used for the positional tests |
|---|---|
| HD 53143 | Twelve STIS rolls from October 2021, GO 16202, plus the ALICE NICMOS F160W and F110W products from 1998 and 2006. The two predicted positions of the bright source lack usable support in the cropped NICMOS maps. |
| HD 38206 | Two ACS/HRC F606W coronagraphic products in GO 10539, November 2005. |
| HD 21997 | Two ACS/HRC F606W coronagraphic products in GO 10539, August 2006. The archive calls the target HR 1082. |
| HD 202628 | Six STIS WEDGEA1.0 rolls from May and September 2015, GO 13786. |

## Measurement and limits

The [input manifest](hst_inputs.json) records the byte length, URL and SHA-256 of each of the 24 FITS products. For every NIRCam F444W threshold peak in the four fields, and for the auxiliary F200W/F1140C triggers in HD 21997 and HD 202628, the merged Gaia propagation code predicts the position of a source fixed relative to the host and of a source at fixed sky coordinates. The propagation includes stellar proper motion and annual parallax using the geocentre. The HST comparison uses fixed projected separation for the companion hypothesis. An orbit can change projected separation over a decade.

The [716 track rows](hst_tracks.csv) give the predicted host-relative position, response at the nearest native pixel, maximum in a 0.15″ search radius, offset to that local maximum, and whether the aperture is supported. ACS positions use the product WCS and propagated host position. STIS positions use the target acquisition reference pixel and the native WCS rotation. For STIS, a detector-space median of other rolls more than 8° away suppresses the stellar halo. We exclude flagged STIS pixels and low-weight ACS pixels and apply a compact difference of Gaussians. Local scores use a 0.22–0.48″ annulus and its median absolute deviation. NICMOS uses the published ALICE SNR map. A supported search radius can contain an invalid central pixel, which leaves the forced score blank.

These values are positional diagnostics. Roll-dependent residuals, optical/infrared colour, WCS and centering systematics, and different PSFs affect them. The HST stacks have no independent completeness calibration here. The native JWST injection trials determine the infrared sensitivity. No optical non-detection supplies a planet exclusion. For example, the HD 21997 F200W trigger's fixed positions have forced ACS scores of −0.4 and +1.4. The 2015 fixed positions of the two HD 202628 F444W triggers have finite forced STIS scores of −0.33 to +0.15 and −0.98 to +0.77, respectively (see the ledger for masked centers).

## Replay

From this directory, run `python fetch_hst_images.py CACHE` and `python test_hst_tracks.py CACHE`. Run `python hst_inventory.py` to requery the archive. The measurement uses the frozen inputs and the neighboring field folders. The STIS visit selections follow [Faramaz et al. 2019](https://doi.org/10.3847/1538-3881/ab3ec1) and [Stark et al. 2023](https://arxiv.org/abs/2304.07370). The earlier processed NICMOS images are from [the ALICE archive](https://archive.stsci.edu/prepds/alice/).
