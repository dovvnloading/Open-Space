# Individual-exposure tests for CEERS-002 and COSMOS-009

Both catalog coincidences fail independent exposure checks for an extragalactic transient. These tests address two individual positions.

## CEERS, two-band brightening at the edge of detector NRCA3

The June 2026 F140M and F210M position (214.813705°, +52.822765°) is 0.180″ from a pre-existing HST galaxy, 3D-HST AEGIS #25530. Its EAZY photo-z is 0.8131 (68% interval 0.790–0.839). The earlier observation is JWST program 7814 observation 029, and the later one is 8559 observation 006. All contributing exposures in the two filters use NRCA3.

| Filter | Early dither 2 | Early dither 3 | Early dither 4 | Late dither 3 | Late dither 4 |
|---|---:|---:|---:|---:|---:|
| F140M, nJy | 1,174 | 796 | 1,084 | 1,169 | 1,174 |
| F210M, nJy | 1,835 | 1,135 | 1,755 | 1,835 | 1,810 |

Each value is sky-subtracted flux measured in a 0.16″ radius aperture after local translation to the corresponding epoch's science mosaic. Non-finite pixels were assigned the stamp median for registration and diagnostic photometry. In a clipped exposure this procedure gives a flux deficit, so the dither-3 values indicate lost signal and cannot be used as total source fluxes.

The earlier third dither cuts across the central galaxy in both filters. Its native detector coordinates are x = 5.95 or 6.12 pixels, putting the trigger 0.187″ or 0.192″ from the physical edge. Native CRF science/DQ planes show 19.5% of aperture pixels flagged DO_NOT_USE and non-finite in F140M and 17.3% in F210M. Both other early dithers have 100% finite aperture pixels. The untruncated early dither 2 reproduces the later fluxes. The apparent two-filter rise depends on the coverage of the earlier data. [10 exposure fluxes](ceers_june-002/dithers.csv) · [Native detector positions and DQ fractions](ceers_june-002/native.csv) · [FITS context and contributing exposures](ceers_june-002/coverage.json) · [Byte receipts](ceers_june-002/receipts.json) · [Host and redshift source](ceers_june-002/host.json)

## COSMOS, a point source with a long-baseline positional shift

COSMOS-009 at (150.1458367°, +2.2836319°) has a pre-existing point-like counterpart. 3D-HST assigns COSMOS #10232 a stellar flag of 1. Registering the historical HST catalog to the 2026 F182M data on 40 nearby bright nonstellar matches gives a displacement of −0.299″ east and −0.149″ north. None of the 40 controls has as large a displacement. The robust per-axis scatters were 0.045″ and 0.022″. Eight 2026 filters give similar shifts. The two JWST epochs alone give shifts comparable to local control scatter. The long-baseline result and source morphology favor a foreground star, potentially blended. The redshift of the 0.30″ COSMOS2020 counterpart cannot be assigned to this source. [HST-to-JWST astrometry](cosmos_medium-009/hst_motion.json) · [JWST-to-JWST local controls](cosmos_medium-009/local_motion.csv) · [3D-HST catalog](https://doi.org/10.1088/0067-0049/214/2/24)

The two filters yield 22 single-exposure measurements in total. Three early dithers per band are on NRCB2. The later eight per band are on NRCB2, NRCB3 and NRCA2 and span seven days. All 22 apertures have 100% finite registered-input coverage within 0.16″. Retaining NRCB2 in both epochs removes most cross-detector differences. Medians from its three dithers per filter show no consistent decline.

| Filter | Radius | Earlier NRCB2, nJy | Later NRCB2, nJy |
|---|---:|---:|---:|
| F162M | 0.16″ | 1,334.6 | 1,352.5 |
| F162M | 0.28″ | 1,420.4 | 1,472.7 |
| F182M | 0.16″ | 864.4 | 865.8 |
| F182M | 0.28″ | 952.9 | 952.5 |

[Individual dither photometry](cosmos_medium-009/dithers.csv) · [FITS context and dates](cosmos_medium-009/coverage.json) · [Byte receipts](cosmos_medium-009/receipts.json)

## Reproduction

Run these commands from this directory after generating the parent search's catalog and mosaic stamp caches. They use the shared archive reader and aperture procedure from the parent search. The context bit map and HDRTAB identify the contributing filenames. The detector-level check uses the original CRF science and DQ extensions.

```sh
python exposure_followup.py ceers_june-002 /path/to/stamps /path/to/ceers-output --dq
python exposure_followup.py cosmos_medium-009 /path/to/stamps /path/to/cosmos-output
python cosmos_astrometry.py /path/to/catalog-cache /path/to/cosmos-output
```

The included [local HST comparison table](cosmos_medium-009/nearby_hst.tsv) preserves the VizieR response and query URL. Its catalog memberships and the neighboring-source astrometric test can be reproduced without another large image download.
