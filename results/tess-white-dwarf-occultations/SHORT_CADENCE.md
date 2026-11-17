# No periodic deep dip in 1,106 short-cadence light curves

The frozen inventory contains 1,693 SPOC target-pixel products on 146 catalogue counterparts. The supplementary search used each product’s associated 20- or 120-second simple-aperture light curve. Exact TIC identity and a SPOC crowding fraction of at least 0.80 left 1,106 products on 86 white dwarfs (52,368,692 quality-zero finite flux samples). The identity cut excluded 13 entries, and the crowding cut excluded 574. The search found zero periodic triggers and 142 isolated template hits in 25 products on 16 stars. The hits occupy 22 star-sector visits. Multiple widths and simultaneous cadences can select the same excursion. [Frozen product inventory](target_pixel_products.csv) · [Photometry receipts, rejection reasons and every hit](short_cadence/roster_screen.json.gz)

| Observation period | Searched 20-s products | Searched 120-s products |
|---|---:|---:|
| Before Cycle 5 | 165 | 333 |
| Cycle 5 onward | 304 | 304 |

The event-seeded screen fits 2–32 minute boxes with a minimum width of two exposures, a 6-hour median baseline and local shoulders. Periodic triggers require three seed events of at least 2% depth and score 4, plus the stacked depth, score and recurrence criteria. Isolated triggers require a 5% shoulder-to-centre contrast and score 10. The larger selection uses pipeline apertures for screening. The primary Sector 109 sample also received a blind 1.5-pixel aperture search on all 16 of its target-pixel products. [Detector](fast_roster.py) · [Sector 109 target-pixel search](short_cadence/pixel_search.json)

## Pixel follow-up

For all 142 hits, the 120-second target-pixel check uses a separate 1.5-pixel aperture centred on the proper-motion-corrected Gaia position. Twelve of the 22 visits pass the full FFI-style Gaia field check. The other ten have a lower worst-case target fraction or at least one Gaia neighbour without a G magnitude. They contribute diagnostic images only. For every supported window the record gives a signed difference image, the source and photocentre locations, the aperture depth, the predicted dilution, and the strongest modelled neighbours. The same spacecraft observed the original and follow-up pixels. [Pixel measurements](short_cadence/roster_pixels.json.gz) · [Extractor](roster_pixels.py)

The 120-second pixels supplied the centre and both shoulders in 131 windows. In every one the centre was within 2.54% of or above the median flux in a wider local window. The large contrasts came from high shoulders or an off-target source. Only GD 1212 and TIC 14802607 had shoulder-relative contrasts above 5% in the 120-second small aperture. GD 1212 is a known pulsating white dwarf. Its observed oscillations persist outside the trigger windows. [Stellar identifications and pulsation measurements](SHORT_LITERATURE.md)

The native 20-second target pixels supply all three windows for all 124 hits originally seen at 20 seconds, including the 11 without complete 120-second coverage. The largest apparent template deficit in TIC 14802607 coincided with an object moving at 39.43 pixels per day. Relative to quiet comparison frames the aperture brightened by 0.54%. A simultaneous fit for the moving source and the white dwarf gave a 0.86% WD loss, sensitive to the fitted profile width. [Fast-pixel results and moving-source fit](FAST_PIXEL_VETTING.md)

## Replay

Run from the survey directory with NumPy, SciPy, Astropy and tess-point installed. The screen fetches the hashed SPOC products directly. Pixel remeasurements prefer the 120-second product when both cadences exist in a visit. They use the archived 3-arcminute Gaia field when available and query VizieR for the remaining fields. The rates in the main README come from separate injections into independently selected FFI light curves. The larger short-cadence screen has no calibrated occurrence denominator.

```bash
python fast_roster.py target_pixel_products.csv tic_match.csv work/fast
python roster_pixels.py target_pixel_products.csv work/fast work/pixelcache work/rosterpixels
python archive_followup.py --root .
```
