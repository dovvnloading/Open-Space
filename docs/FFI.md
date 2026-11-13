# Coordinate-only TESS full-frame image search

The two pilot fields were chosen before measuring their images. One field is centered at ICRS (17.69°, −60.34°), near the known eclipsing binary BB Tuc. The second is at (17.97°, −60.28°), offset from the first field without reference to a catalog or a light curve. Both use 21 × 21 TESS pixels in Sectors 1 and 2. The known object supplies a position and period check *after* blind extraction; the measurement program accepts neither a target ID nor a period. The adjacent field samples only a small region of sky.

## Sky coverage and detector

The input is a position on the sky, a square width and two increasing sectors. The code downloads the MAST TESScut SPOC product for each sector and checks the detector and observation metadata. TESScut generates FITS checksums with request timestamps, so separately hash the science arrays, WCS, observation metadata and ordered parent-frame identifiers. The output preserves both this stable scientific digest and the exact byte checksum of the retrieved FITS file. Repeated cutouts from one detector can share parent frames. [MAST TESScut](https://mast.stsci.edu/tesscut/) · [MAST cutout example](https://spacetelescope.github.io/mast_notebooks/notebooks/TESS/interm_tesscut_requests/interm_tesscut_requests.html) · [Archive reader](../open_space/ffi_archive.py)

1. Retain finite time samples with quality = 0 and at least 98% valid image data. Fit diffuse background gradients, common illumination and first-order translations of the median sky image using the faintest 80% of available pixels. Reject large residuals when refitting the nuisance parameters. Remove residual trends using a running baseline within contiguous time segments.
2. Evaluate 0.18–3 day sinusoidal models at every eligible interior pixel. The two-pixel border provides local comparison pixels. At least 99% of measurements for an eligible pixel must be usable. Keep a peak with ≥8% explained variance and an approximate amplitude/noise ratio ≥8. Group adjacent detections under the strongest local amplitude-weighted peak. A two-minute target list, source catalog or preselected ephemeris is unnecessary.
3. At each peak's preferred period, fit a **second instrument model** using background pixels outside its 7 × 7 footprint. Project each pixel's signal onto the peak's phase and use the positive local flux to find the variable's sky centroid. Flag a strong opposite-sign lobe, diffuse signal, high correlation with detector regressors, or disappearance after withholding the source. This step was developed on the BB Tuc control after the first fit leaked light into the pointing correction. The adjacent field had already been inspected at that time.
4. Search the other sector from scratch and match by sky position. A recurrent detection requires a unique spatial pair within 31.5 arcsec, periods within 2.5%, and acceptable spatial diagnostics in both images. A single-sector detection, frequency mismatch and veto remain in the archive evidence. These are screening thresholds. The search does not provide a survey-wide false-positive probability or resolve components within a TESS pixel. [Detector code](../open_space/ffi.py)

## Blind recovery and controls

Thirty-minute Sectors 1 and 2 supplied 2,527 distinct parent frames for the two pilot locations on the same detector, producing 5,054 cutout time rows. The 17 × 17 eligible grid per cutout yielded 1,156 location/sector period searches. In the development field the dominant peak landed 5.80 and 4.30 arcsec from BB Tuc; the two centroids agreed to 2.38 arcsec. The independent TESS image periods were 0.371645 and 0.371627 days. Gaia's published short-timescale period is 0.371490 days. This supports recovery of its established variability. The neighboring field had zero first-sector peaks and two second-sector peaks. [Pixel search](../results/ffi-fields.json) · [BB Tuc independent Gaia measurements](PIXELS.md)

Synthetic Gaussian PSFs, each carrying one imposed sinusoid, were added to the calibrated images in the neighboring field before processing. Signals of 1.0, 2.5 and 5.0 electrons/s peak-pixel semi-amplitude at 0.411, 0.832 and 1.69 days were recovered at the expected period in both sectors, with offsets of 0.65, 0.24 and 0.09 arcsec. Independently permuting the two sectors' whole-frame time order produced no above-threshold peaks. These controls use a single sky location, three joint amplitude/period combinations and one pair of shuffles. They sample a narrow part of sensitivity and false-alarm behavior. [Control measurements](../results/ffi-controls.json)

## Six untouched fields in the Large Magellanic Cloud

After finalizing the detector on the two pilot fields, the next six centers were set on a regular sky-coordinate grid, with no reference to a star catalog. The two rows have declinations −69.4° and −69.2° at right ascensions 77.0°, 77.4° and 77.8°. The program examined 3,468 location/sector pairs, using 4,270,842 eligible pixel/cadence measurements. These six fields share the same camera/CCD in each sector, supplying 2,527 distinct detector-frame files for their 15,162 cutout time rows. Taken together with the two original fields on another camera, the original two-sector search uses **5,054 distinct full-frame detector files**. [Selection](../data/ffi-lmc-grid.json) · [Every image peak and veto](../results/ffi-lmc-grid.json)

| Outcome on first-sector peaks | Count |
| --- | ---: |
| Recurrent by image position, period and shape | 13 |
| Image artifact veto | 4 |
| Single-sector peak | 2 |

OGLE published period and position checks were applied to these 13 locations after the image search. The crossmatch includes 125,131 Large Magellanic Cloud variables in six OGLE-IV catalog families. A match requires ≤30 arcsec plus ≤2% agreement of the measured period with half, once or twice a catalog period. **Twelve locations have a consistent counterpart within 5.6 arcsec.** The closest consistent counterparts include six eclipsing binaries, an ellipsoidal variable, four classical Cepheids and an RR Lyrae star. Two locations also admit another consistent object within 30 arcsec. All neighbors and period discrepancies are preserved. At some locations the strongest sinusoid runs at twice the catalog orbital or pulsation frequency. [OGLE catalog](https://ogle.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/) · [Individual matches and digests](../results/ffi-ogle.json) · [Reproducible catalog check](../examples/check_ffi_ogle.py)

The remaining location at ICRS (77.05502°, −69.23710°) showed 1.11221 and 1.11694 day modulation, with independent centroids 6.30 arcsec apart. The OGLE-IV sample contains RRc OGLE-LMC-RRLYR-05589 at 13.41 arcsec (0.2476053 days) and ellipsoidal OGLE-LMC-ELL-01519 at 22.38 arcsec (220.436685 days). Neither of these catalog periods satisfies the check. This signal calls for independent imaging and resolved ground-based photometry. Other catalog classes, aliases, amplitudes varying across observations and multiple stars within a TESS resolution element remain relevant.

## Reproduce and search new fields

```sh
python -m pip install -e '.[archive]'
export OPENBLAS_NUM_THREADS=1
python -m open_space.ffi --plan data/ffi-fields.json --cache runs/ffi --output runs/ffi-pilot.json
python -m open_space.ffi_controls --cache runs/ffi --output runs/ffi-controls.json
python -m open_space.ffi --plan data/ffi-lmc-grid.json --cache runs/ffi --output runs/ffi-lmc.json
python examples/check_ffi_ogle.py --cache runs/ogle --evidence runs/ffi-lmc.json --output runs/ffi-ogle.json
```

New field plans use the same JSON structure. Set a unique field ID, ICRS sky center, a field role and width, and sectors in increasing order. A first run can use product records with only a sector. Copy the reported science digests to the plan to pin a replay. Use moderate tile sizes and limit simultaneous cutout requests. For automatic dispatch import `open_space.ffi.run_field(field, cache_dir)`; its `sources` list contains sky coordinates, per-sector frequencies, instrument vetoes and the recurrence decision. The second sector's unmatched peaks are also recorded.

## Independent images of the 1.11-day location

We directed the same field search at Sectors 3, 27 and 28, then evaluated the original 1.112209-day frequency at the original sky position. Sector 3 contributed 1,077 quality-accepted thirty-minute images. Sectors 27 and 28 contributed 3,351 and 3,449 accepted ten-minute images about two years later. Within each sector all detector pixels received the original blind period search. A separate targeted fit used the frequency measured in Sector 1; this second fit allows the amplitude and phase to vary by sector. The independent peak and nearest coordinate pixel can sample different parts of a blended PSF. [Frozen coordinates, frequencies, witness and products](../data/ffi-source-revisit.json) · [All five sectors and phase bins](../results/ffi-source-revisit.json) · [Later blind field search](../results/ffi-lmc-revisit.json)

| Sector | Independently strongest local period (days) | Centroid distance from S1 (arcsec) | Fixed frequency variance fraction at nearest coordinate pixel |
| ---: | ---: | ---: | ---: |
| 1 | 1.11221 | 0.0 | 0.014 |
| 2 | 1.11694 | 6.3 | 0.469 |
| 3 | 1.11065 | 5.0 | 0.376 |
| 27 | 1.09528 | 12.4 | 0.325 |
| 28 | 1.12022 | 12.8 | 0.434 |

The signal recurs in the independent sectors without a preset period, on two CCDs and at two different cadences. The sampling window power at the original signal frequency is <0.001 in every sector. The fixed coordinate's nearest Sector 1 pixel has substantial correlation with detector regressors (0.79), while the neighboring discovery peak retains 0.35 explained variance after estimating the instrument model independently of the source footprint. In the four subsequent sectors the nearest coordinate pixel's correlations with those regressors are 0.12–0.15. This checks several specific instrumental explanations. The detector frequency resolution, phase jitter and crowding still limit a physical interpretation.

The known same-field Cepheid OGLE-LMC-CEP-1027 offers a control on the later images. It reappears at 1.2185 and 1.2145 days in Sectors 27 and 28 and its centroid shifts by only 2.4 and 2.6 arcsec. It fails the original blind spatial criterion in Sector 3; the fixed-frequency measurement there remains strong. The much larger late shift of the 1.11-day location is under one 21-arcsec TESS pixel, and requires consideration of a blend or observing-geometry changes. The two cataloged neighbors are roughly 11–22 arcsec from the individual sector centroids. Current TESS data identify a recurring location and timescale. Sub-pixel source identity and physical classification require independent evidence.

```sh
python examples/revisit_ffi_source.py --plan data/ffi-source-revisit.json \
  --cache runs/ffi --output runs/ffi-source-revisit.json
```
