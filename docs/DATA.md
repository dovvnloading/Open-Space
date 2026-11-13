# Observational data

## First target

The bundled observation is an SDSS BOSS spectrum of a starburst galaxy at RA 214.316940°, Dec +29.275123°. Plate 6497, MJD 56329 and fiber 139 identify the observation. SDSS DR17 pipeline `v5_13_2` reports a good plate and a galaxy redshift of 0.08566851 with pipeline uncertainty 0.00000748 and zero redshift-warning bits. This emission-line galaxy provides positive controls for wavelength handling and line recovery. A population search requires a separate, prespecified target selection. [Original FITS](https://data.sdss.org/sas/dr17/eboss/spectro/redux/v5_13_2/spectra/lite/6497/spec-6497-56329-0139.fits) · [CSV](../data/sdss-6497-56329-0139.csv) · [Provenance and hashes](../data/sdss-6497-56329-0139.json)

The coadd contains 4,624 observed-frame vacuum-wavelength pixels spanning 3565.33–10337.13 Å. Of those, 4,277 pass the project's pixel mask and inverse-variance checks. Four familiar features provide a direct cross-check. The observed column below gives the largest eligible flux pixel within ±4 Å of the prediction. Its precision is limited by pixel sampling and the instrumental profile. [SDSS wavelength reference](https://www.sdss4.org/dr17/spectro/spectro_basics/)

| Transition | Rest vacuum (Å) | Predicted from survey redshift (Å) | Observed peak pixel (Å) |
| --- | ---: | ---: | ---: |
| Hβ | 4862.683 | 5279.26 | 5279.59 |
| [O III] | 4960.295 | 5385.24 | 5385.18 |
| [O III] | 5008.239 | 5437.29 | 5437.51 |
| Hα | 6564.614 | 7126.99 | 7126.89 |

## Recreate the input

Run these commands from the checkout in a Python environment with Astropy. The bundled CSV itself can be read offline with the standard library. The conversion was verified with Astropy 8.0.2.

```sh
python -m pip install astropy
python -m open_space.archive --output work/replay
```

To replay from previously downloaded bytes, add `--fits path/to/spec-6497-56329-0139.fits`. The converter checks the pilot source SHA-256 before reading FITS data. It verifies the observation identity in two HDUs, the reduction version, and the wavelength medium. The JSON records the checksums of both the original FITS and the exact CSV. Replaying identical inputs in the same directory is supported. Conflicting outputs require a fresh destination.

For another DR17 BOSS/eBOSS observation using this reduction, supply its plate, MJD and fiber together. Each command retrieves one lite coadd. The first retrieval establishes a digest. For future executions use the recorded digest, on both downloaded and local copies. For local files the URL is constructed from the supplied identifiers, so compare against a digest from an independently retrieved archive product to establish byte-level equivalence.

```sh
python -m open_space.archive --output work/observations \
  --plate "$PLATE" --mjd "$MJD" --fiber "$FIBER"
python -m open_space.archive --output work/replay \
  --plate "$PLATE" --mjd "$MJD" --fiber "$FIBER" \
  --fits "$FITS_PATH" --sha256 "$SOURCE_SHA256"
```

## Data contract

| CSV column | Source and interpretation |
| --- | --- |
| `wavelength` | `10**COADD.loglam`, in observed-frame vacuum Å on the SDSS heliocentric wavelength scale. |
| `flux` | Calibrated flux density in 10⁻¹⁷ erg cm⁻² s⁻¹ Å⁻¹. No continuum or extinction correction is applied at ingestion. |
| `ivar` | Inverse flux variance in reciprocal squared flux units. Zero weight excludes the pixel. |
| `mask` | The analysis exclusion flag, copied from `and_mask`. Any nonzero value excludes a pixel. |
| `and_mask` | SDSS SPPIXMASK bits flagged in every exposure contributing to that pixel. |
| `or_mask` | SPPIXMASK bits flagged in at least one contributing exposure. Retained for inspection and sensitivity checks. |

The converter writes all wavelengths in their original order, without interpolation. Nonfinite flux is replaced by a zero-flux, zero-weight placeholder. Nonfinite or nonpositive inverse variance is set to zero. SDSS documents that inverse variance already suppresses its worst pixels and recommends the AND mask for coadds. Requiring an entirely clean OR mask would remove the peaks of all four positive controls in this spectrum. A flagged contributing exposure is still relevant evidence. For any candidate, rerun with OR flags excluded, inspect the combined sky and per-exposure residuals, and account for the exposure count and instrumental resolution. Individual exposures, dispersion, sky levels and pipeline fit tables are available in the original full products. [Pixel mask definitions](https://www.sdss4.org/dr17/algorithms/bitmasks/#SPPIXMASK) · [Coadd masks and wavelengths](https://www.sdss4.org/dr17/spectro/spectro_basics/) · [Spectral products](https://www.sdss4.org/dr17/spectro/spectro_access/)

The manifest uses the `SPALL` solution without quasar templates (`*_NOQSO`) for objects whose main class is GALAXY and the default `SPALL` solution for other classes. Inspect `zwarning`, `redshift_column`, and plate quality before using a new target's redshift for a catalog comparison. For repeat observations, de-duplicate by sky position and inspect observation metadata. Two reductions of the same exposures provide a reproducibility check, while independent visits test astrophysical persistence. [SDSS redshift caveats](https://www.sdss4.org/dr17/spectro/caveats/) · [SDSS catalog selection](https://www.sdss4.org/dr17/spectro/catalogs/)

## Credit

The observational measurements are from the Sloan Digital Sky Survey. Cite the DR17 release and the survey's applicable technical papers, and include the survey acknowledgement with scientific outputs. Project software licensing applies to project code. Refer to SDSS for the public data policy. [Abdurro’uf et al. (2022), ApJS 259, 35](https://doi.org/10.3847/1538-4365/ac4414) · [SDSS citation guide](https://www.sdss4.org/collaboration/citing-sdss/) · [SDSS public data policy](https://www.sdss.org/collaboration/image-use-policy/)

## Bounded autonomous intake

`cohort plan` runs one public SkyServer query for one specified plate and MJD. It selects science-primary BOSS main-program galaxies with compatible pipeline versions, clean galaxy-only redshift flags, and redshifts between 0.02 and 0.10 inclusive. Both pipeline class fields must identify a galaxy. Input order is deterministic by fiber. The default batch contains up to 10 targets, with a hard maximum of 100. The selected IDs, exact query and hash of the returned CSV are frozen in a plan before any spectrum is downloaded. An empty selection is preserved as an empty cohort. Use a new plan filename for each selection. [DR17 catalog](https://www.sdss4.org/dr17/spectro/catalogs/) · [DR17 data access](https://www.sdss4.org/dr17/data_access/tools/)

```sh
python -m open_space.cohort plan --plate 6497 --mjd 56329 --limit 3 \
  --output work/first-plan.json
receipt=$(python -m open_space.cohort fetch --plan work/first-plan.json --output work/first-cohort)
python -m open_space enqueue --receipt "$receipt"
```

The live query above returned fibers 88, 96 and 100. The fetch step downloaded all three lite spectra and verified the galaxy selection against the FITS metadata. Its receipt lists the per-target inputs and hashes and links to the plan by SHA-256. Re-fetching into the same output directory checks source and derivative digests against existing records. Rerun after a transient archive error to complete an interrupted batch. To advance, pass `--after-fiber 100` with the same plate and MJD, writing to a new plan file. For parallel campaigns, assign distinct plate/MJD visits to workers. The sample is conditioned on one plate and a targeted spectroscopic survey. It cannot supply a sky-wide prevalence estimate.

A pinned target CSV ships with the repository. It lets a worker reconstruct this first cohort during a catalog outage as long as the spectrum archive is reachable. Use this mode for preselected IDs; the plan records their origin as an operator list. The converter verifies any supplied source hashes and the spectral selection on the downloaded files.

```sh
python -m open_space.cohort plan --plate 6497 --mjd 56329 --limit 3 \
  --input-csv data/cohort-6497-56329.csv --output work/pinned-plan.json
receipt=$(python -m open_space.cohort fetch --plan work/pinned-plan.json --output work/pinned-cohort)
python -m open_space enqueue --receipt "$receipt"
```

Other lists must contain `plate,mjd,fiber`, optionally `source_sha256`, and belong to the requested visit and fiber range. The command rejects duplicates and over-limit lists. Distinct fibers on the same plate can share sky or calibration residuals. Reobservations on independent visits or instruments are needed to test those effects. [Pinned target list](../data/cohort-6497-56329.csv)
