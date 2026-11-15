# Corrected sky positions for the SMGPS supernova-remnant image index

The VizieR image table for Anderson et al. (2025) assigns displaced Galactic longitudes to the SARAO MeerKAT Galactic Plane Survey cutouts. The two coordinates labeled as the image center in the index differ from the center computed from the native FITS world coordinate system (WCS). `coordinates.csv` supplies the native center for each named file and identifies index positions that fall outside that same image. The affected product is `J/A+A/693/A247/list`. The original radio images and their native WCS remain usable. The corrected coordinates apply to file centers. The five separate source tables contain their own source-center measurements.

## Magnitude and cause

Across **693 of 693 files**, the median displacement is **88.776′** (1.480°). The largest is **179.976′**, almost three degrees. In **649 images (93.7%)** the listed center falls outside that image’s pixel boundaries. A five-arcminute cone centered on the native WCS center geometrically excludes its own indexed center for **673 of 693 files**. Vela Jr. (G266.2−1.2) is indexed 157.990′ away at longitude 268.81667°. Its native central longitude is 266.182932°.

The index coordinates agree with a simple empirical expression. With the FITS one-based image center at `x=(NAXIS1+1)/2` and `y=(NAXIS2+1)/2`, its longitude is `CRVAL1 − (x−CRPIX1) CDELT1` and its latitude is `CRVAL2 + (y−CRPIX2) CDELT2`. The native longitude uses the opposite sign for the pixel offset and the full spherical SIN projection. All 693 published positions agree with the reflected-linear expression to within 0.02545″ (median 0.01218″). This identifies the numerical transformation visible in the public data. The empirical calculation does not identify which software generated the index.

## Independent checks

- **Original FITS headers.** `headers.csv` contains the SARAO native WCS keywords, observed image sizes and listed center positions for every filename. The SHA-256 in each row covers the native FITS header cards from SIMPLE through END inclusive. The release uses diagonal, nonslanted Galactic SIN axes with zero rotation. The repair applies inverse orthographic projection at the geometric pixel center, and uses the forward native projection to test whether the listed point lies inside the pixel boundaries. `check_wcslib.py` independently computes every center with Astropy/WCSLIB. The maximum difference after CSV rounding is 0.000025″.
- **Image pixels and an independent telescope.** In six northern fields, we queried the NRAO VLA Sky Survey at the native WCS positions. We kept sources of at least 10 mJy with reported deconvolved major axis at most 45″ and centers at least 55 MeerKAT pixels from the image edges. Twelve sources passed. Seven lie within 6″ and ten within 12″ of a MeerKAT peak. Four translations of the comparison coordinates by 75″ along the pixel axes gave one match in 48 trials at 6″ and two in 48 at 12″. We used 8 robust-sigma local maxima and did not fit a registration shift. These checks locate the actual radio emission in the six images. They do not establish subarcsecond absolute astrometry throughout the survey. `nvss-markers.csv` retains every selected source and separation, including misses. `pixel-check.csv` records per-field counts and digests. One of the six fields had zero qualifying sources.
- **Position-based retrieval.** Live five-arcminute VizieR cones at the native centers of G043.942+01.442, G051.746+01.384 and G057.868+00.221 each returned zero image rows. Repeating those queries at the mislabeled centers returned exactly the named file in every case. The six queries and their exact URLs are in `cone-check.csv`. VizieR documentation explicitly labels GLON and GLAT as coordinates of the image center. A positional search of this index selects indexed centers; a filename lookup still identifies the intended file.

## Duplicate in the new-candidate table

Anderson et al. reported 237 new remnant candidates. Table A.4 records 152 and 153 both give **G354.353−00.018**, Galactic longitude 354.353°, latitude −0.018°, radius 2.3′ and reliability class II. We found 236 distinct names and 236 distinct full scientific records after excluding record number. The new-candidate subset of the image list also contains 236 distinct names, all matching the table's distinct names. Both identical lines occur in the authors’ original `MeerKAT_candidates.tex`. The released records thus support 236 individually listed candidates. This count cannot distinguish an inadvertent duplication from an intended candidate omitted before publication.

## Reproduce

The measurements were retrieved on 2026-11-15. Run commands from this directory. No network or third-party packages are needed for the primary replay.

```sh
python audit.py
```

To fetch the public inputs again, calculate the native positions, compare all coordinates with WCSLIB, repeat the sky queries, verify the six full raster images against NVSS, and check the candidate counts, use the following commands. The auxiliary scripts require Astropy, NumPy and SciPy. Downloaded ranges and images are cached locally. To force a fresh retrieval, run in a new directory containing the committed files only.

```sh
python harvest.py
python audit.py
python check_wcslib.py
python cone_check.py
python pixel_check.py
python count_candidates.py
```

`coordinates.csv` uses Galactic degrees, FITS one-based pixel coordinates and arcminutes for the great-circle index displacement. `listed_inside=0` flags a listed center outside its own image. `model_residual_arcsec` compares the empirical reflected-linear formula to the index. Use the native FITS WCS for arbitrary positions within an image.

## Source products

- Anderson et al., [Supernova remnant candidates discovered by the SARAO MeerKAT Galactic Plane Survey](https://doi.org/10.1051/0004-6361/202451038), including [the original LaTeX source](https://arxiv.org/src/2409.16607v1).
- SARAO, [original 1.3 GHz SMGPS release and FITS cutouts](https://doi.org/10.48479/0n8c-5q84). The authors describe each cutout as centered on its source.
- CDS, [VizieR catalog and image list](https://doi.org/10.26093/cds/vizier.36930247).
- Condon et al., [NRAO VLA Sky Survey](https://www.cv.nrao.edu/nvss/), an independent 1.4 GHz survey with a 45″ beam.
