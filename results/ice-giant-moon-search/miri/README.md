# Neptune MIRI two-visit search

No two-visit source passed the specified search in either F2100W or F2550W. On 45,180 trial orbits per filter, 28,648 F2100W tracks and 26,444 F2550W tracks met the coverage and data-quality criteria. Sign-inverting the science images produced zero passing tracks in either band. The strongest positive scores reached 0.639 and 0.555 of the full selection threshold, respectively.

Of 126 public MIRI calibrated exposures, 11 put Neptune on the projected detector and none puts Uranus on it. All 36 Uranus frames are parallel images from MRS apertures. Their nearest geometric boundary lies at least 17.20 arcsec from Uranus. None covers an equatorial circular orbit at the tested radii of 50,000, 70,000, 100,000 or 150,000 km. At 250,000 km, at most 9% of sampled phases enter any one frame. The catalog instrument and target alone would give a substantially different impression of inner-moon coverage.

| Planet / proposal | Exposures | On-detector planet | Geometry |
|---|---:|---:|---|
| Uranus 1248 | 36 | 0 | Three MRS-parallel visits |
| Neptune 1249 | 27 | 0 | MRS parallels and offset fields |
| Neptune 1272 | 36 | 8 | Two direct Triton visits contain Neptune, four exposures each |
| Neptune 1604 | 27 | 3 | Three F1000W frames put Neptune about 14 pixels inside the right edge |

For every exposure, the input list joins the frozen MAST group and file manifests. The reduction reads and hashes the real primary and SCI headers and queries Horizons from JWST at the FITS mid-time. It maps all four sensor edges through the native SIP model at 129 points per edge and tests the planet center against this sky polygon. Each radius has 720 equally spaced phases in the equatorial plane. Fractions refer to the projected **detector rectangle**. They are upper bounds on usable phase coverage before saturation, bad pixels, planet and rings enter the mask.

GO 1272 imaged Triton directly on December 2 and 5, 2022. In each 13.875-second F2100W exposure the known satellite shows positive residual flux within a 4-pixel aperture at its JWST Horizons position. The two visits differ in structured background. The F2550W control fails in the second visit. The measurements fit a spatially quadratic background in a clipped 9–20 pixel annulus, reject nonzero DQ and unknown/nonpositive ERR, and omit neighbors within 12 pixels of Neptune. Fluxes are background-subtracted finite-aperture sums without aperture correction. Local contrast at the predicted pixel uses PSF-scale smoothing and a robust scatter estimate.

| Filter and visit | Two aperture sums (mJy) | Two predicted-pixel local contrasts |
|---|---:|---:|
| F2100W, December 2 | 0.085, 0.083 | 4.76, 5.0 |
| F2100W, December 5 | 0.062, 0.060 | 1.6, 1.8 |
| F2550W, December 2 | 0.197, 0.188 | 2.1, 1.8 |
| F2550W, December 5 | −0.122, −0.135 | −1.1, −1.1 |

Recompute the footprint and targeted satellite checks from the parent result directory.

```sh
python miri_audit.py manifests /data/miri-header-and-pixel-ranges /data/horizons miri
```

[Exposure-by-exposure geometric support](footprints.csv) · [Native pixel Triton measurements](triton_control.csv) · [Source FITS byte-range hashes](byte_ranges.csv) · [MAST JWST program 1272](https://www.stsci.edu/jwst/science-execution/program-information?id=1272)


## Orbit search and injection depth

The search uses the four calibrated exposures in each band, two dithers per independent visit, separated by about 3.09 days. A native 431-pixel field around the tracked Triton direction includes Neptune in every exposure. Each cutout retains its detector SCI, ERR, DQ and shifted SIP WCS. Explicit Horizons Neptune centers anchor a 301×301 planet-centered grid at 0.11 arcsec per pixel. The procedure masks planet glare, projected rings and all 16 cataloged Neptunian satellites at nine samples during each exposure. Only fully supported pixel rectangles with all nonnegligible bilinear contributors valid enter the images. DQ values of zero and corrected jumps are admitted. Gaussian PSFs use sigma 2.8 and 3.3 pixels. These are approximations to wavelength-dependent MIRI PSFs.

Circular equatorial templates include the Neptune J2/J4 mean motion, with 50,000–300,000 km radii at 1,000 km intervals and 2° spacing in phase. The scoring requires at least 5 per visit, 7 pooled and 2 after dropping either dither. Score units depend on local residual normalization. The 7-pixel planet/ring/known-source padding and 90% PSF support gate apply in both bands. This search does not quantify irregular, eccentric or inclined orbit populations.

The depth test planted a positive PSF on each of 32 off-grid circular orbits drawn uniformly in phase and in semimajor axis from 180,000–182,000 km. One synthetic moon at a time entered all four real registered images before background fitting and scoring. The same 32 positions were tested at every flux. Recovery requires a passing template within four registered pixels of the planted location in every exposure of both visits. Missing/low-support orbits count as misses. Each level uses the 540-template search in that annulus. This is registered-image, Gaussian-PSF completeness for that population. Detector calibration uncertainty, real MIRI PSF wings and broader-domain trial multiplicity remain outside the injection test.

| Filter | Injected flux density (mJy) | AB magnitude | Recovered / 32 |
|---|---:|---:|---:|
| F2100W | 0.2 | 18.15 | 0 |
| F2100W | 0.6 | 16.95 | 11 |
| F2100W | 2 | 15.65 | 24 |
| F2100W | 6 | 14.45 | 27 |
| F2550W | 2 | 15.65 | 0 |
| F2550W | 6 | 14.45 | 14 |
| F2550W | 20 | 13.15 | 22 |
| F2550W | 60 | 11.95 | 23 |

Neither filter reached 90% observed completeness. The measured fraction crosses 50% between 0.6 and 2 mJy in F2100W, and between 6 and 20 mJy in F2550W. Sampling uncertainty is appreciable with 32 plants; exact binomial intervals accompany the individual numbers. The measured depths imply no specific moon diameter without a mid-infrared thermal model. Fluxes use the instrument's calibrated flux-density convention and the 3631 Jy AB zero point.

From the parent result directory, use `search_inputs.py` to retrieve the eight calibrated-image regions and generate both orbit-search cubes. The broad scan stores a compact checkpoint per radius; the committed summaries give the tested and covered template counts per radius.

```sh
python miri/search_inputs.py manifests/products.csv miri/footprints.csv /data/miri-science /data/horizons miri
for f in F2100W F2550W; do
    python broad_scan.py /data/miri-science/${f}.npz /data/miri-science/${f}-broad --grid 50000 300000 1000 2
    python miri/plants.py ${f} /data/miri-science miri/${f}_plants
done
```

[Cutout and original data hashes](cutouts.csv) · [Native FITS byte receipts](cutout_ranges.csv) · [Ephemeris sources](horizons_sources.csv) · [F2100W support and thresholds](F2100W_radii.csv) · [F2550W support and thresholds](F2550W_radii.csv) · [F2100W depth and intervals](F2100W_plants/depth.json) · [F2550W depth and intervals](F2550W_plants/depth.json)

The scan receipts bind the historical `scan_sources/prepare.py` and `scan_sources/broad_scan.py` bytes used for the recorded MIRI scans. The parent-directory scripts were subsequently extended for other observations; rerunning them is a separate reproduction check.
