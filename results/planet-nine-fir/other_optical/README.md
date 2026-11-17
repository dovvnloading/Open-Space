# Optical tests of the two Chen candidates and the Rowan-Robinson candidate

**No linked optical point source survived on the tested trajectories.** Each Chen direction comes from one far-infrared observing day. We searched 792 discrete circular trajectories for each direction, with 400–900 au in 50 au steps and unknown heading sampled every 5°. The Rowan-Robinson model has 504 trajectories, with 210–240 au in 5 au steps and the same heading increments. The Rowan-Robinson check below covers one northern sector of its much larger annulus. There is no posterior or measured exclusion fraction for the continuous allowed families. [Candidate definitions and ephemerides](../other/README.md)

## DECam

| Published target | Usable at ≥2 r epochs | Epochs UTC | Injected at r=22.5 | Injected at r=23 | Injected at r=23.5 |
|---|---:|---|---:|---:|---:|
| Chen 1, July 2007 | 679/792 | 2014-11-20, 2016-11-03, 2025-02-15 | 244/280 | 228/280 | 159/280 |
| Chen 2, July 2006 | 622/792 | 2024-01-19, 01-22, 01-28 | 219/280 | 168/280 | 23/280 |

The Chen 1 exposures lasted 90, 90 and 15 s. Each Chen 2 exposure lasted 70 s. All trial positions use the actual science-header midtimes, full DECam focal-plane CCD astrometry, and the epoch-specific seasonal parallax. At every usable direction the code fit a local empirical PSF, subtracted a local annular sky estimate, and combined inverse-variance-weighted fluxes across independent visits along the trial trajectory. The noise scale came from displaced blank-sky measurements on the same CCD. The search required a combined score of 5.5 and at least two individual scores of 2.5. The numbers above count synthetic planets passing the initial flux-and-link test. The 280 trajectories were sampled without replacement from each 792-row grid before masks. Final confirmation and identity vetting have a different selection function.

Each CCD used point-source calibrators and photometric scale measurements from Legacy DR10 r. The pipeline rejected inadequately calibrated CCDs and used the original data quality masks. Post-2022 visits also excluded positions within 2.2″ of a significant DR10 counterpart. PSF calibration does not assume that a Gaussian describes the camera. The injections use empirical point sources, subpixel coordinates and an approximate gain of 3.8 e⁻/ADU. They add planets to local science stamps and repeat the flux measurement with fixed calibrations and existing masks. Planet self-contamination of the calibration and saturation were not tested. [Per-CCD calibrations for Chen 1](Chen1_calibration.csv) · [Chen 2 calibrations](Chen2_calibration.csv) · [Flux-stack procedure](chen_pixels.py)

Each field yielded one initial trigger.

- **Chen 1, track 162.** Its three positional scores were 2.91 (2025), −0.17 (2016) and 7.45 (2014). At the 2014 sky position, a separate visit in 2016 measured 0.857 nanomaggies when the trial planet should have moved, versus 0.708 nanomaggies in 2014. The 2025 frame has only one marginal residual at the forecast location. The sequences fail as evidence for a moving counterpart. [Off-date photometry](Chen1_stationarity.json) · [Fixed-position image comparison](Chen1_162_vet.jpg)
- **Chen 2, track 412.** All three positive measurements fall along the same extended, edge-on galaxy in the original January 2024 frames and in the independent DR10 r coadd. It produces detectable flux at all three trial positions on all three nights. The nearest-catalogue-center veto misses its extended light. No moving point source is identified. [Nine crossed epoch/position measurements](Chen2_stationarity.json) · [Registered images](Chen2_412_vet.jpg)

Neither trigger passed this static-sky check. No unresolved source was advanced to a minor-planet identification. The test samples exact model directions. It does not search the spaces between all headings, the covariance of the single-day centroids, eccentric orbits, or the complete DECam archive. [Chen 1 measurements](Chen1_result.json) · [Chen 2 measurements](Chen2_result.json) · [Selected archive products](selected.json)

## ZTF difference images

We used the survey-provided science-minus-reference difference FITS, difference PSF and matching science mask. The first-stage statistic and two-visit requirement match the DECam test. An empirical noise estimate from up to 1,000 local random stamps accounts for background and matched-filter correlations. Synthetic sources were added in the difference-image plane at every tested eligible track, with the archived PSF and header `MAGZP`. These are ZTF natural r-band estimates at the zero-color convention. This procedure does not test reference self-subtraction, source-dependent subtraction, source Poisson noise, nonzero color terms, or any contribution from a pixel-level injection into the original reference construction. [IRSA ZTF image products](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html) · [ZTF data processing](https://arxiv.org/abs/1902.01872)

| Target | Usable difference-image visits | Trajectories eligible out of all sampled | Linked triggers | Recovered at r=21 | Recovered at r=21.5 |
|---|---:|---:|---:|---:|---:|
| Chen 1 | 2, September 2026 | 236/792 | 0 | 217/792 | 34/792 |
| Chen 2 | 3, 2025 | 260/792 | 0 | 252/792 | 80/792 |
| Rowan-Robinson | 3, August 2025 | 11/504 | 0 | 11/504 | 4/504 |

The 11 Rowan-Robinson directions are a small part of the annulus in a single ZTF quadrant. Center uncertainty at the infrared epoch, including the approximately 0.5° quoted by Rowan-Robinson, could move a trajectory far outside these tests. Some requested IRSA files returned temporary service errors. The table includes only sets in which the difference image, difference PSF and matching mask were retrieved and verified. [Per-visit IDs, dates and injections](ztf_results.json) · [Selected ZTF products](ztf_selection.json)

## Reproduction

The `other/` and `orbits/` siblings provide the literature inputs and ephemeris routines. Python 3.12 with NumPy, SciPy, Astropy, PyERFA and Matplotlib is sufficient. The frozen receipts include the archive parameters, returned files, timestamps and SHA-256 hashes. The downloader fetches only the data used in these reductions. A live TAP query can return rows in a different order or reflect a later database revision. [Input receipts](inputs.json)

```bash
python download.py
python chen_pixels.py Chen1
python chen_pixels.py Chen2
python ztf_tracks.py
```

Review of source identity used the original public images and cross-epoch sky measurements. `vet.py` repeats the crossed photometry for both published flags. The `query.py`, `calibration.py`, `static_sources.py`, `select.py`, `ztf_meta.py` and `ztf_pixels.py` scripts give the wider archive discovery and selection queries. The two survey tables overlap in the hypotheses tested and cannot be multiplied to estimate a combined exclusion.
