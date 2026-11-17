# SO/PHI-FDT: a held-out radial-proxy persistence screen

A fixed prediction for one magnetic-field proxy patch **failed** its two-scan holdout rule. The selected patch g42 increased by 8.659 G between two August 13, 2026 development scans. The frozen test required the same positive change, relative to its reference and control median, to exceed 3.236 G in each of two August 15 scans. The first held-out score was −9.846 G and failed; the second was +5.296 G and passed. Both were required, so the registered result is a non-replication.

This is a quality-limited comparison of `B_LOS/μ` medians, a **radial-field proxy** in gauss under an assumed radial field. It is not a measurement of true magnetic flux emergence or evidence for a new physical mechanism. Changing viewing direction, transverse fields, registration, solar evolution and inversion artifacts remain plausible explanations. No calibrated global false-alarm probability or population inference is claimed.

## Measurement and test

The observing unit is one SO/PHI full-disc telescope v03 L2 acquisition. Development used 2026-08-13 10:25 and 16:25 UTC, followed by held-out acquisitions on 2026-08-15 16:25 and 19:25 UTC. Each is a roughly 4½-minute scan, not an instantaneous image. Per-image acquisition timestamps and nonzero processing-operand IDs are disjoint; the scans share the instrument, calibration and inversion pipeline. The [input inventory](inventory.json) contains exact filenames, public archive URLs and SHA-256 hashes. `LAST_PRODUCT` URLs are mutable, so the hashes identify the bytes used.

The pre-pixel protocol defined a 5×5 grid. Each patch is a fixed 317-node solar-surface lattice with a 0.02-radian angular radius; all sampled nodes and bilinear stencils had to pass four-map visibility, μ ≥ 0.65, image-edge and moving-dust checks. This is a finite lattice, not proof that every point of a continuous disk is valid. Twenty-four patches were usable; g13 was excluded by the August 13 dust footprint. No patch had an invalid science stencil. Among the 24, the registered code selected the largest absolute development change, compared it with the median of the other patches, and required five times their robust MAD scale (with a 1 G floor) and unchanged sign under ±1-pixel registration shifts. These patches share calibration and can contain real solar evolution; their MAD is not a calibrated noise-only or false-alarm distribution.

| Quantity | g42 median B_LOS/μ (G) | Change from early development (G) |
| --- | ---: | ---: |
| Aug 13 10:25 development reference | −70.762 | 0 |
| Aug 13 16:25 development | −62.103 | +8.659 |
| Aug 15 16:25 holdout | −80.729 | −9.967 |
| Aug 15 19:25 holdout | −65.588 | +5.174 |

The development control median was −0.122 G and control MAD scale 1.079 G. The frozen signed rule subtracts that median, giving holdout scores −9.846 G and +5.296 G versus the 3.236 G threshold. The [prediction](evidence/prediction.json), [primary result](evidence/primary-result.json), and the full [development patch audit](evidence/patch-audit.csv) preserve the actual values and 25 dispositions.

The registered fixed-detector-coordinate control samples a different solar location as the Sun rotates. Its holdout changes were +99.595 and +83.789 G relative to its development reference; these are diagnostics, not a separate detection. In secondary ±1-pixel checks, the first held-out patch's change relative to the early reference ranged from −9.418 to −4.596 G; the second ranged from +3.803 to +12.708 G. Relative to the later development scan, the second held-out shifts ranged from −4.855 to +4.049 G, crossing zero. All 24 eligible patches were usable in both held-out scans. The 23 other-patch change MADs relative to the early reference were 1.687 and 1.367 G. These post-selection diagnostics do not rescue the failed primary rule or quantify independent solar noise. Full values and shifted measurements are in [secondary.json](evidence/secondary.json) and [secondary-patch-audit.csv](evidence/secondary-patch-audit.csv).

## Instrument and interpretation limits

The instrument team describes all four products as on-board averages of three datasets. Its exact quality-note row for the later development map reports rings and stripes in Bmag; `B_LOS` is derived from Bmag and inclination, so an effect on this proxy is possible but not quantified. The displayed note cells also repeat cumulative labels, making their attribution to individual exposures ambiguous. `ACCSTAT2=64` is present in all four files, but its meaning was not established. The [quality memo](evidence/instrument-quality-notes.md) preserves the row-attachment caveat, HTML source hashes and timing details; the [team release notes](https://www.mps.mpg.de/solar-physics/solar-orbiter-phi/data-releases) describe dust and other calibration artifacts.

During development the two-spacecraft viewing separation from Earth was close to quadrature. HMI and SO/PHI do not measure the same line-of-sight component at that geometry. No HMI, EUI or target-thumbnail pixels were used. `DATE_EAR` and `DATE_SUN` in the FDT headers are corrected *start* times; comparisons with other instruments must also account for scan duration, time scale and observer geometry. The result establishes neither persistence nor absence of physical evolution.

## Evidence and reproducibility

[PROTOCOL.md](PROTOCOL.md), [PACKET.json](PACKET.json), the geometry, selection and validation scripts, and [inventory.json](inventory.json) retain the exact pre-development bytes from commit `24ec48d083e4214d5e2cb3eb5c17eaf3ccf169b2`. The [development authorization](evidence/development-authorization.json) and [holdout authorization](evidence/holdout-authorization.json) were local records written before their respective analyses; they are not an externally registered preregistration or cryptographic timestamp. Original commands, logs and run receipts are in [evidence](evidence/). The [secondary script](evidence/secondary_diagnostics.py) was written after the primary result for predeclared diagnostics and did not choose the primary result or its threshold.

Before inspecting solar arrays, [synthetic geometry](evidence/predev-test_geometry.log), [selection and injection](evidence/predev-test_development.log), and [end-to-end pipeline](evidence/predev-test_pipeline.log) checks passed. In the simple Gaussian toy, 40 null trials had zero selections and 40 of 40 extended injections were recovered; these trials are implementation checks, not calibration on real solar backgrounds.

A separately coded [manual interpolation check](evidence/root-independent-holdout-review.json) reproduced the two held-out medians from the same pixels; it is a computation check, not an additional observation. Its [script](evidence/root-independent-holdout-check.py) records the local raw-data paths used.

The compact [g42 stencils](evidence/g42-stencils.npz) contain only the four source pixels per node and the fractional coordinates/μ needed to recompute each of the four reported g42 medians. The **post-holdout** [stencil exporter](posthoc_stencils.py) did not participate in selection or prediction. The offline check independently interpolates those pixels, recomputes the development selection, controls and frozen holdout rule, and verifies the packet hashes. It requires Python 3.14 and NumPy; it does not rederive image calibration, WCS or the full 25-patch pixel extraction:

```sh
python results/solar-phi-fdt-pilot/check_offline.py
```

Raw FITS remain external (about 58 MB for these four products). To replay the original event-level protocol, use a separate, clean checkout at exact commit `24ec48d083e4214d5e2cb3eb5c17eaf3ccf169b2`, the input hashes above, and Python 3.14 with the exact versions in PACKET.json. Set `PYTHONPATH` to the SunPy environment; use the recorded [development command](evidence/run-development.ps1), then compare its output with the deposited development record. The frozen validation code checks that Git HEAD and the packet match; it expects a separately frozen prediction and local authorization before reading held-out arrays. Use a **new** output path and preserve the original receipts. SunPy/Astropy coordinate transforms and native standard-L2 inversion results are dependencies of this analysis, not independent validation of them.

These are external FITS inputs, not jobs in the Open-Space worker ledger; no ledger job or report ID is claimed. Original data come from the public [Solar Orbiter Archive](https://soar.esac.esa.int/soar/) and software dependencies retain their respective licenses.
