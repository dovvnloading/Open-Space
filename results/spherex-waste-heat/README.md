# SPHEREx DR1 search for stellar infrared excess

No new stellar excess or observational null has been established. This draft records a targeted follow-up of published candidates and a separate blind search within 300 pc. Source confusion, photospheric structure, companions and dust are competing hypotheses for every excess.

## Completed

Verified the November 2026 DR1 all-sky cube release (102 spectral channels, 0.75–5 µm) and read the repository’s inclusion bar. [IRSA DR1](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/overview_dr1.html)

## Pending

- Obtain and validate epoch-propagated identities, existing follow-up and passive control positions.
- Measure DR1 spectra, inspect local backgrounds and measure WISE epoch and catalog fluxes. Assess whether excess originates on the moving star.
- Freeze the blind field, measure selection attrition, compare matched controls and inject 100–600 K excesses. Record missing coverage and measure sensitivity before interpreting non-detections.

The observing unit is a Gaia star at the observation epoch with an explicitly identified infrared counterpart. An offset mid-infrared centroid or a resolved background emitter can falsify the proposed association. A featureless low-S/N spectrum cannot establish a thermal emitter’s origin.

## Published contamination controls

Candidates D and E are known background-galaxy blends resolved in [JWST/MIRI imaging and spectroscopy](https://arxiv.org/abs/2607.09460); candidate G has a [background radio AGN](https://arxiv.org/abs/2501.05152). Treat them as known confusion controls, not new stellar-excess discoveries. Other candidate associations remain to be assessed individually.

`wise/fetch.py` requires an explicit `--output` campaign directory; cached response/receipt pairs must match the requested URL, query, byte count, and SHA-256. Partial, corrupt, or changed-query pairs and preserved errors stop the request; use a separate new campaign for changed inputs. A zero-row successful response remains a valid receipt.

For a new D: campaign, fetch the frozen response pairs and write derived tables to a separate new directory:

```powershell
python results/spherex-waste-heat/wise/fetch.py --output D:/Open-Space/campaigns/spherex-wise-20261117/archive
python results/spherex-waste-heat/wise/measure.py --archive D:/Open-Space/campaigns/spherex-wise-20261117/archive --output D:/Open-Space/campaigns/spherex-wise-20261117/derived
```

The NEOWISE single-exposure request is restricted to MJD 60291–60524; it is not a full 2010–2024 motion baseline.
