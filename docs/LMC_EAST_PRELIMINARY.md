# LMC east: provisional grid extension and source-attribution follow-up

**Evidence status: provisional, not a published independently reproducible result.**
This page transcribes the reported outputs of an external extension and follow-up
run, as recorded on 2026-11-13. The raw science inputs, measurements, detailed
archive receipts and catalog-query responses are **not in this repository**.
The report has passed arithmetic/schema checks only; this checkout cannot yet
independently verify the observational findings. Do not cite this page as a new
discovery or a survey-wide false-positive rate. The original, separately
published six-field experiment remains documented in [FFI.md](FFI.md).

A structured transcription, including the exact limitations and filenames of the
unavailable artifacts, is in
[`results/ffi-lmc-east-preliminary.json`](../results/ffi-lmc-east-preliminary.json).

## Design and reported outcomes

The reported extension followed the existing [coordinate-only runbook](DISCOVERY.md#run-a-blind-image-survey-on-your-compute): a frozen 3-column by 2-row grid
at ICRS right ascensions **78.2°, 78.6°, 79.0°** and declinations **−69.4°, −69.2°**.
Each field requested 21 × 21 TESS pixels in Sectors 1 and 2. The code searches
every eligible interior pixel rather than only stars in a selected target list.

| Reported search result | Count |
| --- | ---: |
| First-sector spatial peaks | 23 |
| Two-sector recurrent sources | 13 |
| Image-artifact vetoes | 6 |
| No matching second-sector spatial peak | 3 |
| Period inconsistency | 1 |
| Recurrent sources repeating in at least one later sector | 12 of 13 |
| Automatic VSX associations | 9 of 13 |
| OGLE period-compatible associations | 11 of 13 |
| Direct Gaia period matches | 8 of 13 |

The reported usable common sky area was **0.04991 square degrees**, and the 13
recurring positions were all at least **104.6 arcseconds** from preselected TESS
two-minute target centers. These are reported measurements, not independently
recomputed from WCS and archive inputs in this commit. Adjacent tiles can overlap;
raw candidate counts and nominal coordinate footprints alone do not provide
unique-source surface densities.

### Reported controls and their interpretation

| Injected peak semiamplitude | Recovered / injected |
| --- | ---: |
| 1.0 electron/s | 13/18 |
| 2.5 electrons/s | 18/18 |
| 5.0 electrons/s | 18/18 |

The reported null experiments used 96 field trials each and 0.799
square-degree-trials of effective area. Independent whole-frame time shuffles
produced **0** in-area recurring peaks; half-day-block time shuffles produced
**13**, or about **16.27 per square-degree-trial** under *that particular null*.
These null procedures destroy different amounts of astrophysical and
instrumental time correlation. They do **not** establish the false-positive
fraction among the 13 candidates. Injections also do **not** establish full
survey completeness: amplitude, crowding, time sampling and sky region matter.

## Source attribution: what still needs observing

| Reported source case | Independent evidence and caveat | New-sector blind recovery | Frozen-frequency phase check |
| --- | --- | ---: | ---: |
| **OGLE-LMC-CEP-1429**, near (79.02466°, −69.34685°) | Catalog period 3.599551 d; third harmonic 1.199850 d is compatible with the ~1.196 d image-grid signal. Reported detections in OGLE and both EROS colors. An eclipsing binary 11.55″ away and a Gaia star only 1.06″ from the Cepheid complicate attribution. | 2/38 | 24/26 |
| **OGLE-LMC-ECL-11330** | Catalog period 7.284297 d; fourth harmonic 1.821074 d is compatible with the ~1.822 d image-grid signal. Reportedly supported by resolved two-band EROS data; useful as a known-binary control, not a new system. | 26/37 | 25/25 |

The reported survey examined **87 archive entries** across these two positions;
entries with unusable detector coverage must remain visible in the denominator.
The blind search and the targeted phase test answer different questions. A
frozen-frequency signal with consistent phase may be real even when it does not
pass the full blind spatial/period selection. It is not independent evidence of
which of two stars 1.06″ apart is modulating. **The strongest remaining
observation is resolved time-series photometry across several 1.20-day cycles**
covering the Cepheid/Gaia close pair, with cadence chosen to distinguish daily
aliases and the fundamental 3.599551-day cycle.

## Rerun the grid from the published code

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[archive]'
export OPENBLAS_NUM_THREADS=1
python -m examples.image_survey plan --origin 78.2 -69.4 --shape 3 2 \
  --step 0.4 0.2 --sectors 1 2 --limit 6 \
  --output runs/lmc-east/selection.json
python -m examples.image_survey pin runs/lmc-east/selection.json \
  --cache runs/lmc-east/ffi --output runs/lmc-east/pinned.json
python -m open_space.discovery --root runs/lmc-east/queue init \
  --ffi-plan runs/lmc-east/pinned.json
python -m open_space.discovery --root runs/lmc-east/queue work \
  --ffi-cache runs/lmc-east/ffi --workers 4 --max-jobs 300
python -m open_space.discovery --root runs/lmc-east/queue export runs/lmc-east/export
python -m examples.image_survey triage runs/lmc-east/export \
  --output runs/lmc-east/followups.json
```

The commands are an **independent rerun**, not a claim that the provisional counts
above will necessarily reproduce. Pin the newly returned scientific digests and
freeze code and catalog versions before comparing. A failed archive request is
not an empty search result. For later-sector analysis, construct a frozen plan
using the schema of [`data/ffi-source-revisit.json`](../data/ffi-source-revisit.json)
and run [`examples/revisit_ffi_source.py`](../examples/revisit_ffi_source.py);
record off-detector cases instead of removing them silently.

## Deposit and audit the external evidence

The prior run produced separately stored archives and tables, including
`LMC-East-Evidence.zip`, `LMC-East-Image-Inputs.zip`,
`LMC-Source-Evidence.zip`, `LMC-Updated-Followups.csv` and
`LMC-Source-Sector-Register.csv`. Their **reported sizes** are in the JSON.
Their SHA-256 hashes were not available for this commit. The original FITS
archives exceed GitHub's 100 MB regular-file limit; preserve the inputs in a
stable external scientific deposit (or release storage that supports the size),
with public URLs, hashes, licensing and archive citations. Do not use Git LFS
alone as a substitute for a public, verifiable data release.

Run the offline metadata audit from a checkout:

```sh
python -m open_space.ffi_release_audit \
  --summary results/ffi-lmc-east-preliminary.json \
  --output runs/lmc-east/audit.json
# After downloading every external artifact into one directory:
python -m open_space.ffi_release_audit \
  --summary results/ffi-lmc-east-preliminary.json \
  --evidence-dir /path/to/downloaded/artifacts --check-zip-crc --require-all \
  --output runs/lmc-east/audit-with-files.json
```

**For publication**, an independent researcher should verify the source image
science hashes, replay the detection and control programs with exactly frozen
inputs, reconcile every first-sector peak with a second-sector disposition,
recompute coverage from the WCS, cross-check catalogs with preserved responses,
and compare source-attribution positions across sectors. Attach these outputs,
their SHA-256 hashes and a machine-readable task register before changing the
JSON `evidence_state` or moving these findings into
[`DISCOVERY_RESULTS.md`](DISCOVERY_RESULTS.md). Passing the offline summary audit
alone never promotes the claims.
