# Running an autonomous archive survey

Start a campaign from a checkout on Linux or macOS with Python 3.10 or later. The measurements use local CPUs and RAM. An optional Chat Completions-compatible model server supplies the research analyst and critic. Keep the SQLite queue on a local filesystem on each worker host. Public archives require outbound HTTPS. The largest reusable local catalog is the TESS rotation catalog and its exact-TIC index. Allow several GB for a pilot with pixel follow-up. Pixel files can exceed 100 MB per target/sector. Scale archive access gradually and reuse downloaded products across workers.

```sh
git clone https://github.com/dovvnloading/Open-Space.git
cd Open-Space
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[archive]'
```

## Advance to unexamined targets

The first three pages published in this release cover ranks 0–151 of the search order. This command starts at the next rank, downloads the corresponding NASA products, records checksums, runs two independent searches per target, schedules pixel and catalog follow-up, exports results and advances to the next page. The page size below is a computation budget. Setting it to another positive number in a new survey preserves the ranked order. An operator can pass many pages in a single invocation.

```sh
python -m open_space.survey --root runs/tess-s1-s2 --sectors 1 2 \
  --start 152 --page-size 32 --pages 1 --workers 4
```

Repeat the command to continue. The driver saves the original sector inventories under its root and keeps every target in its denominator. It stops at the requested page budget, inventory exhaustion, an incomplete page or a failed measurement. A job budget that ends partway through a page yields `in_progress`. The next invocation resumes it. If the result says `needs_retry`, inspect the failed task and its recorded error. Correct a transient archive, disk or model configuration problem, then explicitly retry that task. The page cursor advances after all required measurements and export have completed.

```sh
python -m open_space.discovery --root runs/tess-s1-s2/pages/00000152 status
python -m open_space.discovery --root runs/tess-s1-s2/pages/00000152 show TASK_ID
python -m open_space.discovery --root runs/tess-s1-s2/pages/00000152 retry TASK_ID
```

## Attach a model server and compute workers

Configure the identifiers served by your inference endpoint, and use a new campaign root when switching model setup. A single GPU server can serve both roles, or the same API gateway can route the two model identifiers to different backends. Supply a credential using `OPEN_SPACE_API_KEY` if your endpoint requires it.

```sh
export OPEN_SPACE_MODEL_BASE_URL=http://127.0.0.1:8000/v1
export OPEN_SPACE_ANALYST_MODEL=your-served-model
export OPEN_SPACE_REVIEWER_MODEL=your-review-model
export OPEN_SPACE_MAX_TOKENS=1600
python -m open_space.survey --models --root runs/model-survey --sectors 1 2 \
  --start 152 --page-size 32 --pages 8 --workers 4 --model-call-ceiling 300
```

The analyst reads measurements and selects among available observation-specific actions. Workers perform those actions on the archived data. The critic receives the original evidence, the analyst's decision, the follow-up results and any missing checks. An offset in a difference image can request a second catalog search at the measured position. The system records exact prompts, responses as parsed decisions, model IDs, request and response digests, reported usage and explicit failures. A response selecting an unavailable action is rejected. The ceiling counts *reserved* calls, including interrupted or unsuccessful attempts. Once exhausted, subsequent planners execute available standard checks and an outstanding critic records the missing review. The ceiling on the survey driver spans its pages. A single standalone discovery queue has its own ceiling. Keep the ceilings of independent shards within the total you want to spend.

## Reproduce the multi-archive campaign

```sh
python -m open_space.discovery --root runs/replay init \
  --tess-plan data/tess-plan.json --tess-plan data/tess-expansion-plan.json \
  --wise-plan data/wise-fields.json --ffi-plan data/ffi-fields.json \
  --ffi-plan data/ffi-lmc-grid.json --catalog data/discovery-planet-catalog.json
python -m open_space.discovery --root runs/replay work --workers 4 --max-jobs 700
python -m open_space.discovery --root runs/replay status
python -m open_space.discovery --root runs/replay export runs/replay/export
```

A worker claims a job under an expiring lease and renews it while executing. Completing a job publishes the result and schedules its children in one transaction. Duplicate proposals for identical work refer to one measurement. Lost owners cannot overwrite a newer result. The code and measurement-library versions are frozen for a running campaign. Inventory hashes, archive product identifiers, individual checksums, catalog versions, original job payloads and parent IDs travel with the export. A live catalog may change after the published snapshot. Compare the recorded query and response hash when reproducing a catalog match.

For a compute cluster, give every array element a different root and pass the same plans with `init --shards N --shard I`, where I runs from 0 through N−1. That partitions targets and maps both optical searches for a star to one node. Give the node-local `work` invocation its own `--workers` and `--max-jobs`. Use `--tess-cache`, `--pixels-cache`, `--ffi-cache`, `--wise-cache` and `--catalog-cache` to select persistent disks. Export every shard separately. For large batches use `OPENBLAS_NUM_THREADS=1` or the equivalent setting for your BLAS implementation to avoid multiplying parallel CPU threads. Pick the worker count according to available memory and archive throughput. The code has been exercised on the published sample. Cluster-wide throughput and optimal worker count depend on the operators' hardware and archive limits.

## Interpreting a candidate

| Measurement | What promoted it | What can still explain it |
| --- | --- | --- |
| TESS transit | A sector-1 box search predicts depths and individual events in another sector at fixed period, epoch and duration. A sign-inverted search, odd/even comparison, secondary search and held-out profile comparison challenge the signal. | Rotational variability, a neighboring eclipsing system, depth dilution, an aliased period, a secondary below the numerical veto. |
| TESS stellar variability | A single-frequency fit explains a measurable fraction of variance with comparable phase, frequency and amplitude in another sector. Whole-block shuffles are recorded for coherent cases. | Other harmonics, a star elsewhere in the aperture, evolving spots, short baseline, already published rotation. |
| TESS full-frame pixels | Independent temporal spectra and spatial peaks in two sky-projected cutouts agree on position and frequency after instrument checks. | Crowding, two stars in one pixel, spatially correlated residuals, known Cepheids, or confusing harmonics. |
| NEOWISE motion | Multi-epoch positions predict late observations better than a stationary model. Gaia crossmatching follows the independent fit. | Cataloged motion, centroid shifts in a blend, mismatched seasonal sources, a motion smaller than the threshold. |

The output `registry.csv` is a machine screening ledger. A flag remains visible after a catalog or image explains it. To evaluate a claim, open its `evidence` member in the export, then inspect the follow-up jobs in `registry.json` via their parent IDs. A failed catalog query is different from a successful empty query. A catalog miss alone gives no evidence of a new planet or star. A strong candidate dossier needs source localization, held-out data, tests of aliases and instrumental effects, relevant catalog and literature history, and a discriminating observation still to obtain. [Public results](DISCOVERY_RESULTS.md) · [Transit selection](TESS_PROTOCOL.md) · [Imaging](PIXELS.md) · [Full-frame search](FFI.md) · [Astrometry](WISE.md)

## Run a blind image survey on your compute

The image survey accepts sky coordinates without a target-star list. Generate a grid, download both sector cutouts for each cell and pin the scientific arrays, then run the normal worker. This example covers six new grid centers east of the published grid. Rows advance in declination and columns in right ascension, with steps specified in degrees. Across a larger grid, use `--start` and `--limit` to partition scan pages. Overlap depends on declination, detector orientation and requested cutout size. The operator chooses fields observed in both sectors, and an unavailable cutout is recorded as an acquisition failure.

```sh
python -m examples.image_survey plan --origin 78.2 -69.4 --shape 3 2 \
  --step 0.4 0.2 --sectors 1 2 --limit 6 --output runs/lmc-east/selection.json
python -m examples.image_survey pin runs/lmc-east/selection.json \
  --cache runs/lmc-east/ffi --output runs/lmc-east/pinned.json
python -m open_space.discovery --root runs/lmc-east/queue init --ffi-plan runs/lmc-east/pinned.json
python -m open_space.discovery --root runs/lmc-east/queue work \
  --ffi-cache runs/lmc-east/ffi --workers 4 --max-jobs 300
python -m open_space.discovery --root runs/lmc-east/queue export runs/lmc-east/export
python -m examples.image_survey triage runs/lmc-east/export --output runs/lmc-east/followups.json
```

The last command produces a machine-readable list of source positions, both measured periods, matched variable stars, original task IDs and remaining follow-up needs. Its automatic rule tests the VSX position within 15 arcseconds and direct, half or double periods in both sectors. An incomplete archive lookup receives its own state. The release replay associated nine of 14 recurring locations across the development and six independent fields with VSX by this rule. Five were referred to later images and resolved photometry. The fuller OGLE and source-level analyses recovered additional associations. For a prioritized source, run the later-sector experiment with a frozen position and frequency; its example plan also includes a nearby known star as a control. [Automatic triage](../results/ffi-automatic-triage.json) · [Five-sector example plan](../data/ffi-source-revisit.json) · [Follow-up measurement script](../examples/revisit_ffi_source.py)

## Other archive entry points

The older SDSS program retrieves spectra and searches recurring wavelengths. The PDS and HEASARC adapters retrieve observational inventories for future analyses, including the data-product URLs exposed by PDS. The released inventory includes three OSIRIS-REx observational products at Bennu and three NuSTAR pointings toward the Crab. Archive enumeration is the completed step for these two services. [Observation receipts](../results/archive-inventory.json) · [SDSS workflow](CAMPAIGNS.md)

```sh
python -m open_space.archive_index --limit 3 pds urn:nasa:pds:context:target:asteroid.101955_bennu
python -m open_space.archive_index --limit 3 heasarc 83.6331 22.0145 --radius 0.25
```
