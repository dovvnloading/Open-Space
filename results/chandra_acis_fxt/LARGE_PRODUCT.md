# Recovering the group-1 observation

Activate these commands in the existing Chandra workflow only if the deposited group-1 record needs repair. ObsID 27099 has a single product, `acisf27099N002_evt2.fits.gz`, 127,158,960 bytes with SHA-256 `40405115dd4e6df3c7dc5d081ded4674c8a1631929cddeb882554b876a76fd24`. Its uncompressed FITS SHA-256 is `158372f2b47dd3eda7bdd67b40207ca1a70717b49f4ecd4c6a5863cda4db3113`.

Use the primary workflow's Python 3.11 patch, pinned requirements, and `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`. Run the two fixture modules in Actions before recovery. Here `$S` denotes `results/chandra_acis_fxt`, `$E` its absolute evidence directory, and `$T` a directory in `RUNNER_TEMP`.

```sh
(cd "$S" && python -m unittest -v test_scan_units test_large_obs)
python "$S/large_obs.py" phase "$E"
```

The original program checkpoints the full raw search before any catalog request. An absent raw scan calls for unit recovery. A present raw scan, numerous sealed response memos and incomplete row classifications identify interrupted follow-up. A fully complete baseline needs no work. The historical public group-1 attempt finished its product download but was cancelled before a raw search was saved. If the primary terminal artifact is genuinely lost, phase checks the immutable plan's original ledger, the terminal missing-artifact audit and an authenticated environment from a sibling worker of the same run.

**Only if the raw scan is absent**, fetch the exact event once into an Actions artifact, prepare the units, run the independent group workers and upload each output with `if: always()`. Use up to eight concurrent jobs. The original unit manifest states the exact group count. `G` selects one group. The event artifact can be restored in subsequent workers. For targeted retries, a work directory with valid existing unit capsules skips those completed calls.

```sh
python "$S/large_obs.py" fetch "$E" "$T/retained/27099"
EVT="$T/retained/27099/acisf27099N002_evt2.fits.gz"
python "$S/scan_units.py" plan "$E/observations/27099.json.gz" "$EVT" "$T/unit-plan.json" --groups 16
python "$S/scan_units.py" work "$T/unit-plan.json" "$EVT" "$T/units-group-$G" --group "$G"
python "$S/large_obs.py" prepare "$E" "$T/prep" --groups 8 --manifest "$T/unit-plan.json" --evt "$EVT" --capsules "$T/all-unit-artifacts"
```

Where a complete original raw search exists, prepare it directly.

```sh
python "$S/large_obs.py" prepare "$E" "$T/prep" --groups 8
```

Distribute catalog groups 0 through 7 with `$T/prep` unchanged. Each completed real or control row receives an atomic receipt. A rerun may restore an interrupted group's folder and finish its other rows. Each group shares identical cone responses within its own cache. The original survey classifications that already have query receipts stay under the standard corrective process. Historical URL-matched HTTP bodies are copied and authenticated from loose or packed original storage. New classifications require the ordinary corrected catalog function and an independent offline replay on the saved bodies. Keep partial responses, failures and the worker's zip digest when restoring artifacts. Start a fresh group output if a file upload ended mid-write, retaining the interrupted artifact for diagnosis.

```sh
python "$S/large_obs.py" catalogs "$E" "$T/prep" "$G" "$T/catalog-group-$G"
```

Use the ordinary terminal follow-up plan and the canonical worker output layout. In the corrective group owning 27099, reserve it for the following command and run all other assigned ObsIDs normally. Keep the same workflow run and attempt as that corrective worker so the finalizer can authenticate its namespace. Restore the prior raw and catalog artifacts under separate temporary directories. Set the usual `CIAO_PREFIX` and `RETAIN_EVENTS`.

```sh
python "$S/large_obs.py" finish "$E" "$T/prep" "$T/all-catalog-artifacts" "$T/followup-plan" "$T/worker-output" --retained "$T/retained"
```

This runs the usual visit processing and real peak calibration. Promotion requires exact EVT2 verification, the original raw-search identity, all 24 controls and 12 injection trials, every trigger’s required catalog evidence and all applicable real photometry/diagnostics. A missing or conflicting group, unsealed row, missing original source, changed bytes or partial calibration remains an explicit gap. The finalizer seals this worker beside the other repairs. It retains the original group ledger, raw unit plan and digests, all prefetch body and metadata files, per-row receipts and the revised science ledger. Keep the raw unit artifacts and temporary retained original event available at least through the completed final verification. Recovering any unit after artifact expiry requires recomputing that unit on the identical archived event bytes.
