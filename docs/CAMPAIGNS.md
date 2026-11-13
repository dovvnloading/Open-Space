# Run a research campaign

A campaign follows one SDSS DR17 plate visit through a finite sequence of fiber IDs. For each page it queries the catalog, freezes the selected targets, retrieves and checks the spectra, runs the research workers and exports one evidence packet per galaxy. The next invocation continues after the last fully exported page. Use a POSIX host with a local filesystem supporting file locks and SQLite WAL. The public archive commands require internet access and Astropy.

```sh
python -m pip install -e '.[archive]'
python -m open_space.campaign --root runs/boss-6497 --plate 6497 --mjd 56329 --limit 3
```

This command runs one page with deterministic analysis and zero model requests. It emits a JSON summary including the last committed fiber, completed jobs, model requests reserved, and whether the visit is exhausted. Running the command again advances the same visit. A resumed invocation reuses completed evidence and any archived inputs from its current page. If an archive download stopped partway through a batch, the downloader refetches and verifies the completed inputs before publishing the receipt. [Target selection and data semantics](DATA.md)

## Model-led and scheduled operation

Configure the analyst and reviewer as described in the model setup, then give the model-enabled campaign its own directory. Use the same flags and model settings every time that directory is resumed. [Model setup](MODELS.md)

```sh
python -m open_space.campaign --root runs/boss-6497-models --plate 6497 --mjd 56329 \
  --models --limit 3 --pages 2 --max-calls 8
```

Here each invocation processes up to two pages of three observations and reserves at most eight provider calls. The worker can use one analyst and one reviewer call per spectrum with detected features. If the call budget expires during a job, the checkpoint records the completed stages. The next invocation starts with a fresh per-invocation budget and resumes that job. For a campaign-wide spending cap shared across invocations, configure a quota on the model provider or limit the number of submitted invocations. The output-token setting, model names, endpoint, code fingerprint, plate/MJD and batch size are pinned for the campaign. Change the model mode, selection or code by starting a separate directory. The page and call limits are per-invocation controls and can be changed on resume.

The process exits at its configured page or call budget, a failed job, or visit exhaustion. A scheduler can invoke it periodically. This example runs at 15 minutes past each hour using a prepared virtual environment in a checkout at `/srv/Open-Space`. Use absolute paths appropriate to your machine. Add the model configuration to the scheduler environment when enabling models.

```cron
15 * * * * cd /srv/Open-Space && .venv/bin/python -m open_space.campaign --root /srv/runs/boss-6497 --plate 6497 --mjd 56329 --limit 3 --pages 1
```

Each new page performs one bounded catalog request and up to `--limit` sequential spectrum downloads. The catalog and model transports have 60-second network timeouts and each FITS download has a 30-second timeout. Set a wall-clock limit in your scheduler if needed. A POSIX file lock prevents two coordinators from advancing the same campaign directory concurrently. Use separate directories for independent visits, with per-visit resource budgets. A repeated invocation after `exhausted` returns without archive or provider requests.

## Inspect and resume

| Output | Meaning |
| --- | --- |
| `campaign.json` | The exact software version, method configuration and visit selection. |
| `state.json` | Current cursor, exhaustion state and campaign hash, written atomically after the entire page completes. |
| `pages/after-NNNN/plan.json` | Query, catalog response digest and selected target IDs. |
| `pages/after-NNNN/observations/` | Source metadata, normalized pixels, per-input checksums and cohort receipt. |
| `recurrence.json` | Cross-target rest-frame groups plus observer-frame controls, recomputed from all completed reports before the page cursor advances. |
| `ledger.sqlite` | Pinned inputs, jobs, model reservations and stage checkpoints. |
| `pages/after-NNNN/reports/` | JSON evidence packets for individual targets, including partial results on a failed or paused job. |

A failed job leaves its error in the exported packet and ledger. Inspect the error, resolve the underlying network, provider or input problem, then repeat the campaign command. Each invocation retries a failed job at most once and reuses its completed stages. A killed worker can leave a running lease for up to ten minutes. Repeat the invocation after that lease expires. Evidence packets can be reviewed independently of the database. The maximum scientific claim from the automated pipeline is a feature for follow-up. [Science and evidence standards](SCIENCE.md)
