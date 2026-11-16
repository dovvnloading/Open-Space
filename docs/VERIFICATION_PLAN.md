# Verification tool: plan and change record

**Status, 16 November 2026.** This plan was written after the tool was built, on a local checkout, without a branch on GitHub. It is a record of what was done and why, so that the work can be reviewed. It is not a plan frozen in advance, and the pull request should be read with that in mind.

## The question

Open-Space holds 193 result folders, most drafted by AI research agents and none peer reviewed. Each folder states numbers in prose and deposits data, checksums and scripts beside them. Nothing checked, by machine, that those parts agree with each other.

The question for this work: **which statements in the repository can be re-checked from the committed bytes alone, with no model and no opinion involved, and do they hold?**

## What the branch adds

New files only, in these places:

| Path | What it is |
| --- | --- |
| `open_space/verify/` | The tool. Standard library only, Python 3.10 or later. Run with `python -m open_space.verify`. |
| `tests/test_verify_*.py` | Its tests, including one that requires the repository to pass its own verification. |
| [VERIFICATION.md](VERIFICATION.md) | The guide: what each check does, how it controls its own error rate, and its limits. |
| `.github/workflows/verify.yml` | A workflow that runs the tool on pushes and pull requests. |
| `results/*/claims.json` (8 files) | Five declare headline numbers as formulas over committed tables, 67 claims in all. Three record a reason for a number that sits close to a different stored value. |

The five checks are: catalog (links, coverage, entry counts), checksums, printed numbers traced to deposited data, declared claims recomputed with independent arithmetic, and, on request, an offline replay of documented commands.

## What it changes in files that already existed

Each row is its own commit, so each can be reviewed, changed or dropped without the others.

| File | Change | Why | Evidence |
| --- | --- | --- | --- |
| `pyproject.toml` | Adds `open_space.verify` to the package list. | The tool would not be installed otherwise. | The package list named only `open_space`. |
| `README.md`, "How these results were made" and "Contribute" | One bullet and one command that point to the tool. | So readers and contributors can find it. | None needed; wording only. |
| `README.md`, catalog | Four entries and the entry counts. | Four merged result folders had no catalog entry: `cme-speed-storm-latency` and `forbush-storm-onset` (pull request 101), `mars-gravity-mro120g-2026` (pull request 103) and `oec-orbital-scale-audit` (pull request 100). | The catalog check lists every folder under `results/` that the README does not link. The labels and wording were chosen by the tool's author from each folder's README and need the owner's review. |
| `results/star-quasar-alignments/sources.json` | Adds a note and a `committed_table` record with the SHA-256 of the committed `gaia.csv`. | The existing entry named `gaia.csv` describes a 53,558,888-byte download that `screen.py` caches. The committed `gaia.csv` is a different, 715-byte table read by `forecast.py`, and it had no digest. | Sizes and digests of both objects. The existing fields are unchanged, so `screen.py` and `plates.py` read the file as before. |
| `results/neowise-nuclear-flares/fetch.py`, `vet.py` | Both now open a positions file whose name ends in `.gz` with `gzip`. | The README documents both commands with `data/2mrs-positions.csv.gz`, and both opened it as plain text and stopped with `UnicodeDecodeError`. | Before: the documented commands fail on the first read. After: `fetch.py` reads 44,599 positions and all 113 flare identifiers have a position. Neither script was run against the network, so the downloads themselves are untested. |

## What it deliberately does not change

- **`AGENTS.md` and `CONTRIBUTING.md`.** An earlier draft added an instruction to run the tool before a pull request. Those files set the project's process, so that decision belongs to the project owner. The edits are left out.
- **Any result's data, numbers or conclusions.** The tool reads results. Apart from the two rows above, it changes none.
- **Replays whose outputs differ.** Of 739 documented commands run offline on Windows, 27 regenerated outputs that differ from the committed ones. The cause is not established and is being looked at separately. Nothing was changed.
- **Results merged after the work started.** Live `main` moved by 197 commits while this was being built. On current `main` the tool reports nine result folders with no catalog entry (`juno-uvs-impact-control`, `maven-2025-sep-events`, `mex-vmc-2025h2-atmospheric-transients`, `mission03-dr15-thermal-ins`, `mission04-o4b-bursts`, `pace-ocean-optical-transitions`, `rubin-alerce-2026`, `rubin-jr80-mutual-motion`, `spherex-dr1-moving-spectra`) and one number that sits one unit from a stored value (`results/juno-uvs-impact-control/CLEAN_NIGHTSIDE_RESULT.md`, line 7: 12,441,537 printed, 12,441,538 stored as an exclusive row bound). These are reported, not acted on.

## How it was checked

- The tool's tests pass on Python 3.10 and 3.14 on Windows.
- The report and ledger are byte-identical under Python 3.10, 3.12 and 3.14.
- The measured result of each check is in the pull request that carries this plan, and the workflow attaches the full report to every run.

## Limits

- A number that is found in the data is consistent with the data. That does not make it correct.
- About one in eight precise numbers in result READMEs has no counterpart in the deposited data. The tool reports these as gaps and does not judge them.
- Checksums of upstream downloads are counted and not re-fetched.
- The replay was run on one Windows machine. Its outcomes depend on the installed packages.
- The workflow had never run when this plan was written. Its first run is the pull request.

## Review requested

1. Whether the tool belongs in the repository, and whether its default run should gate pull requests.
2. The four catalog entries: labels and wording.
3. The two changes inside existing result folders.
4. Whether `AGENTS.md` and `CONTRIBUTING.md` should mention the tool at all.
