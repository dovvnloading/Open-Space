# Open verifier warnings: plan and change record

**Status, 16 November 2026.** Written before any file was changed, on branch `validation/open-warnings`, cut from `main` at `f02ee9b`. The "Outcome" column is filled in as each item is worked, one commit per item or group.

## The question

A full run of `python -m open_space.verify --rev HEAD --replay` on `main` at `f02ee9b` reports 0 failed checks and 25 warnings. **For each warning: what causes it, and is it a defect in the result, a defect in the verifier, or an expected difference that should be recorded?**

The run was made on one machine: Windows 11, Python 3.14, NumPy 2.5.4, SciPy 1.18.1.

## Rules for this branch

- The cause is established first. Nothing is changed for a warning whose cause is not known.
- No result's printed number or conclusion is changed. If a rerun shows a printed number is wrong, that is reported in the pull request for the owner and left as printed.
- A script is changed only to remove an accident of the machine (file order, round-off written where a value is zero by construction). Every such change is checked by regenerating the output.
- A difference that comes from the numerical libraries is documented as a tolerance in the result's own README and recorded as an acknowledgement in its `claims.json`. Its outputs are left untouched.
- If the warning is the verifier's mistake, the verifier is fixed and a test is added. The result is left alone.
- Catalog wording is taken from each result's own README. Labels and wording need the owner's review.
- `AGENTS.md` and `CONTRIBUTING.md` are not changed.

## The 25 warnings

### A. Nine result folders with no catalog entry

`juno-uvs-impact-control`, `maven-2025-sep-events`, `mex-vmc-2025h2-atmospheric-transients`, `mission03-dr15-thermal-ins`, `mission04-o4b-bursts`, `pace-ocean-optical-transitions`, `rubin-alerce-2026`, `rubin-jr80-mutual-motion`, `spherex-dr1-moving-spectra`.

Plan: read each folder, establish whether it is a finished result or work in progress, and add a catalog entry only for finished ones. Raise the entry counts to match.

### B. One number one unit from a stored value

`results/juno-uvs-impact-control/CLEAN_NIGHTSIDE_RESULT.md`, line 7: prints 12,441,537; the evidence stores 12,441,538 as `last_row_exclusive`.

Plan: confirm from the evidence that the printed number is the last row and the stored one is the exclusive bound, then record that as an acknowledgement.

### C. Ten replays whose outputs differ

| Result | Command | Size of the difference |
| --- | --- | --- |
| `mission05-swot-rivers` | `python congo/profile/compare_d.py` | 88 of 183 numbers, up to 180% |
| `mission04-o4b-bursts` | `python noise_controls.py noise_controls.json` | 4 of 3,357 numbers, up to 100% |
| `mission03-dr15-thermal-ins` | `python analysis/ledger.py results overlap cache` | different text |
| `tic-43472931` | `python timing.py` | 2 of 48 numbers, up to 2.2e-04 |
| `tic-96100077` | `python ground-check/eclipse_controls.py` | 2 of 34 numbers, up to 2.3e-05 |
| `koi-7865-eclipse-source` | `python ztf_clock.py` | 7 of 49 numbers, up to 1.3e-05 |
| `smc-04602` | `python delve04602.py data data/iv-04602-I.dat` | 1 of 64 numbers, 8.6e-06 |
| `v612-sct-clock` | `python replay.py`, `python figure.py` | 1 of 101 numbers, 2.9e-06 |
| `rubin-jr80-mutual-motion` | `python sensitivity.py` | 8 of 271 numbers, up to 1.9e-06 |

Plan: for each, compare the regenerated file with the committed one, name the values that move, and trace them to the line of the script that produces them. The three large ones come first, because a difference that size may reach a printed number.

### D. Two commands that fail

`results/pace-ocean-optical-transitions/environment`: `python regional_hydrography.py` and `python black_sea_reference.py` stop with `TypeError: argument of type 'NoneType' is not a container or iterable`.

Plan: read the traceback. Decide whether the script needs an input that is not deposited (then the verifier should say so), or has a defect.

### E. Three documented scripts that do not exist

`results/spherex-dr1-moving-spectra`: `fit_psg_alternatives.py`, `fit_psg_feature.py`, `swan_broadening_review.py`.

Plan: search the repository and its history for the scripts. If they exist under another path, correct the documented command. If they were never deposited, report that and change nothing.

## Outcome

Filled in on 16 November 2026, after the work. Each row is one commit unless it says otherwise.

| Item | Cause | Change | Left open |
| --- | --- | --- | --- |
| A. Eight unlisted results | Merged without a catalog entry. | Eight entries in `README.md`; counts raised from 176 to 184 (86 measured, 44 candidate, 20 forecasts, 20 corrections, 14 null). | Labels and wording need the owner's review. |
| A. `mex-vmc-2025h2-atmospheric-transients` | The folder holds a search plan and no result. | Recorded as an acknowledgement in a new `claims.json`. No catalog entry. | It needs an entry once the search is run. |
| B. `juno-uvs-impact-control`, 12,441,537 | The sentence prints the last row, inclusive. The evidence stores the exclusive bound, 12,441,538. First row 12,431,751 and 9,787 rows agree with both. | Acknowledgement in a new `claims.json`. | Nothing. |
| C. `mission05-swot-rivers`, `compare_d.py` | The regenerated table holds the same 16 rows in a different order. The script sorts by revision with an unstable sort, so tied rows come out in a library-dependent order. | The script sorts the rows by pair, reach, bin and node before writing; `d-node-common.csv` regenerated in that order. No value changes. | Nothing. |
| C. `mission04-o4b-bursts`, `noise_controls.py` | The output records the NumPy and SciPy versions it ran under. Those two strings are the only values that differ. | Acknowledgement in a new `claims.json`. | Nothing. |
| C. `mission03-dr15-thermal-ins`, `analysis/ledger.py` | Verifier mistake. The command reads a download cache that is not deposited, finished anyway and wrote empty tables. | Verifier: a command that finishes with different outputs remains a drift warning; missing-input or refused-network probes are reported as context because their effect is uncertain. Test added. The ledger now stops before writing when its required counterpart cache is absent or empty. | Cached inputs are not deposited, so a full replay remains pending. |
| C. Seven commands in `koi-7865-eclipse-source`, `rubin-jr80-mutual-motion`, `smc-04602`, `tic-43472931`, `tic-96100077`, `v612-sct-clock` | An iterative fit or linear solve whose last digits depend on the installed NumPy and SciPy. The largest change is 1.3 parts in 100,000, apart from one interval bound near zero. The printed-number check finds the same numbers on the regenerated outputs as on the committed ones in all six folders. | One commit: a short note beside each replay command, an acknowledgement in each `claims.json`, and updated `SHA256SUMS` entries for the three READMEs they list. Committed outputs untouched. | The differences remain; they are documented, not removed. |
| D. `pace-ocean-optical-transitions`, two commands | Verifier mistake. Both load the system Blosc library through `ctypes`. It is not installed on the test machine, and the README says it is required. | Verifier: an actual loader error naming an unavailable shared library is reported as a missing dependency; optional or unreachable lookups do not prevent replay. Test added. The regional reader now attempts a named library rather than loading the process when lookup fails. | Not replayed with Blosc installed; no measured values were changed. |
| E. `spherex-dr1-moving-spectra`, three scripts | The scripts are deposited under `examples/`, and the reports named them without a path. Two match the digests recorded in the folder's hash files. The third, `spherex_swan_broadening.py`, takes the same arguments and writes the same keys as `swan_broadening_review.json`, and was added in the same commit. | The three commands name `../../examples/<script>`. The recorded digests of the two edited reports are updated. | All three still need inputs that are not deposited, so none was replayed to the end. |

## How it will be checked

- The full verifier run, with replay, on the branch tip.
- The verifier's own tests.
- The workflow on the pull request, which replays the touched results on Linux.

## Limits

- The replay runs on one Windows machine and, for touched results, on Linux in the workflow. A command that reproduces on both may still differ elsewhere.
- A recorded tolerance says the difference is understood. It does not say the printed number is right.
