# Checking the repository against itself

Open-Space publishes many results quickly, and most were drafted by AI research agents. A reader should not have to trust the prose. One command re-checks, from the committed files alone, the statements that have exactly one right answer:

```sh
python -m open_space.verify
```

It needs Python 3.10 or later and nothing else. It reads the last commit straight from Git, takes about half a minute, prints a summary, and exits with status 1 if any machine-checkable statement in the repository is false.

The verifier does not judge whether a scientific interpretation is right. It answers narrower questions. Did the bytes change after they were hashed? Does the number in the README exist in the data? Does the stated statistic follow from the committed table? Does the documented command still run and give the committed output?

## The five checks

| Check | The question | A failure means |
| --- | --- | --- |
| Catalog | Do the links resolve, is every result folder listed, and do the stated entry counts equal the labelled entries? | The top-level README and the files behind it disagree. |
| Checksums | Does every SHA-256 recorded for a committed file equal the hash of the committed bytes? | A file changed after it was hashed, or the wrong hash was recorded. |
| Numbers | Does each number a README prints exist in the folder's deposited evidence, at the precision printed? | Nothing. An unmatched number is a gap, reported as a note. A number that almost equals a stored value is a warning. |
| Claims | Where a folder declares its key numbers in `claims.json`, does each one recompute from committed data? | The README states a number that the committed data do not give. |
| Replay | Does each `python` command a README documents run offline and reproduce the committed outputs? | Run on request. Drift and errors are warnings. |

The first four run by default. Add `--replay` for the fifth, which executes the result scripts and takes longer.

### Catalog

Every relative link in every Markdown file must lead to a file in the repository. Every folder under `results/` should be linked from the top-level README. The sentence "The catalog below has N entries: …" must state the same counts as the labelled bullets beneath it.

### Checksums

Recorded digests are found in `SHA256SUMS`-style manifests, in JSON, in tables with a digest column and in prose. Each one is tied to the committed file it speaks about. The name of the key decides this: `"sha256"` next to `"file": "x.csv"` describes `x.csv`, while `"raw_sha256"` or `"archive_sha256"` describes an upstream object. Digests of files inside committed ZIP archives and of gzip payloads count as committed bytes.

Bytes are read from the Git object store, not from the checkout, so line-ending conversion on Windows cannot produce false alarms.

A mismatch is never reported bare. The check first tries to explain it:

- the digest matches the file with different line endings, a byte-order mark or a final newline;
- the digest matches an earlier commit of the file, which was then changed without updating the record (the report names the commit);
- the record describes a different, larger upstream file that happens to share the committed file's name. This is a warning only while the committed file has no digest of its own. Once its digest is recorded and verified, the upstream entry is counted as an external receipt.

Digests of upstream downloads that are not in the repository are counted as external receipts. Checking those means fetching the files again.

### Numbers

A printed number is a claim about an interval: "12.99" asserts a value between 12.985 and 12.995. The check reads every number in a result folder's Markdown, works out that interval, and searches the folder for a machine-readable value inside it. Catalogue designations, dates, years and table or section labels are recognised and set aside.

The evidence is searched in four pools, in order:

1. values stored in JSON and in small tables;
2. facts computed from the files: row counts, counts per category and per pair of categories, column minima, maxima, sums, means and medians, list lengths, numbers inside short text;
3. the cells of large tables;
4. numeric constants in the folder's scripts.

So "843 of 944 windows" is confirmed by a table with 944 rows of which 843 carry one label, without anyone declaring it.

**A match is only evidence if a wrong number would not also have matched.** For every match the check measures how crowded the evidence is around the printed value: of the 100 neighbouring rounding cells, what fraction contain some value? That fraction is the probability that a number of the same precision, picked nearby at random, would have matched by accident. It is added up over every pool searched before the match.

- **Confirmed**: found, and a wrong number would have matched at most 5% of the time.
- **Found**: a value rounds to it, but the evidence is crowded there, so the match alone proves less.
- **Unmatched**: nothing in the folder rounds to it. This is not an error. The number has no machine-readable counterpart yet. Many are literature values quoted for context.

One minus the average chance over the matched numbers is the **sensitivity** of the check: how often a mistyped number would have been caught. The report states it for the whole repository and per folder.

A **mismatch** warning is raised only when coincidence is implausible: the printed number sits a few last-place units from a stored value, or differs from one by a single digit or a swapped pair, and the evidence is so sparse at that precision that an unrelated number would land there less than once in a thousand times.

The same search checks the top-level catalog: every precise number in a catalog bullet must appear in the prose or the evidence of the result it links to.

Only numbers with three or more significant digits count as precise. Low-precision numbers ("two stars", "about 30") match too easily to mean anything.

### Claims

A result folder may carry a `claims.json` that ties its key numbers to the data. Each claim gives the text as printed in the README and an expression that recomputes the number from the folder's own files:

```json
{
  "schema": "open-space/claims/1",
  "let": {
    "profiles": "table('profiles.csv')",
    "burst": "profiles[profiles['gbm_id'] == 211211549]"
  },
  "claims": [
    {
      "id": "hard-band-correlation",
      "printed": "0.96857",
      "statement": "Pearson correlation of STIX 50-84 keV counts with the Fermi template.",
      "compute": "pearson(burst['stix_50_84'], burst['gbm_channels_2_4_rate'])"
    },
    {
      "id": "delay-agrees-with-geometry",
      "expect": true,
      "statement": "The delay equals the geocentric delay to within the light-time across Fermi's orbit.",
      "compute": "abs(-dot(event['solar_orbiter_km'], unit(event['ra'], event['dec'])) / (c / 1000) - 30.6996) <= 7000 / (c / 1000)"
    }
  ]
}
```

A claim holds when the recomputed value rounds to exactly what is printed. The printed text must still be in the document, so editing a number in the README without updating the data makes the claim fail as stale.

The arithmetic is the verifier's own. Correlations, fits, tail probabilities, periodograms, angular separations and time scales are implemented in [`open_space/verify/stats.py`](../open_space/verify/stats.py) from the textbook definitions, with the standard library only. The result's own scripts are never imported, so a bug in them cannot hide itself.

Fields of a claim:

| Field | Meaning |
| --- | --- |
| `id` | Unique name within the file. |
| `printed` | Text exactly as it appears in the document, containing one number. Use this or `expect`. |
| `expect` | A value the expression must equal: `true`, a whole number, or a decimal (compared within a relative 1e-12 unless a tolerance is given). |
| `compute` | The expression. |
| `document` | Markdown file the text is in, relative to the claims file. Default `README.md`. |
| `let` | Named sub-expressions, for this claim or, at the top level, for all of them. |
| `tolerance` | Optional: `{"abs": x}`, `{"rel": x}` or `{"ulps": n}` (units in the last printed place). |
| `statement` | One plain sentence saying what the number is. |

The expression language is Python's expression grammar with a fixed vocabulary. There is no attribute access, no import and no file access outside the folder.

- **Data**: `table('file.csv')` gives a table; `t['column']` a column; `t[mask]` the rows where the mask holds. `json('file.json')` gives the document and `json('file.json', '/a/b/0')` one value. Gzipped files are read transparently. Blank cells are missing values and are skipped by every statistic. Whole numbers stay exact, so long catalogue identifiers compare correctly. Text columns can be tested with `startswith`, `endswith` and `contains`.
- **Columns**: arithmetic and comparisons work element by element. Combine masks with `&`, `|` and `~`. `count(mask)`, `len`, `sum`, `mean`, `median`, `min`, `max`, `std`, `pstd`, `var`, `rms`, `mad`, `quantile(v, q)`, `unique`, `nunique`, `first`, `last`, `diff`, `sorted`, `isnan`.
- **Statistics**: `pearson`, `spearman`, `slope`, `intercept`, `wmean`, `chi2(observed, model, error)`, `norm_sf`, `norm_cdf`, `norm_isf`, `chi2_sf(x, k)`, `binom_pmf`, `binom_sf(k, n, p)` (at least k), `binom_cdf`, `poisson_pmf`, `poisson_sf(k, mu)` (at least k), `poisson_cdf`, `comb`, `factorial`.
- **Periodicity**: `ls_power(t, y, period)` and `ls_peak(t, y, period_min, period_max, samples)` (Lomb-Scargle).
- **Sky, geometry and time**: `angsep(ra1, dec1, ra2, dec2)` in degrees, `unit(ra, dec)`, `dot`, `cross`, `norm`, `tt_minus_utc(date)`, `tai_minus_utc(date)`, `jd(date)`, `mjd(date)`.
- **Mathematics**: `sqrt`, `exp`, `log`, `log10`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `hypot`, `radians`, `degrees`, `floor`, `ceil`, `round`, `abs`, `erf`, `erfc`, `gamma`.
- **Constants** (SI): `pi`, `tau`, `c`, `G`, `GM_sun`, `GM_earth`, `GM_jupiter`, `au`, `pc`, `ly`, `day`, `jyear`, `arcsec`, `mas`, `deg`, `k_B`, `h_planck`, `R_sun`, `R_earth`, `R_jupiter`. The name `e` is Euler's number.

Five folders carry claims today, 67 in all: the [gamma-ray burst in the STIX flare list](../results/grb-211211a-solar-orbiter/claims.json), the [69.184-second catalog clock](../results/nicer-magnetar-microbursts/catalog-clock/claims.json), the [white-dwarf orbital decay](../results/eclipsing-dwd-gr-chirps/claims.json), the [star stream built on Reticulum II stars](../results/ngc1261b-reticulum-ii/claims.json) and the [wrong star IDs in a white-dwarf table](../results/wide-white-dwarf-triples/claims.json). They recompute, among others, a correlation from the committed light curves, TT minus UTC from the leap-second table, the Peters quadrupole formula from the tabulated masses and periods, and identifier mismatches compared as exact 19-digit integers.

### Replay

```sh
python -m open_space.verify --replay --python .venv/bin/python
```

Every `python …` command found in a result folder's Markdown is run in a scratch directory that holds an exact copy of the committed tree. Inside it the network is switched off, so the outcome does not depend on an archive being reachable. After each command the files it changed are compared with the committed ones, numbers by value and text by character, and the sandbox is restored. Gzip files are compared by their content. Images and binary containers such as `.npz` and `.fits` are noted as rewritten and not compared, because their bytes carry timestamps.

| Outcome | Meaning |
| --- | --- |
| reproduced | A documented check passed, or committed outputs regenerated identically. Numbers that agree to one part in a million count as identical. |
| ran | Exit code 0, but the command changed no committed file, so there was nothing to compare. |
| drift | Outputs regenerated with different contents. The report names the file, how many of its numbers moved and by how much. |
| failed | Non-zero exit for a reason other than the three below. |
| needs the network | The command tried to download something. The sandbox counts every refused connection, so this also covers a script that swallows the refusal and stops later on the data it never received. |
| needs a missing package | A package the script imports is not installed in the interpreter used. |
| needs inputs that are not deposited | It asked for a file or an argument that is not in the repository. This includes a file named on the command line, and a file, folder or file pattern the script searched for inside the tree and did not find before it stopped. |
| names a script that does not exist | The README documents a script that is not there. |

A stopped command is explained only from what the sandbox saw it do: the connections it attempted and the paths it looked for. Where a command sits in a document is not used, so a crash on deposited data stays a failure even when an earlier step of the same block needed the network.

`--python` chooses the interpreter that runs the scripts, so the replay can use an environment with NumPy, SciPy and Astropy while the verifier itself needs nothing.

## Reading the result

- **Failed**: a machine-checkable statement is false. The exit status is 1. Fix the file or the record.
- **Warning**: not provably wrong, but someone should look. `--strict` makes warnings fail too.
- **Note**: a gap in what can be checked mechanically, such as an unmatched number or an external receipt.

A warning that turns out to be intended can be acknowledged, with its reason, in the folder's `claims.json`. It then stays on record as a note:

```json
"acknowledged": [
  {"rule": "trace/mismatch", "key": "0.3604030",
   "reason": "Catalogue period quoted for comparison; the stored value is our own fit."}
]
```

Failures cannot be acknowledged.

## Before you open a pull request

```sh
python -m open_space.verify --worktree --result your-folder-name
```

`--worktree` checks the files as they would be committed, including new ones, without touching the Git index. To raise a folder's score:

1. Deposit the numbers the README states in a JSON or CSV file. A number that only ever appeared on a terminal cannot be traced.
2. Record a `sha256` for each deposited input next to its file name. Name upstream digests differently, for example `raw_sha256`.
3. Declare the headline numbers in `claims.json`, recomputed from the deposited tables.
4. Document a command that runs offline from the deposited files, ideally with a `--check` mode.

## Reports and continuous integration

```sh
python -m open_space.verify --report VERIFICATION.md --json verification-ledger.json
```

The Markdown report lists what was checked, everything that needs attention and one row per result folder. The JSON ledger holds every finding with its file, line and numbers, the tree digest (the SHA-256 of the tree's own `SHA256SUMS` listing, so anyone can recompute it) and the digest of the verifier's source. The same tree and the same verifier always give the same report.

The ordinary test run includes [`tests/test_verify_repository.py`](../tests/test_verify_repository.py), which fails if any statement the verifier can check is false. The [Verify workflow](../.github/workflows/verify.yml) publishes the report with every push and pull request, and replays the documented commands of the result folders a pull request touches.

## What this does not establish

- That a measurement is correct. A number can be consistent everywhere in a folder and still be the wrong number for the sky.
- That upstream archives still serve the bytes that were hashed. External receipts are counted, not re-fetched.
- That a result is new. Prior-work searches are documented in each folder and are not machine-checkable.
- That an unmatched number is wrong, or that a matched one is right for the right reason. The sensitivity figure says how much a match is worth.

Agreement between the prose, the deposited data and an independent recomputation removes transcription slips, stale text and arithmetic errors. It is a floor under the science, not a substitute for independent observation.
