# Run and reproduce

The NANOGrav and EPTA archives are downloaded inside the study workflow. It saves input checksums, reduces the pulsars with PINT and fits the individual and common-signal models. Four reduction partitions share two concurrent runner slots. The individual search uses the same partitioning limit. The array fit runs once. The report job verifies its evidence commit and runs the repository tests before pushing the commit to the research branch.

After registration on the default branch, use the dispatch command below. For the pre-review run the maintainer also enabled a push trigger restricted to the research branch and this workflow file. The first run downloaded both release archives. The corrected file selection can reuse those immutable Actions inputs until they expire, after which the workflow downloads the same checksummed releases. [Initial run](https://github.com/dovvnloading/Open-Space/actions/runs/44980795983) · [GitHub documentation](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow)

```sh
gh workflow run nanograv15-transients.yml --repo dovvnloading/Open-Space --ref research/nanograv15-blind-transients
```

The runner uses Python 3.11 on Ubuntu 24.04. The workflow records the installed versions and input provenance. These commands show the unpartitioned equivalent of its science steps. Paths are relative to the repository root. Run these steps inside the workflow when reproducing this study.

```sh
study=results/nanograv15-transients
python -m pip install -r "$study/requirements.txt" pandas matplotlib
python "$study/prepare.py" fetch --work raw --out receipts
for release in nb wb epta; do
  python "$study/prepare.py" reduce --dataset "$release" --work raw --out reduced
done
python "$study/epta_modern.py" --work raw --out reduced --shard 0
python "$study/search.py" --data reduced --out evidence
python "$study/array_search.py" --input reduced --output evidence/array
python "$study/coherence_check.py" --data reduced --evidence evidence --out evidence/coherence
python "$study/dm_control.py" --data reduced --evidence evidence --out evidence/dm-control
python "$study/known_controls.py" --data reduced --evidence evidence
python "$study/figure.py" --evidence evidence
cp receipts/inputs.json reduced/coverage-*.json reduced/runtime-*.json evidence/
python -m pip freeze > evidence/software-0.txt
python "$study/report.py" --evidence evidence --root .
```

The workflow verifies its evidence commit with the repository's two requested commands.

```sh
python -m open_space.verify
python -m unittest discover -s tests
```

## Timing-file compatibility

The first pass found 68 pulsars in each NANOGrav reduction and could select a unique file pair for 62. Six sources have both combined and site-specific solutions in the archive. The revised selection uses the combined solution for J1600−3053, J1643−1224, J1713+0747, J1903+0327, J1909−3744 and B1937+21. The receipt identifies the specific files and SHA-256 values. Wideband processing applies the release flags before building a set of paired arrival-time and DM measurements. Its receipt counts removed and unpaired rows separately.

The independent reduction interprets EPTA's TEMPO2 comment convention and records observatory clock coverage. It holds the pulse counts fixed when refitting the released parameters after PINT's TCB-to-TDB conversion. The receipt retains the initial and refitted WRMS and the original published TRES. The same acceptance window applies to that original TRES. A rejected reduction remains visible in the coverage table.

Each EPTA timing file has its own TEMPO2 time-offset and uncertainty directives. The adapter records their line numbers and effective values. When the full timing pair fails a quality gate, the workflow also tries the collaboration's published DR2new timing and parameter files from the same archive. DR2new covers a shorter modern-backend observing span. The same phase and residual checks apply, and only one reduction per pulsar enters the external comparison. File hashes, release selection, observing span and failed attempts are in the EPTA coverage receipts. [EPTA DR2 selection](https://doi.org/10.1051/0004-6361/202346844)

For a post-discovery check, the workflow can reuse retained artifacts from an earlier run of this workflow. It checks the run, completed jobs, artifact availability, commands and code revision. A change confined to the EPTA reader reuses the NANOGrav residuals and repeats the EPTA reduction and validation. A complete primary-search replay also checks the individual and coherent search code. The run receipt identifies the source runs. A dispatch with both optional run IDs blank reproduces every stage from the archives.

The final run repeated both NANOGrav reductions from the archived timing inputs. It recorded 5,394 primary crossings. The preceding run recorded 5,395. The difference was a 14-day Doppler template for J1640+2224 at MJD 57122.5, whose earlier survey-adjusted p was 0.00999853 at a 0.01 threshold. It already failed the observing-day deletion and EPTA checks. The full rerun retained the same 77 internal common-pattern fits and 368 calibrated individual limit cells. The published tables and file hashes refer to the final run. [Final run](https://github.com/dovvnloading/Open-Space/actions/runs/44986027784) · [Preceding run](https://github.com/dovvnloading/Open-Space/actions/runs/44984927588)
