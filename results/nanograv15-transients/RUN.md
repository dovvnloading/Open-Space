# Run and reproduce

The NANOGrav and EPTA archives are downloaded inside the study workflow. It saves input checksums, reduces the pulsars with PINT and fits the individual and common-signal models. Four reduction partitions share two concurrent runner slots. The individual search uses the same partitioning limit. The array fit runs once. The report job verifies its evidence commit and runs the repository tests before pushing the commit to the research branch.

After registration on the default branch, use the dispatch command below. For the pre-review run the maintainer also enabled a push trigger restricted to the research branch and this workflow file. The first run downloaded both release archives. The corrected file selection can reuse those immutable Actions inputs until they expire, after which the workflow downloads the same checksummed releases. [Initial run](https://github.com/dovvnloading/Open-Space/actions/runs/44980795983) · [GitHub documentation](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow)

```sh
gh workflow run nanograv15-transients.yml --repo dovvnloading/Open-Space --ref research/nanograv15-blind-transients
```

The runner uses Python 3.11 on Ubuntu 24.04. The workflow records the installed versions and input provenance. These commands show the unpartitioned equivalent of its science steps. Paths are relative to the repository root. Run these steps inside the workflow when reproducing this study.

```sh
study=results/nanograv15-transients
python -m pip install -r "$study/requirements.txt"
python "$study/prepare.py" fetch --work raw --out receipts
for release in nb wb epta; do
  python "$study/prepare.py" reduce --dataset "$release" --work raw --out reduced
done
python "$study/search.py" --data reduced --out evidence
python "$study/array_search.py" --input reduced --output evidence/array
cp receipts/inputs.json reduced/coverage-*.json reduced/runtime-*.json evidence/
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
