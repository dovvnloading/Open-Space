# Run and reproduce

The NANOGrav and EPTA archives are downloaded only by the manual workflow. It saves the versioned input checksums, reduces individual pulsars with PINT and fits the single-pulsar and common-signal models. Four reduction partitions share two concurrent runner slots. The individual search uses the same partitioning limit. The array fit runs once. Raw downloads, residual matrices and intermediate artifacts have finite retention. The report job commits the evidence used for the catalog entry to the research branch after both repository checks pass.

GitHub requires a manually triggered workflow to exist on the default branch before its first dispatch. A branch option selects which revision runs after registration. [GitHub documentation](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow)

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
