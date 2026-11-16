# Reproduce the O4b searches

Use Python 3.12 on a machine with outbound HTTPS and curl. The original run used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, h5py 3.16.0 and Matplotlib 3.10.8. LALSuite 8.0.2 supplies GWF support for auxiliary follow-up. Commands below run from this research directory. STORE is a writable data cache outside the checkout. The holdout streams 250 large HDF5 files and can release each pair after a verified checkpoint. Four workers can be used on a machine with approximately 6 GiB free RAM. Change worker count for available memory.

```bash
python -m pip install numpy==2.3.5 scipy==1.17.0 h5py==3.16.0 matplotlib==3.10.8 lalsuite==8.0.2
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export STORE="$PWD/cache"
mkdir -p "$STORE"
python -m unittest -v test_core
python noise_controls.py noise_controls.json
```

## First ten epochs and the point control

The committed manifest protects the exact data selection from later catalog changes. To replay those files, copy the frozen manifest and fetch each object with the official checksum. A second route, `python acquire.py "$STORE"`, reconstructs the manifest from live archive metadata and downloads it. The run-time catalogue vetoes always use the pinned `public_events.csv`.

```bash
cp results/search_selection.json "$STORE/search_selection.json"
python - <<'PYCODE'
import json, os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from acquire import fetch
cache = Path(os.environ["STORE"])
m = json.loads((cache/"search_selection.json").read_text())
def get(item):
    name, record = item
    fetch(record["url"], cache/name, record["md5"])
with ThreadPoolExecutor(max_workers=6) as workers:
    list(workers.map(get, m["manifest"].items()))
PYCODE
python survey.py "$STORE" --tail-counts
python sensitivity.py "$STORE"
python detector_context.py "$STORE" results
python plot_background.py
```

Run the per-month auxiliary and line checks with the original 64-second public frame for each epoch's loudest foreground tile.

```bash
python - <<'PYCODE'
import json, os, subprocess, sys
from pathlib import Path
dest = Path("results/diagnostics")
dest.mkdir(exist_ok=True, parents=True)
for p in sorted(Path("results/epochs").glob("*.json")):
    e = json.loads(p.read_text())
    if e["role"] != "search":
        continue
    t = e["top"][0]
    subprocess.run([sys.executable, "diagnose.py", os.environ["STORE"],
                    str(dest/p.name), *map(str,(t["gps"],t["duration"],*t["band"]))], check=True)
PYCODE
```

## Disjoint 125-frame extension

```bash
python holdout.py "$STORE/stream" "$STORE/checkpoints" --workers 4
python common_exposure.py
python holdout_followup.py calibration "$STORE"
python holdout_followup.py follow "$STORE"
python audit.py --results results --pilot-cache "$STORE" --output results/audit_results.json
python plot_campaigns.py
```

To resume an interrupted extension, repeat the same command and checkpoint path. The runner validates the configuration fingerprint and CRC before reusing per-frame results. To rebuild the monthly ZIPs and summary without downloading any strain, run `python holdout.py "$STORE/stream" "$STORE/checkpoints" --aggregate-only`. To verify the calendar selection independently, fetch the H1 and L1 timeline text and MD5 index from the recorded URLs under the frozen manifest's `inputs`. Then run `python holdout_data.py "$STORE" "$STORE/regenerated-selection.json"` and compare its SHA-256 to the committed file.

The 10 monthly ZIPs contain a folder per GPS start. Each holds the detector masks and 129 shift summaries in `epoch.json`, the exact unshifted score for every searchable second in `foreground.csv`, all shifted and unshifted scores ≥3 in `tail.csv`, and every assigned sensitivity trial in `injections.json`. Each offset row identifies both detectors' GPS epochs. The campaign JSON lists the monthly SHA-256 hashes. The first-stage ledgers additionally retain the winning tile metadata for each positive foreground score. Source strain remains available from GWOSC via the checksummed manifests.

## Independent checks

The five core regressions cover localization, polarity, physical delays, excluded seconds, nonwrapping slide exposure, rejected inputs, amplitude normalization and crop equivalence. For the repository suite, install `.[archive]` at the repository root and run `python -m unittest discover -s tests`. To refresh the three public-trigger sources for comparison in a new directory, run `python vetting.py "$STORE/vetting" "$STORE/refreshed-vetoes"`. Archives and GraceDB may revise public metadata, so keep the published snapshot when reproducing the result.
