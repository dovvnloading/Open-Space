# Cross-spectrum recurrence

Run the recurrence search on completed campaign evidence to rank groups of same-sign, unassigned features at similar rest wavelengths in different sky targets. The same calculation in the observer frame provides a control for stationary spectral features. Both outputs carry the status `screening_only`. The command reads the saved scans and writes its results to standard output.

```sh
python -m open_space.recurrence runs/campaign --tolerance-kms 150 --min-targets 2 > runs/recurrence.json
```

Pass a campaign root, one or more directories of reports, or explicit evidence JSON files. Multiple campaign directories can be supplied in one invocation. For a campaign root, the reader selects JSON files under its `reports/` directories. A standalone directory should contain only exported evidence packets. Incomplete jobs are counted and skipped. Completed packets must agree on software hash, detector name and schema, catalog version and threshold. Mixed configurations stop the comparison. The command operates offline.

## How to read a group

- The separation limit applies to the entire group. At the default 150 km/s, the most widely separated representative wavelengths must satisfy `c × ln(max / min) ≤ 150 km/s`. A series of nearby pairs cannot connect features farther apart than this limit. The threshold is configurable and is a screen for possible follow-up.
- Every vote comes from a distinct sky position. Repeated reviews of the same FITS digest count once. Positional neighbors within two arcseconds are conservatively merged into one target, including transitive positional matches. Set `--target-radius-arcsec` to change the aperture. Close pairs of galaxies can merge in a crowded field. Large coordinate shifts can separate repeat observations of the same physical source.
- Within a wavelength window each target contributes its highest-scoring feature. Windows rank by target count, median detector score and narrower wavelength span. Overlapping windows in each sign and frame compete for the same features. Only the preferred one is retained. Each representative records its evidence-packet feature ID and observation ID. The number of distinct plate/MJD visits describes those representatives. Absence of a count for some other visit does not measure absence of the line there.
- Compare `rest_frame_groups` with `observer_frame_control_groups`. A group concentrated at one observed wavelength calls for sky and detector checks. A small rest wavelength span across very different observed wavelengths motivates a search for a known transition. Rank scores lack survey-wide statistical calibration. Multiple tested wavelengths, missing spectral coverage, redshift errors and shared plate conditions affect apparent recurrence. Neither list alone estimates a false-discovery rate.

## Result on the first three galaxies

At the defaults, the three archived scans contained 114 features outside the bundled catalog. They produced 25 same-sign rest-frame groups involving at least two galaxies, including six groups involving all three. The observer-frame control produced eight groups involving two galaxies and none involving all three. These galaxies, from a single plate and observing date, share the same reduction. Comparing the two raw group counts does not account for chance coincidences or differing coverage in the two frames.

The second-ranked rest-frame group is absorption with median wavelength 3934.71 Å. Its members span 122.3 km/s in the rest frame and 166.43 Å in the observer frame. An investigation of the NIST reference and neighboring hydrogen lines identified Ca II K as a plausible known explanation. The original detector labels remain available for reproducibility. [Three-galaxy investigation](FIRST_INVESTIGATION.md) · [Independent measurement](INDEPENDENT_REVIEW.md)

To reproduce these counts, first fetch the three pinned observations as described in the investigation. From the checkout, enqueue those manifests, run the detector and export one evidence file per job. This example uses a separate ledger and an explicitly bounded worker.

```sh
mkdir -p runs/recurrence-reports
python -m open_space --db runs/recurrence.sqlite enqueue runs/calcium-inputs/sdss-*.json > runs/recurrence-jobs.json
python -m open_space --db runs/recurrence.sqlite work --max-jobs 3 --max-calls 0
python - <<'PY'
import json
from pathlib import Path
import subprocess
import sys

for job in json.loads(Path("runs/recurrence-jobs.json").read_text()):
    command = [sys.executable, "-m", "open_space", "--db", "runs/recurrence.sqlite", "report", job]
    packet = subprocess.check_output(command)
    (Path("runs/recurrence-reports") / f"{job}.json").write_bytes(packet)
PY
python -m open_space.recurrence runs/recurrence-reports > runs/recurrence.json
```
