"""Extend the first snapshot in place with targeted follow-up products."""
import csv
import importlib
import json
import os
import traceback
from concurrent.futures import ThreadPoolExecutor
from common import OUT, coverage, events, receipts, save, status


def main():
    if os.environ.get("GITHUB_ACTIONS") != "true":
        raise SystemExit("Use the dedicated GitHub Actions workflow for archive collection.")
    with (OUT / "timeline.csv").open() as f:
        for row in csv.DictReader(f):
            for k in ("mjd","end_mjd"):
                if row[k]:
                    row[k]=float(row[k])
            events.append(row)
    with (OUT / "coverage.csv").open() as f:
        coverage.extend(csv.DictReader(f))
    receipts.extend(json.loads(line) for line in (OUT / "inputs.jsonl").read_text().splitlines())
    old = json.loads((OUT / "snapshot.json").read_text())
    (OUT / "collection-history.json").write_text(json.dumps({
        "archive_pass":{k:old.get(k) for k in ("collection_revision","workflow_run","timeline_rows")},
        "supplement_revision":os.environ["GITHUB_SHA"], "supplement_run":os.environ["GITHUB_RUN_ID"],
    }, indent=2)+"\n")
    def attempt(name):
        try:
            importlib.import_module(name).run()
        except Exception as exc:
            status(name, "unavailable", f"Supplement stopped early: {exc}")
            (OUT / (name+"-error.txt")).write_text(traceback.format_exc())
        print("Finished "+name,flush=True)
    names = ["phot_extra","xray_extra","spec_extra","lit_extra"]
    with ThreadPoolExecutor(max_workers=len(names)) as pool:
        list(pool.map(attempt,names))
    save()


if __name__ == "__main__":
    main()
