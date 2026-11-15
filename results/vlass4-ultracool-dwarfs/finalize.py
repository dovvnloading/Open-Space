"""Measure all retrieved cutouts and make the public, deterministic receipt ledger."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import pandas as pd
from astropy.io import fits
from survey import measure
from summarize import summarize


def record(row, cache):
    path = cache / "fits" / f"{row.id}.fits"
    if not path.exists():
        return dict(id=row.id, status="archive unavailable")
    try:
        m = measure(path.read_bytes(), row)
        with fits.open(path, memmap=False) as h:
            method = "byte range" if h[0].header.get("OSRANGE") else "SODA"
        return dict(id=row.id, method=method, status="measured", **m)
    except (ValueError, OSError) as err:
        return dict(id=row.id, status=str(err))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path("."))
    a = p.parse_args()
    targets = pd.read_csv(a.cache / "targets.csv")
    with ThreadPoolExecutor(max_workers=3) as e:
        rows = list(e.map(lambda r: record(r, a.cache), targets.itertuples()))
    result = targets.merge(pd.DataFrame(rows), on="id")
    result.to_csv(a.cache / "screen-final.csv", index=False)
    a.out.mkdir(exist_ok=True, parents=True)
    result.to_csv(
        a.out / "cutout-ledger.csv.gz",
        index=False,
        float_format="%.9g",
        compression=dict(method="gzip", compresslevel=9, mtime=0),
    )
    report = summarize(result, pd.read_csv(a.cache / "ucd.csv", low_memory=False))
    report["retrieval_status"] = result.status.value_counts().to_dict()
    report["retrieval_method"] = result.method.value_counts().to_dict()
    (a.out / "coverage.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
