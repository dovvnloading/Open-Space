"""Search each whole local shard, preserving every individual S/N trigger."""
import argparse
from collections import Counter
import csv
from dataclasses import asdict
import json
from pathlib import Path

from desi_shard import targets
from laser_search import scan_target


def scan(cache, survey, program, healpix, output):
    output = Path(output)
    output.mkdir(exist_ok=True, parents=True)
    tag = f"{survey}-{program}-{healpix}"
    counts = Counter()
    dest = output / f"scan-{tag}.jsonl"
    tmp = dest.with_suffix(".partial")
    with tmp.open("w") as file:
        for meta, arms in targets(cache, survey, program, healpix, counts):
            ids = [r["expid"] for r in meta["exposures"]]
            arm_results = {}
            for arm, records in arms.items():
                res = scan_target(records, ids)
                arm_results[arm] = asdict(res)
                counts['searched_arms'] += 1
                counts['triggers'] += sum(len(h) for h in res.exposure_hits.values())
                counts['persistent'] += len(res.candidates)
            print(json.dumps(dict(meta=meta, arms=arm_results), allow_nan=False, separators=(',', ':')),
                  file=file, flush=True)
    tmp.replace(dest)
    (output / f"count-{tag}.json").write_text(json.dumps(counts, indent=2) + "\n")
    print(tag, dict(counts), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache")
    p.add_argument("survey")
    p.add_argument("program")
    p.add_argument("healpix", type=int)
    p.add_argument("output")
    a = p.parse_args()
    scan(a.cache, a.survey, a.program, a.healpix, a.output)
