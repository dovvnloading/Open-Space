"""Freeze disjoint O4b calendar slots using only GWOSC DATA and the archive index.

python holdout_data.py DATA_CACHE results/holdout_selection.json
acquire_holdout(cache, block) downloads and verifies one H1/L1 pair.
"""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from acquire import BASE, DURATION, FS, META, RELEASE, RUN_END, fetch

RUN_START = 1396796418
GPS_MINUS_UNIX = 18 - 315964800  # constant for O4b
INPUT_HASHES = {
    "H1_DATA.txt": "b25f3d612babb02b6cca5f20238c7ea6d90a313e82c7aafdadd8585ab5f26733",
    "L1_DATA.txt": "15b83906df1907f3af64fef79c9176bb240b826d029ada2faa4f41e620d9176c",
    "MD5.txt": "b3d6d9fa4f176df7b796995a28e70d7cb910a68d1785c602405180551a7eb232",
}
PILOT = (1397198848, 1399799808, 1402445824, 1404825600, 1408577536,
         1410396160, 1414307840, 1415680000, 1418264576, 1420955648)
CONTROL = 1420877824


def utc(gps):
    return datetime.fromtimestamp(gps-GPS_MINUS_UNIX, timezone.utc)


def freeze(cache):
    cache.mkdir(parents=True, exist_ok=True)
    for name, digest in INPUT_HASHES.items():
        fetch(META[name], cache/name)
        if hashlib.sha256((cache/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Pinned input hash mismatch: {name}")
    timelines = {d: np.loadtxt(cache/f"{d}_DATA.txt", dtype=int)[:,:2] for d in ("H1","L1")}
    index = {}
    for line in (cache/"MD5.txt").read_text().splitlines():
        md5, path = line.split()
        det, directory, name = path.split("/")
        gps = int(name.rsplit("-",2)[1])
        index[(gps,det)] = {"name":name, "url":f"{BASE}/archive/data/{RELEASE}/{directory}/{name}",
                            "md5":md5}
    excluded = set(PILOT) | {CONTROL}
    date = utc(RUN_START).replace(hour=0,minute=0,second=0)
    monday = date-timedelta(days=date.weekday())
    blocks, empty = [], []
    while monday < utc(RUN_END):
        for day,stop_day in ((0,2),(2,4),(4,6),(6,7)):
            start_date = monday+timedelta(days=day)
            lo = max(int(start_date.timestamp())+GPS_MINUS_UNIX,RUN_START)
            hi = min(int((monday+timedelta(days=stop_day)).timestamp())+GPS_MINUS_UNIX,RUN_END)
            if lo >= hi:
                continue
            start = (lo+DURATION-1)//DURATION*DURATION
            slot = {"target_utc":start_date.isoformat(),"slot_gps":[lo,hi]}
            for gps in range(start,hi-DURATION+1,DURATION):
                if gps in excluded:
                    continue
                if not all(np.any((s[:,0]<=gps)&(gps+DURATION<=s[:,1])) for s in timelines.values()):
                    continue
                sources = {d:index[(gps,d)] for d in ("H1","L1")}
                # The Virgo receipt is available to post-ranking diagnostics
                # without affecting H/L selection or consuming bulk bandwidth.
                if (gps,"V1") in index:
                    slot["virgo_source"] = index[gps,"V1"]
                blocks.append({**slot,"gps":gps,"role":"holdout","month":utc(gps).strftime("%Y-%m"),
                               "files":{d:r["name"] for d,r in sources.items()},"sources":sources})
                break
            else:
                empty.append(slot)
        monday += timedelta(days=7)
    starts = [b["gps"] for b in blocks]
    assert len(starts) == len(set(starts)) and not excluded.intersection(starts)
    return {"release":RELEASE,"duration":DURATION,"fs":FS,"bounds_gps":[RUN_START,RUN_END],
            "pilot_exclusions_gps":list(PILOT),"control_exclusion_gps":CONTROL,
            "inputs":{k:{"url":META[k],"sha256":d} for k,d in INPUT_HASHES.items()},
            "blocks":blocks,"unfilled_slots":empty}


def acquire_holdout(cache, block):
    """Return detector->local Path, after both complete 4 kHz HDFs match official MD5.

    The caller owns the per-block files and may unlink them after durable
    output is written. The only filenames ever published are verified files.
    """
    cache = Path(cache)
    cache.mkdir(parents=True,exist_ok=True)
    def one(det):
        r = block["sources"][det]
        path = cache/r["name"]
        if path.exists():
            # Verify files even if a previous interrupted run left them.
            with path.open("rb") as stream:
                if hashlib.file_digest(stream,"md5").hexdigest() == r["md5"]:
                    return det,path
        # Unique same-filesystem temp directory allows concurrent readers
        # and retries without exposing partially received files.
        with TemporaryDirectory(dir=cache,prefix=".strain-") as tmp:
            local = Path(tmp)/path.name
            fetch(r["url"],local,r["md5"])
            local.replace(path)
        return det,path
    with ThreadPoolExecutor(max_workers=2) as pool:
        return dict(pool.map(one,("H1","L1")))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    data=freeze(args.cache)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2)+"\n")
    print(f"{len(data['blocks'])} selected, {len(data['unfilled_slots'])} unfilled")
