#!/usr/bin/env python3
"""Reproduce the large-stamp pilot and hash-selected deep-field injections."""
import argparse
import csv
import hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from replay_guard import frozen_tiles,validate_cutouts,fresh_output,absolute_external,seal

LARGE = {"102020065","102022484","102044820","102159484"}


def select(rows):
    # The EDF-F large-stamp center contains no usable pixels; its visibility is zero.
    # Exclude it explicitly from the pilot. The hash-based arm includes empty stamps.
    tasks = [(t,12,"nonempty pilot") for t in sorted(LARGE-{"102044820"})]
    for field in ("EDF-F","EDF-N","EDF-S"):
        available = [r["tile"] for r in rows if r["field"]==field and r["tile"] not in LARGE]
        available.sort(key=lambda t:hashlib.sha256(("wake-centres-1/"+t).encode()).hexdigest())
        tasks.extend((t,2,"hashed field sample") for t in available[:8])
    return tasks


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tiles",type=Path,required=True)
    p.add_argument("--input-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--jobs",type=int,default=3)
    a=p.parse_args()
    a.input_dir=absolute_external(a.input_dir)
    frozen_tiles(a.tiles)
    with a.tiles.open() as f:
        tasks = select(list(csv.DictReader(f)))
    validate_cutouts(a.input_dir,[t for t,_,_ in tasks])
    a.out=fresh_output(a.input_dir,a.out)
    with (a.out/"injection_sample.csv").open("w") as f:
        w = csv.writer(f,lineterminator="\n")
        w.writerow(["tile","trials_per_bin","selection"])
        w.writerows(tasks)
    def execute(task):
        t,n,_ = task
        im = a.input_dir/f"{t}_VIS.fits"
        mask = a.input_dir/f"{t}_rms.fits"
        out = a.out/f"{t}_inject.json"
        if not (im.exists() and mask.exists()):
            raise FileNotFoundError(t)
        if out.exists():
            raise FileExistsError(out)
        cmd = [sys.executable,str(Path(__file__).with_name("wake_search.py")),
               "inject",str(im),"--mask",str(mask),"--mask-type","rms","--ignore-zero",
               "--scale",".1","--zeropoint","24.6","--trials",str(n),"--seed",t,
               "--mu","23,24,25,26,27","--out",str(out)]
        subprocess.run(cmd,check=True)
    with ThreadPoolExecutor(a.jobs) as pool:
        list(pool.map(execute,tasks))
    seal(a.out,"calibration")


if __name__ == "__main__":
    main()
