#!/usr/bin/env python3
"""Run the Euclid screen on downloaded cutouts, with checkpoints per tile and sign."""
import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor


def run(tile, args):
    image = args.input_dir / f"{tile}_VIS.fits"
    mask = args.input_dir / f"{tile}_rms.fits"
    if not (image.is_file() and mask.is_file()):
        return
    for name in ("search","invert"):
        dest = args.out / f"{tile}_{name}.json"
        if dest.exists() and not args.replace:
            continue
        cmd = [sys.executable, str(Path(__file__).with_name("wake_search.py")), "search",
               str(image), "--mask", str(mask), "--mask-type", "rms", "--ignore-zero",
               "--scale", ".1", "--lengths", "6,10", "--widths", ".2",
               "--threshold", "6", "--min-length", "3", "--out", str(dest)]
        if name == "invert":
            cmd.append("--invert")
        subprocess.run(cmd,check=True,env={**os.environ, "OPENBLAS_NUM_THREADS":"1"})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--jobs",type=int,default=4)
    p.add_argument("--replace",action="store_true")
    p.add_argument("--tiles",type=Path,help="optional frozen CSV with a tile column")
    args = p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    if args.tiles:
        with args.tiles.open() as f:
            tiles = [r["tile"] for r in csv.DictReader(f)]
    else:
        tiles = sorted(f.name.split("_")[0] for f in args.input_dir.glob("*_VIS.fits"))
    with ThreadPoolExecutor(args.jobs) as pool:
        list(pool.map(lambda t:run(t,args),tiles))
