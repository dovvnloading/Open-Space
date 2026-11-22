"""Freeze the exact numerical inputs for one shard of the period survey."""
import argparse
from pathlib import Path
import shutil
import tarfile


def freeze(work, out, shard):
    out.mkdir(exist_ok=True, parents=True)
    receipt = work / f"ingest-{shard}.jsonl"
    shutil.copy(receipt, out)
    curves = sorted((work / "curves").glob("*.npz"))
    with tarfile.open(out / f"curves-{shard:02}.tar.gz", "w:gz", compresslevel=1) as archive:
        for path in curves:
            info = tarfile.TarInfo("curves/" + path.name)
            info.size = path.stat().st_size
            with path.open("rb") as handle:
                archive.addfile(info, handle)
    print(f"Archived {len(curves)} light curves for shard {shard}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("shard", type=int)
    args = parser.parse_args()
    freeze(args.work, args.out, args.shard)
