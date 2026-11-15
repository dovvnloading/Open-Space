"""Check archived inputs and replay against recorded measurements."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from replay import ROOT, Image, predicted


def main():
    receipts=json.loads((ROOT/"sources.json").read_text())["inputs"]
    with (ROOT/"epochs.csv").open() as f:
        measures=list(csv.DictReader(f))
    for source, row in zip(receipts, measures, strict=True):
        input_path=ROOT/source["local"]
        assert hashlib.sha256(input_path.read_bytes()).hexdigest()==source["derived_fits_sha256"]
        image=Image(input_path)
        assert abs(image.fit(predicted(image.date))-float(row["forced_mjy"])) < 1e-6
    tracks=np.loadtxt(ROOT/"translated-tracks.csv.gz",delimiter=",",skiprows=1)
    assert tracks.shape==(8260,6)
    target=np.array([float(r["forced_mjy"]) for r in measures])
    assert np.array_equal((tracks[:,2:] >= target).sum(axis=0),[0,6931,5,6])
    assert not np.any(np.all(tracks[:,[2,4,5]]>=target[[0,2,3]],axis=1))
    print("Four FITS hashes, target measurements and 8,260 translated tracks verified.")


if __name__=="__main__":
    main()
