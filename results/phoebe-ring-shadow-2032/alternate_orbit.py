"""Regenerate the sampled PH20 minus SAT441 offset with the published PH20 SPK."""
import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
from jplephem.spk import SPK

ROOT = Path(__file__).resolve().parent
SHA256 = "a507317c8b9a03a8596068e5020a8172fb818b5ae1a4e688be7fccb366dfb89d"


def compare(path, output):
    if hashlib.sha256(path.read_bytes()).hexdigest() != SHA256:
        raise ValueError("PH20 kernel checksum mismatch")
    text = (ROOT / "checks/phoebe-saturn-barycenter.txt").read_text()
    table = csv.reader(text.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines())
    rows = []
    with SPK.open(str(path)) as spk:
        for row in table:
            t = float(row[0])
            segment = next(s for s in spk.segments if s.start_jd <= t < s.end_jd
                           and s.center == 6 and s.target == 609)
            xyz = segment.compute(t)[:3]
            rows.append([t, *(xyz - np.array(row[2:5], dtype=float))])
    np.savetxt(output, rows, fmt="%.9f", delimiter=",", header="jd_tdb,dx_km,dy_km,dz_km", comments="")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("kernel", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    compare(args.kernel, args.output)
