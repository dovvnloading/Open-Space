"""Extract DE440s barycentric states used to initialize the planetary integrations."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from jplephem.spk import SPK

SHA256 = "c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2"
AU = 149597870.7
GMS = [132712440041.27942,22031.86855,324858.592,403503.235625,42828.375816,
       126712764.1,37940584.8418,5794556.4,6836527.10058]
EPOCHS = (2459317.5,2458780.5,2461200.5,2461000.5)


def extract(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != SHA256:
        raise ValueError("DE440s SHA-256 mismatch")
    eps = np.deg2rad(23.439291111111)
    co,si = np.cos(eps),np.sin(eps)
    rotate = np.array([[1,0,0],[0,co,si],[0,-si,co]])
    states = {}
    with SPK.open(str(path)) as kernel:
        for epoch in EPOCHS:
            vectors = []
            for body, gm in zip([10,1,2,3,4,5,6,7,8],GMS):
                pos,vel = kernel[0,body].compute_and_differentiate(epoch)
                vectors.append([gm,*(rotate@pos/AU),*(rotate@vel/AU)])
            states[str(epoch)] = vectors
    return states


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("kernel",type=Path)
    p.add_argument("destination",type=Path)
    args = p.parse_args()
    args.destination.write_text(json.dumps(extract(args.kernel),indent=2)+"\n")
