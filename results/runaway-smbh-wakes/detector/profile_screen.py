#!/usr/bin/env python3
"""Annotate previously detected segments against their original science plane."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from wake_image import ImagePlane
from wake_quality import extended_profile


def main():
    p = argparse.ArgumentParser()
    p.add_argument("catalog",type=Path)
    p.add_argument("--input-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a = p.parse_args()
    d = json.loads(a.catalog.read_text())
    par = d["params"]
    x0,y0,w,h = d["window_xywh"]
    def load(product):
        m = d["inputs"][product]
        path = a.input_dir / m["name"]
        if sha256(path.read_bytes()).hexdigest() != m["sha256"]:
            raise ValueError(f"unexpected {product} hash")
        r = ImagePlane(path,m["hdu"])
        im = r.read(x0,y0,w,h)
        r.close()
        return im
    image = load("science")
    if par["invert"]:
        image = -image
    valid = np.isfinite(image)
    if par["ignore_zero"]:
        valid &= image != 0
    if "mask" in d["inputs"]:
        m = load("mask")
        valid &= np.isfinite(m)
        if par["mask_type"] == "bad":
            valid &= m == 0
        else:
            valid &= m > 0
            if par["mask_type"] == "rms":
                valid &= m < 1e8
    for row in d["detections"]:
        local = dict(row)
        for key in ("x0","x1"):
            local[key] -= x0
        for key in ("y0","y1"):
            local[key] -= y0
        row.update(extended_profile(image, valid, local))
    a.out.parent.mkdir(exist_ok=True,parents=True)
    a.out.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
