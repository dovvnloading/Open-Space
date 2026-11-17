"""Run the original F444W blind RDI and injection search on the GO 3989 associations."""
import argparse
import csv
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
ASSOCIATIONS = {"HD21997": "jw03989-c1007", "HD202628": "jw03989-c1009"}
# Bright, AAA-quality 2MASS hosts. Blackbody colour gives an approximate F444W
# denominator only; all absolute limits in microJy are independent of this choice.
HOSTS = {"HD21997": ("03315364-2536509", 8300.), "HD202628": ("21182725-4320048", 5750.)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("merged", type=Path, help="results/jwst-nearby-companions")
    parser.add_argument("cache", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("field", choices=ASSOCIATIONS)
    args = parser.parse_args()
    sys.path.insert(0, str(args.merged.resolve() / "nircam"))
    import nircam_reduce

    assoc = next(args.cache.glob(ASSOCIATIONS[args.field] + "*_asn.json"))
    products = json.loads(assoc.read_text())["products"]
    members = next(p["members"] for p in products if "f444w" in p["name"].lower())
    science = [m["expname"] for m in members if m["exptype"].lower() == "science"]
    refs = [m["expname"] for m in members if m["exptype"].lower() == "psf"]
    assert len(science) == 1 and len(refs) > 0
    ident, temperature = HOSTS[args.field]
    cats = list(csv.DictReader((HERE / "catalogs" / (args.field + "_2MASS.csv")).open()))
    kmag = float(next(row["k_m"] for row in cats if row["designation"] == ident))
    # Isophotal approximation. 2MASS Ks nominal zero point is 666.7 Jy.
    wave_k, wave_l = 2.159, 4.433
    ratio = (wave_k / wave_l)**3 * np.expm1(14387.77 / (wave_k * temperature)) / np.expm1(14387.77 / (wave_l * temperature))
    host_jy = float(666.7 * 10**(-0.4 * kmag) * ratio)
    with tempfile.TemporaryDirectory() as staged:
        p = Path(staged)
        for name in science + refs:
            f = args.cache / name
            with fits.open(f) as hd:
                assert hd[0].header["TARGNAME"].replace(" ", "") == args.field if name in science else hd[0].header["TARGNAME"].replace(" ", "") != args.field
            (p / name).symlink_to(f.resolve())
        for f in list(args.cache.glob("jw01537028*calints.fits")) + [args.cache / n for n in ["g191b2b_mod_012.fits", "F444W.vot", "jwst_nircam_psfmask_0019.fits"]]:
            (p / f.name).symlink_to(f.resolve())
        nircam_reduce.FIELDS[args.field] = ("jw03989", science[0], len(refs))
        sys.argv = [nircam_reduce.__file__, "--inputs", staged, "--out", str(args.out), "--field", args.field, "--host-jy", str(host_jy)]
        nircam_reduce.main()


if __name__ == "__main__":
    main()
