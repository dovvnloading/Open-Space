"""Recompute the two difference-photometry tables from the downloaded cutouts."""
import argparse
import gzip
from pathlib import Path
import sys
import numpy as np


HERE=Path(__file__).resolve().parent
FIELDS="jd x_sub y_sub path_x path_y segment quality diff_mad psf_flux aperture_flux first_ref_zeros n_refs median_sky".split()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cutouts",type=Path,help="Parent of stamps29489_20 and stamps3177_39")
    p.add_argument("--out",type=Path,default=HERE)
    p.add_argument("--locked-code",type=Path,default=HERE/"locked",help="Directory containing the frozen image.py")
    args=p.parse_args();args.out.mkdir(exist_ok=True,parents=True)
    sys.path.insert(0,str(args.locked_code.resolve()))
    from image import get,reduce
    for number,sector in [(29489,20),(3177,39)]:
        path=args.cutouts/f"stamps{number}_{sector}"
        frames=get(path)
        if not frames:
            raise ValueError(f"No input cutouts in {path}")
        rows=reduce(frames)[0]
        with gzip.open(args.out/f"photometry_{number}.csv.gz","wt") as f:
            np.savetxt(f,rows,delimiter=",",header=",".join(FIELDS),comments="",fmt="%.17g")


if __name__=="__main__":
    main()
