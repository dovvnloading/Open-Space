"""Quantify changes in blank-sky time-series noise independent of the target aperture."""
import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter

def blank_sky(root,curves,rec):
    with fits.open(root/f"{rec['source_id']}-s{rec['sector']}.fits") as ff:
        with np.load(curves/rec["lc_file"]) as lc:
            t,rows=lc["time"],lc["pixel_row"]
        cube=np.asarray(ff[1].data["FLUX"][rows],float)
    yy,xx=np.indices(cube.shape[1:])
    edge=(xx==0)|(yy==0)|(xx==xx.max())|(yy==yy.max())
    pix=cube[:,edge]
    diff=np.diff(pix,axis=0,prepend=pix[:1])
    center=np.nanmedian(diff,axis=1)
    frame_noise=1.4826*np.nanmedian(abs(diff-center[:,None]),axis=1)/np.sqrt(2)
    # 43 exposures ~ 2.4 h, chosen to compare local noise with the full sector.
    series=median_filter(frame_noise,size=43,mode="nearest")
    return t,series

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs","curves","images","out"):parser.add_argument(arg,type=Path)
    args=parser.parse_args()
    images=json.loads(args.images.read_text())
    specs=json.loads((args.curves/"ffi_measurements.json").read_text())
    records=[]
    for row in specs:
        if row["status"]!="included":continue
        t,noise=blank_sky(args.inputs,args.curves,row)
        im=next(r for r in images if (r["tic"],r["sector"])==(row["tic"],row["sector"]))
        out=dict(tic=row["tic"],sector=row["sector"],full_visit_blank_sky_noise_e_s=float(np.median(noise)),events=[])
        for ev in im["image_checks"]:
            index=int(np.argmin(abs(t-ev["epoch_btjd"])))
            local=noise[index]
            out["events"].append(dict(epoch_btjd=ev["epoch_btjd"],diagnostic=ev["diagnostic"],
                            blank_sky_noise_e_s=float(local),factor_over_visit_median=float(local/np.median(noise)),
                            percentile_within_visit=100*np.mean(noise<=local)))
        records.append(out)
    args.out.write_text(json.dumps(records,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":main()
