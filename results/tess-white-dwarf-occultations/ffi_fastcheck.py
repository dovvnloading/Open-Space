"""Test each Sector 98 FFI dip on the 20-s SPOC target pixels, aligned in spacecraft time."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from ffi_images import event_masks

def validate(root, images):
    records=[]
    for r in images:
        if not r["n_isolated"] or r["sector"]!=98:continue
        p,=root.glob(f"tess*-s0098-{int(r['tic']):016d}-*_fast-tp.fits")
        with fits.open(root/f"{r['gaia']}-s98.fits") as ffi, fits.open(p) as targ:
            ffd=ffi[1].data
            ft=ffd["TIME"]
            dat=targ[1].data
            tx=dat["TIME"]
            use=np.isfinite(tx+dat["TIMECORR"])
            t=tx[use]-dat["TIMECORR"][use]
            aperture=(targ[2].data.astype(int)&2)>0
            if not aperture.any():raise ValueError("missing pipeline aperture")
            min_t=min(e["epoch_btjd"] for e in r["image_checks"])
            max_t=max(e["epoch_btjd"] for e in r["image_checks"])
            local=(t>min_t-1)&(t<max_t+1)
            index=np.flatnonzero(use)[local]
            t=t[local]
            quality=dat["QUALITY"][index]
            flux=dat["FLUX"][index][:,aperture].sum(axis=1)
            mask=(quality==0)&np.isfinite(flux)
            if abs(float(targ[1].header["TIMEDEL"])*86400-20)>.05: raise ValueError("not 20 s")
            rows=[]
            for ev in r["image_checks"]:
                center=ev["epoch_btjd"]
                i=int(np.argmin(abs(ft-center)))
                center-=float(ffd["TIMECORR"][i])
                duration=ev["duration_min"]/1440
                inside,left,right=event_masks(t,center,duration,None)
                base=mask & (abs(t-center)<.125)&(abs(t-center)>2*duration)
                ref=np.median(flux[base]) if base.any() else np.nan
                item=dict(epoch_btjd_ffi=ev["epoch_btjd"],epoch_tess_clock=center,
                          duration_min=ev["duration_min"],ffi_depth=ev["depth"],n20_raw=int(inside.sum()),
                          n20_clean=int((inside&mask).sum()),n20_shoulder_clean=int(((left|right)&mask).sum()),
                          flagged_rows=int(((~mask)&inside).sum()),
                          n20_spoc_bit4096=int(((quality & 4096 != 0)&inside).sum()),
                          spoc_quality_values_counts={str(int(val)):int(count) for val,count in
                                                      zip(*np.unique(quality[inside],return_counts=True))})
                if item["n20_clean"] and item["n20_shoulder_clean"] and base.any():
                    norm=flux/ref
                    middle=float(np.mean(norm[inside&mask]))
                    shoulder=float(np.mean(norm[(left|right)&mask]))
                    item.update(spoc20_fraction_vs_6h=1-middle,spoc20_contrast=shoulder-middle,
                                spoc20_flank_excess=shoulder-1)
                rows.append(item)
        records.append(dict(tic=r["tic"],filename=p.name,sha256=hashlib.file_digest(p.open("rb"),"sha256").hexdigest(),
                            archive_uri="mast:TESS/product/"+p.name,events=rows))
    return records

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ("input","images","out"):parser.add_argument(k,type=Path)
    args=parser.parse_args()
    args.out.write_text(json.dumps(validate(args.input,json.loads(args.images.read_text())),indent=2,allow_nan=False)+"\n")
