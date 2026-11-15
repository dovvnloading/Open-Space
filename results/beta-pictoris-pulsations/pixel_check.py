"""Compare source-light and slow-signal relative column contributions on two detectors."""
import argparse
import csv
from pathlib import Path

from astropy.io import fits
import numpy as np
from scipy.ndimage import gaussian_filter1d

from measure import MODES, ROOT, coefficients, design, verified_product
from slow_signal import FREQUENCY, write_csv


def measure(tp, lc, sector):
    with fits.open(tp, memmap=True) as pixels, fits.open(lc) as target:
        d = pixels[1].data
        ref = target[1].data
        good = (d["QUALITY"] == 0) & (ref["QUALITY"] == 0)
        good &= np.isfinite(ref["PDCSAP_FLUX"]) & np.isfinite(d["TIME"])
        time = d["TIME"][good]
        cube = d["FLUX"][good]
        finite = np.all(np.isfinite(cube), axis=0)
        avg = np.nanmedian(cube, axis=0)
        keep = finite & (avg > 100)
        rows = []
        # Collapse the saturated bleed direction to a column sum.
        columns = np.sum(np.where(keep, cube, 0),axis=1).astype(float)
        segs = np.split(np.arange(len(time)),np.flatnonzero(np.diff(time)>.05)+1)
        selected = []
        for ids in segs:
            if len(ids) <= 600:
                continue
            coarse = np.arange(0,len(ids),10)
            sm = gaussian_filter1d(columns[ids[coarse]],36,axis=0,mode="nearest")
            for j in range(columns.shape[1]):
                columns[ids,j] -= np.interp(np.arange(len(ids)),coarse,sm[:,j])
            selected.extend(ids[300:-300])
        freqs = np.r_[FREQUENCY,MODES]
        t=time[selected]
        x=design(t-np.mean(t),freqs)
        b=coefficients(x,columns[selected])
        z=b[:len(freqs)]-1j*b[len(freqs):]
        for i in [0,1,2,3]:
            a=z[i]
            # Project each column on the phase of the summed signal to preserve sign.
            proj=np.real(a*np.exp(-1j*np.angle(a.sum())))
            centroid=np.dot(np.arange(len(a)),proj)/np.sum(proj)
            rows.append(dict(sector=sector,frequency=freqs[i],column_centroid=centroid,
                             signed_column_amplitudes=" ".join(f"{v:.3f}" for v in proj)))
        return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path,required=True)
    parser.add_argument("--out",type=Path,default=ROOT/"pixel_checks.csv")
    args=parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    with (ROOT/"pixel_manifest.csv").open() as f:
        stamps=list(csv.DictReader(f))
    with (ROOT/"manifest.csv").open() as f:
        lcs={r["sector"]:r for r in csv.DictReader(f)}
    rows=[]
    for row in stamps:
        tp=verified_product(row,args.cache)
        lc=verified_product(lcs[row["sector"]],args.cache)
        rows.extend(measure(tp,lc,int(row["sector"])))
    write_csv(args.out,rows)


if __name__=="__main__":
    main()
