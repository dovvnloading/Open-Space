"""Replay frozen FFI products for signed flux, detector peers, and eclipse recurrence."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

from acquire import fetch
from confirm import unpack
from scan import BIN, Curve, Binning, templates, mad, write_csv
from vet import average, read_rows, vet_rows


def recurrence(t, y, e):
    """Explore a period anchored on the event; require two *additional* recovered dips.

    Try submultiples of spacings to other local matched-amplitude maxima, then
    require the predicted timestamps to fall on actual matching dips in >=70%
    of observable cycles. This is a follow-up flag with no trials correction.
    """
    h,f,peak,dep = (float(e[k]) for k in ['duration_h','ingress_fraction','peak_btjd','depth'])
    _,_,idx,kernel = next(v for v in templates() if v[0]==h and v[1]==f)
    d = h/24
    heights = np.full(len(y),np.nan)
    heights[-idx[0]:len(y)-idx[-1]] = np.correlate(y,kernel,mode='valid')
    finite = np.isfinite(heights)
    if finite.sum() < 30:
        return {}
    cutoff=max(.003, .25*dep, 4*mad(heights[finite]))
    ix=find_peaks(np.where(finite,heights,-np.inf),height=cutoff,distance=max(2,round(.5*d/BIN)))[0]
    minimum=max(2*d,.35)
    best=None
    for dt in np.abs(t[ix]-peak):
        for j in range(1,1+int(dt/minimum)):
            p=dt/j
            if p > np.ptp(t)/2:
                continue
            cycles=np.arange(np.ceil((t[0]-peak)/p),np.floor((t[-1]-peak)/p)+1)
            offsets=peak+cycles[cycles!=0]*p
            count=0
            matched=[]
            for tt in offsets:
                ci=round((tt-t[0])/BIN)
                width=max(1,round(.15*d/BIN))
                local=heights[max(0,ci-width):min(len(t),ci+width+1)]
                if len(local) and np.isfinite(local).all():
                    count+=1
                    if np.max(local)>=cutoff:
                        matched.append(np.max(local))
            if count>=2 and len(matched)/count>=.7 and len(matched)>=2:
                score = np.sqrt(len(matched))*np.median(matched)/cutoff*len(matched)/count
                if best is None or score>best[0]:
                    best=score,p,count,len(matched),np.median(matched)
    return dict(period_trial_d=best[1],additional_covered=best[2],
                additional_matching=best[3],median_repeat_depth=best[4]) if best else {}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--products',type=Path,required=True)
    ap.add_argument('--triggers',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    st={}
    def get(r):
        raw=fetch(r['url'])
        if hashlib.sha256(raw).hexdigest()!=r['sha256']:
            raise ValueError(f"FITS hash changed for {r['tic']}")
        head,table=unpack(raw)
        meta=dict(camera=int(head['CAMERA']),ccd=int(head['CCD']),ra_obj=float(r['ra']),dec_obj=float(r['dec']))
        c=Curve(table['TIME'].astype(float),table['DET_FLUX'].astype(float),table['QUALITY']==0,{},meta)
        b=Binning(c)
        aux={k:average(b,c,table[k].astype(float)) for k in ['DET_FLUX','SAP_FLUX','SAP_BKG','SAP_X','SAP_Y',
             'DET_FLUX_SML','DET_FLUX_LAG']}
        return (r['tic'],r['sector']),(meta,b.centers,b.bin_flux(b.flux),aux)
    with ThreadPoolExecutor(8) as pool:
        for key,star in pool.map(get,[r for r in read_rows(args.products) if r['status']=='ok']):
            st[key]=star
    triggers=read_rows(args.triggers)
    out=vet_rows(st,triggers)
    for r,e in zip(out,triggers):
        _,t,y,aux=st[(e['tic'],e['sector'])]
        r.update(recurrence(t,y,e))
        d,p,f=[float(e[k]) for k in ['duration_h','peak_btjd','ingress_fraction']]
        d/=24
        sel=(t>=p-d*f)&(t<=p+d*(1-f))
        q = np.isfinite(y) & ~sel
        r['sector_corrected_flux_median']=np.median(y[q])
        r['event_corrected_flux_min_bin']=np.nanmin(y[sel])
        r['event_corrected_flux_median']=np.nanmedian(y[sel])
    write_csv(args.out,out,sorted(set().union(*(r.keys() for r in out))))


if __name__=='__main__':
    main()
