"""Test whether the simultaneous source variations predict KOI-5275's eclipse depth."""
import argparse
import csv
import hashlib
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits

HERE=Path(__file__).resolve().parent
SOURCE,TARGET=6367628,6367663
QUARTERS={14:"2012277125453",15:"2013011073258"}
PERIOD=3.77990376255
# Time convention BKJD. A deeper eclipse of the actual source will be selected
# in Q14, then held fixed through Q15. Transit phase of KOI-5275 is not fitted.


def get_lightcurve(cache,kic,stamp):
    name=f"kplr{kic:09}-{stamp}_llc.fits"
    file=cache/name
    if not file.exists():
        file.write_bytes(urlopen(f"https://stpubdata.s3.amazonaws.com/kepler/public/lightcurves/{str(kic).zfill(9)[:4]}/{kic:09}/{name}",timeout=90).read())
    with (HERE/"input_manifest.csv").open() as handle:
        manifest = {row["filename"]: row["sha256"] for row in csv.DictReader(handle)}
    actual_hash = hashlib.sha256(file.read_bytes()).hexdigest()
    if actual_hash != manifest[name]:
        raise ValueError(f"Kepler product checksum mismatch: {name}")
    with fits.open(file) as h:
        data=h[1].data
        good=(data.SAP_QUALITY==0)&np.isfinite(data.TIME)&np.isfinite(data.PDCSAP_FLUX)
        flux=np.asarray(data.PDCSAP_FLUX[good],dtype=float)
        return data.CADENCENO[good],data.TIME[good],flux/np.median(flux)-1


def main(cache):
    output=[]
    center=None
    for quarter,stamp in QUARTERS.items():
        sc,_,source=get_lightcurve(cache,SOURCE,stamp)
        tc,t,target=get_lightcurve(cache,TARGET,stamp)
        _,i,j=np.intersect1d(sc,tc,return_indices=True)
        t,src,tgt=t[j],source[i],target[j]
        phase=(t/PERIOD)%1
        if center is None:
            edges=np.linspace(0,1,101)
            med=[np.median(src[(phase>=a)&(phase<b)]) for a,b in zip(edges[:-1],edges[1:])]
            center=(edges[np.argmin(med)]+.005)%1
        phi=(phase-center+.5)%1-.5
        primary=abs(phi)<.047
        base=(abs(phi)>.15)&(abs(phi)<.35)
        orbit=np.floor(t/PERIOD-center+.5).astype(int)
        cycles=np.unique(orbit)
        # Train a linear dilution on the opposite quarter below. Here measure event-local
        # contrasts for both stars and fit their orbit-demeaned light curves.
        pairs=[]
        for o in cycles:
            g=orbit==o
            if (g&primary).sum()<2 or (g&base).sum()<4:
                continue
            a,b=src[g]-src[g].mean(),tgt[g]-tgt[g].mean()
            pairs.append([a@b,a@a,src[g&base].mean()-src[g&primary].mean(),
                          tgt[g&base].mean()-tgt[g&primary].mean()])
        arr=np.array(pairs)
        output.append((quarter,arr))
    rng=np.random.default_rng(5275)
    rows=[]
    for (q,obs),(_,train) in [output,output[::-1]]:
        slope=train[:,0].sum()/train[:,1].sum()
        expected=slope*np.mean(obs[:,2]);measured=np.mean(obs[:,3])
        res=[]
        for _ in range(3000):
            a=train[rng.integers(len(train),size=len(train))]
            b=obs[rng.integers(len(obs),size=len(obs))]
            pred=a[:,0].sum()/a[:,1].sum()*b[:,2].mean()
            res.append((b[:,3].mean()-pred)*1e6)
        lo,hi=np.percentile(res,[2.5,97.5])
        rows.append(dict(quarter=q,source_phase=center,orbits=len(obs),slope_other=slope,
                         observed_ppm=measured*1e6,predicted_ppm=expected*1e6,
                         residual_ppm=(measured-expected)*1e6,resid_low_ppm=lo,resid_high_ppm=hi))
    with (HERE/"cross_quarter_test.csv").open('w') as f:
        writer=csv.DictWriter(f,rows[0],lineterminator='\n');writer.writeheader();writer.writerows(rows)
    for row in rows: print(row)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=HERE/'.cache')
    cache=parser.parse_args().cache;cache.mkdir(exist_ok=True)
    main(cache)
