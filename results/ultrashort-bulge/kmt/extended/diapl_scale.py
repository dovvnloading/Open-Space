"""Transfer the released pySIS I scale to DIAPL differential counts on chip BLG31M."""
import argparse
import csv
import hashlib
import io
import json
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd

from sightlines import HERE,KMT,PORTAL,read,catalog,urlencode

# Same subfield first, then high signal controls on the same physical CCD.
CAL=['2412','0273','1694','1096','3131']


def evaluate(item):
    n,cache=item;name='KMT-2024-BLG-'+n;r=catalog(name)
    c=cache/name;c.mkdir(parents=True,exist_ok=True)
    lcid='BLG'+r['Field'][3:]+'.'+r['Star ID']
    d,rec=read(c,'cal-pos',PORTAL+urlencode(dict(action='download_lc',lcid=lcid)))
    a,ar=read(c,'cal-archive',f'https://kmtnet.kasi.re.kr/ulens/event/2024/data/KB24{n}/pysis/pysis.tar.gz')
    with (KMT/'pySIS-2024.tsv').open() as f:
        frozen=next(z for z in csv.DictReader(f,delimiter='\t') if z['event']==name)
    if hashlib.sha256(a).hexdigest()!=frozen['sha256']:raise ValueError(name)
    frame=pd.read_csv(io.BytesIO(d));t0,te,u0=[float(r[k]) for k in ('t_0','t_E','u_0')]
    data=[];sites={'A':'SSO','S':'SAO','C':'BLG'};rng=np.random.default_rng(6328)
    with tarfile.open(fileobj=io.BytesIO(a)) as arx:
        for m in arx:
            if not m.name.endswith('31_I.pysis'):continue
            s=sites[m.name[3]];aa=np.loadtxt(arx.extractfile(m),usecols=(0,3,4))
            p=pd.DataFrame(aa,columns=['t','mag','magerr']);p.t-=2400000
            p=p[(p.magerr>0)&(p.magerr<.5)&(p.mag>10)&(p.mag<26)].sort_values('t')
            q=frame[(frame.site==s)&(frame.year==2024)&(frame.err>0)&(frame.chi2>0)&(frame.see>0)].sort_values('date')
            # Match repeated exposures within 2.6 minutes. Collapse them by integer night.
            w=pd.merge_asof(p,q,left_on='t',right_on='date',direction='nearest',tolerance=.0018).dropna(subset=['date'])
            w['f20']=10**(-.4*(w.mag-20));w['night']=np.floor(w.t)
            z=w.groupby('night').agg(t=('t','median'),flux=('flux','median'),f20=('f20','median'))
            t,df,y=z.to_numpy().T if len(z) else (np.array([]),)*3
            peaks=abs(t-t0)<te*max(.5,u0);off=abs(t-t0)>2.5*te
            out=dict(site=s,matched_exposures=len(w),n_peak=int(sum(peaks)),n_base=int(sum(off)))
            if sum(peaks)>=3 and sum(off)>=8:
                on,base=np.flatnonzero(peaks),np.flatnonzero(off)
                def contrasts(ii,jj):
                    return np.median(df[ii])-np.median(df[jj]),np.median(y[ii])-np.median(y[jj])
                di,py=contrasts(on,base)
                bs=np.array([contrasts(rng.choice(on,len(on)),rng.choice(base,len(base))) for _ in range(700)])
                sign=np.mean((bs[:,0]>0)&(bs[:,1]>0))
                out.update(delta_diapl=round(float(di),2),delta_py_f20=round(float(py),4),bootstrap_fraction_positive=round(float(sign),3))
                if sign>.95 and di>0 and py>0:
                    ratio=bs[(bs[:,0]>0)&(bs[:,1]>0)]
                    out.update(counts_at_I20=round(di/py,2),bootstrap_p16_p84=np.percentile(ratio[:,0]/ratio[:,1],[16,84]).round(2).tolist())
            data.append(out)
    return dict(event=name,field=r['Field'],catalog_u0=u0,lc_receipt=rec,archive_receipt=ar,streams=data)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--cache',type=Path,required=True)
    c=ap.parse_args().cache
    with ThreadPoolExecutor(2) as pool:rows=list(pool.map(evaluate,[(n,c) for n in CAL]))
    transfer=[]
    prof=json.loads((HERE/'tier2-results.json').read_text())['events'][0]
    peak=prof['tests'][0]['position_fits'][0]['sites'];R=15.;amax=np.sqrt(1.+4/R**2)-1
    for site,data in peak.items():
        z=[s['bootstrap_p16_p84'][1] for r in rows for s in r['streams'] if s['site']==site and 'counts_at_I20' in s]
        if z:
            transfer.append(dict(site=site,fixed_shape_peak_counts=data[0],uniform_center_Aminus1=round(amax,6),
                                 minimum_source_counts_at_centroid=round(data[0]/amax,1),
                                 largest_calibrator_p84_counts_I20=max(z),required_ratio_to_comparison_I20=round(data[0]/amax/max(z),2)))
    (HERE/'diapl-scale.json').write_text(json.dumps(dict(calibrators=rows,transfer=transfer),indent=2)+'\n')


if __name__=='__main__':main()
