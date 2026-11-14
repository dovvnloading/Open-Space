"""Independent historic fixed-window check of OGLE-LMC-ECL-13803."""
import gzip
import json
from pathlib import Path

import numpy as np

root=Path(__file__).parent
period,epoch=1.9964371,2460001.023
rng=np.random.default_rng(13803)
def compare(t,m,e,maxerr):
    usable=np.isfinite(t+m+e)&(e>0)&(e<=maxerr)
    t,m,e=t[usable]+2450000,m[usable],e[usable]
    ph=np.abs(((t-epoch)/period+.5)%1-.5)
    results={}
    for w in (.03,.05):
        days=[np.array([np.median(m[mask&(t.astype(int)==n)]) for n in np.unique(t[mask].astype(int))])
              for mask in [(ph>.16)&(ph<.34),ph<w,ph>.5-w]]
        b,*core=days
        row={'baseline_nights':len(b)}
        for k,a in zip(['phase_zero','phase_half'],core):
            row[k]={'nights':len(a)}
            if len(a)>3:
                diffs=np.median(rng.choice(a,(5000,len(a)),replace=True),axis=1)-np.median(rng.choice(b,(5000,len(b)),replace=True),axis=1)
                row[k].update(depth_magnitude=np.median(a)-np.median(b),bootstrap95=list(np.percentile(diffs,[2.5,97.5])))
        results[str(w)]=row
    return {'used':int(usable.sum()),'hjd_range':[float(t.min()),float(t.max())],'by_core_half_width':results}
with gzip.open(root/'inputs'/'eros-target.tsv.gz', 'rt') as source:
    tab = [line for line in source if line.startswith('lm0010m13681')]
data=np.array([[float(v.strip() or 'nan') for v in line.rstrip('\n').split('\t')[1:]] for line in tab])
assert len(data)==654
t,r,er,b,eb=data.T
out={'EROS_association':'lm0010m13681, 0.420 arcsec, 654 paired exposures. Time released as HJD minus 2450000. Error cut <=0.5 mag. Separate passbands. No zero point or dilution transfer across surveys.', 'OGLE_III':compare(*np.loadtxt(root/'inputs'/'ogle3-I.dat').T,.1),'EROS_R':compare(t,r,er,.5),'EROS_B':compare(t,b,eb,.5)}
(root/'lmc-early.json').write_text(json.dumps(out,indent=2)+'\n')
