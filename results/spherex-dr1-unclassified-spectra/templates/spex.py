"""Compare coadd continua to independently observed stellar spectra, with telluric gaps masked."""
import argparse
import csv
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote
from urllib.request import urlopen

import numpy as np

HERE=Path(__file__).resolve().parent
ARCHIVE='https://irsa.ipac.caltech.edu/data/IRTF_Extra/SpectralLibrary/spectra/spectrallibrary/'
STARS=('F5V_HD218804 G0V_HD109358 G5V_HD165185 K0V_HD145675 K3V_HD219134 '
       'K7V_HD237903 M0V_HD19305 M3V_Gl388 M5V_Gl51 M8V_Gl752B K0III_HD100006 '
       'K5III_HD181596 M5III_HD175865 L3.5_2MASSJ0036+1821 T4.5_2MASSJ0559-1404').split()
TELLURIC=[(1.35,1.45),(1.80,1.95),(2.50,2.85),(4.15,4.50)]

def write(file,rows):
    with file.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def fetch(grid):
    def one(star):
        u=ARCHIVE+quote(star+'.txt')
        with urlopen(u,timeout=100) as r:raw=r.read()
        wl,f,err=np.loadtxt(io.BytesIO(raw),unpack=True)
        # Normalized Fnu is proportional to lambda squared times Flambda.
        good=(wl>0)&(f>0)&(err>=0)&np.isfinite(f)&(f!=-999)
        vals=[]
        for a,b in grid:
            use=(wl>=a)&(wl<=b)&good
            val=np.nan
            if use.sum()>=6 and not any(a<h and b>l for l,h in TELLURIC):
                w=wl[use]
                if (w[-1]-w[0])>=.85*(b-a) and np.diff(w).max()<.2*(b-a):
                    val=np.trapezoid(f[use]*w*w,w)/(w[-1]-w[0])
            vals.append(val)
        vals=np.array(vals);vals/=np.nanmedian(vals)
        return vals,dict(template=star,url=u,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
    with ThreadPoolExecutor(5) as pool:all_stars=list(pool.map(one,STARS))
    (HERE/'spex_inputs.json').write_text(json.dumps([x[1] for x in all_stars],indent=2)+'\n')
    rows=[]
    for i,(a,b) in enumerate(grid):
        rows.append(dict(channel=i+1,**{s:format(v[i],'.8g') if np.isfinite(v[i]) else '' for s,(v,_) in zip(STARS,all_stars)}))
    write(HERE/'spex_grid.csv',rows)

def compare(field):
    import pandas as pd
    specs=pd.read_csv(field/'spectra.csv')
    stars=pd.read_csv(HERE/'spex_grid.csv').set_index('channel')
    rows=[]
    for src,g in specs.groupby('source_id'):
        g=g.sort_values('channel')
        y=g.flux.to_numpy();sig=g.sigma.to_numpy()
        good=np.isfinite(y+sig)&(sig>0)&(g.flagged.to_numpy()==0)
        # Fixed comparison support across all templates for this source. If scarce,
        # leave the source indeterminate rather than rewarding templates with gaps.
        common=good&np.isfinite(stars.to_numpy()).all(axis=1)
        if common.sum()<8:continue
        errs=[]
        for name in stars:
            v=stars[name].to_numpy()[common]
            scale=np.median(y[common]/v)
            if scale<=0:continue
            ratio=np.abs(y[common]-scale*v)/np.maximum(np.abs(scale*v),sig[common])
            errs.append((np.median(ratio),name,scale))
        if not errs:continue
        loss,name,scale=min(errs)
        rows.append(dict(source_id=src,template=name,n_common=common.sum(),median_abs_fraction=round(loss,6),
                         scale=round(scale,8)))
    if rows:write(field/'catalogs'/'spex_fits.csv',rows)
    return len(rows)

if __name__=='__main__':
    import pandas as pd
    p=argparse.ArgumentParser();p.add_argument('--fetch',action='store_true');p.add_argument('--wide',action='store_true');a=p.parse_args()
    field=HERE.parent/'wide' if a.wide else HERE.parent
    if a.fetch:
        g=pd.read_csv(field/'spectra.csv').drop_duplicates('channel').sort_values('channel')
        fetch(g[['wave_lo_um','wave_hi_um']].to_numpy())
    print('empirical comparisons', compare(field))
