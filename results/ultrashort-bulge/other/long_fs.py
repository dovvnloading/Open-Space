"""Independent light curves for KMT-2024-BLG-1516, a weeks-long finite-source flag.

Use the common continuous fitter for each stream separately. Compare the KMT
geometry with a point source at KMT's best t0,tE,u0 and score unrelated seasons.
"""
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
# Load this directory's survey helpers before fit imports OGLE's lens module.
from crossmatch import coords, extern
from fetch import SHA, get, write_csv
from run_grid import read_curve
from refine import fit, local_fit
from concurrent.futures import ProcessPoolExecutor

KMT = 'KMT-2024-BLG-1516'
PEAK = 2460507.911
EVENTS = ['MOA-2024-BLG-113','PRIME-2024-BLG-109','OGLE-2024-BLG-0859']
URL = 'https://www.astrouw.edu.pl/ogle/ogle4/ews/2024/blg-0859/phot.dat'


def jobs():
    phot = {r['event']:r for r in csv.DictReader((HERE/'photometry.csv').open())}
    for event in EVENTS[:2]:
        r = phot[event]
        yield (event,Path('.survey-cache/curves')/(event+'.npz'),float(r['page_t0_jd']),15.),r['phot_url'],r['input_sha256']
    p = Path('.ogle-cache/raw/2024/blg-0859.dat')
    if not p.exists():
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_bytes(get(URL))
    yield (EVENTS[2],p,2460508.281,12.3),URL,SHA(p.read_bytes())


def analyze(item):
    job,url,checksum = item
    event = job[0]
    # Published local solution seeds derived from six KMT 2024 I-band streams;
    # each other survey independently optimizes its geometric and flux parameters.
    seed = {'ps_t0':2460508.185587,'ps_te':13.544575,'ps_u0':0.956148,
            'fs_t0':2460507.911049,'fs_te':11.496478,'fs_u0':1.597735,'fs_rho':1.843634}
    out = local_fit(job,seed,False)
    t,f,e,g,diff=read_curve(job)
    t,f,e,g=fit.clean(t,f,e,g)
    order=np.argsort(t)
    t,f,e=t[order],f[order],e[order]
    if diff:
        e=fit.reported_error_floor(e,np.zeros(len(t),int))
    off=abs(t-PEAK)>max(5,min(2*job[3],100))
    med=np.median(f[off])
    e=e*max(1.,np.median(abs((f[off]-med)/e[off]))/.67448975)
    cap=np.inf if diff else 1.15*np.median(f)
    w=1/e**2
    def chi(a):
        mean=np.average(a,weights=w)
        mflux=np.average(f,weights=w)
        centered=a-mean
        scale=np.clip(np.sum(w*centered*(f-mflux))/np.sum(w*centered**2),0,cap)
        return np.sum(w*(f-mflux-scale*centered)**2)
    tau=(t-seed['ps_t0'])/seed['ps_te']
    ps=fit.lens.magnification(tau**2+seed['ps_u0']**2)-1
    vb=fit.engine()
    t_off=t-seed['fs_t0']
    fs=np.asarray(vb.ESPLLightCurve([np.log(seed['fs_u0']),np.log(seed['fs_te']),0.,np.log(seed['fs_rho'])],t_off.tolist())[0])-1
    out.update(input_sha256=checksum,input_url=url,kmt_fixed_ps_chi2=round(chi(ps),3),
               kmt_fixed_fs_chi2=round(chi(fs),3),kmt_fixed_ps_minus_fs=round(chi(ps)-chi(fs),3))
    baseline=[]
    year=np.array([datetime.fromtimestamp((a-2440587.5)*86400,timezone.utc).year for a in t])
    for yr in np.unique(year):
        mask=(year==yr)&(abs(t-PEAK)>65)
        if mask.sum()<10: continue
        tt,ff,ee=t[mask],f[mask],e[mask]
        # Require a full night with at least three points and measure night-to-night
        # scatter from actual medians, including correlated and crowding noise.
        nights=np.unique(np.floor(tt))
        bins=[(n,np.median(ff[np.floor(tt)==n])) for n in nights if (np.floor(tt)==n).sum()>=3]
        if len(bins)<4:continue
        ij,m=np.array(bins).T
        med=np.median(m);s=1.4826*np.median(abs(m-med))
        if s<=0:continue
        baseline.append(dict(event=event,year=int(yr),outside_65d_epochs=int(mask.sum()),n_nights_three_epochs=len(m),
                             median=round(med,5),night_mad_sigma=round(s,5),brightest_night_jd=int(ij[np.argmax(m)]),
                             brightest_night_sig=round((max(m)-med)/s,3),faintest_night_sig=round((min(m)-med)/s,3)))
    return out,baseline


def catalogs():
    rows=[]; objects=[]
    for event,ra,dec in [(KMT,'17:58:38.76','-27:24:57.82'),('KMT-2024-BLG-3349','18:04:47.04','-27:07:24.89')]:
        x,y=coords(ra,dec)
        for name in ['Gaia DR3','VVV DR2']:
            receipt,found=extern((event,x,y,name))
            rows.append(receipt)
            objects.extend(found)
    write_csv(HERE/'long_fs_catalog_queries.csv',rows)
    write_csv(HERE/'long_fs_catalog_neighbours.csv',objects)


if __name__ == '__main__':
    with ProcessPoolExecutor(3) as proc:
        records=list(proc.map(analyze,jobs()))
    write_csv(HERE/'long_fs_independent.csv',[dict(row,source_fractions=json.dumps(row['source_fractions'])) for row,_ in records])
    write_csv(HERE/'long_fs_baselines.csv',[r for _,rr in records for r in rr])
    catalogs()
