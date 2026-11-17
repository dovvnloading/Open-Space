"""Select follow-up from the observed multi-survey screens, preserving all distance tests."""
import csv
import glob
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from screen_desi import cartesian


def rows(path):
    with open(path) as f: return list(csv.DictReader(f))

def save(path,table):
    keys=list(dict.fromkeys(k for r in table for k in r))
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(table)

def number(v):
    try:return float(v)
    except (TypeError,ValueError):return np.nan


def build(path):
    gaia={r['source_id']:r for p in sorted(glob.glob('gaia/rows_*.csv')) for r in rows(p)}
    bj={r['source_id']:r for r in rows('bj.csv')};scr={r['SOURCE_ID']:r for r in rows('desi_screen.csv')}
    los={r['source_id']:r for r in rows('desi/los/leads.csv')}
    with fits.open(path) as f:
        a=f[1].data
        for r in a[np.isin(a['SOURCE_ID'],[int(k) for k in los])]:
            scr[str(r['SOURCE_ID'])]={k:r[k] for k in a.names}
    for r in rows('desi/los/gaia.csv'):gaia[r['source_id']]=r
    def source(sid,rv,erv,channel,sp=np.nan,esp=np.nan):
        g=gaia[sid];s=scr[sid]
        return dict(source_id=sid,targetid=s['TARGETID'],channel=channel,rv=rv,e_rv=erv,spec_d=s['DIST'],spec_e=s['DISTERR'],zp=s['PARALLAX_ZPC'],sp_d=sp,sp_e=esp,**{k:v for k,v in bj.get(sid,{}).items() if k!='source_id'})
    queue=[];checks=rows('desi_checks.csv')
    for r in checks:
        sid=r['source_id'];g=gaia[sid]; b=bj.get(sid,{})
        d=number(b.get('r_med_photogeo'))/1000
        v=np.linalg.norm(cartesian(number(g['ra']),number(g['dec']),d,number(g['pmra']),number(g['pmdec']),number(r['rv']))[1]) if d>0 else np.nan
        r['bj_photogeo_speed']=v
        if r['rv_ok']=='1' and (r['nearby']=='1' or v>=550):
            queue.append(source(sid,r['rv'],r['rv_err'],'nearby' if r['nearby']=='1' else 'fast_spec',r['spdist_d'],r['spdist_sigma']))
    for sid,r in los.items():
        if sid in scr and number(gaia[sid]['ruwe'])<1.4:
            queue.append(source(sid,r['vrad'],r['vrad_err'],'radial'))
    save('screen_review.csv',checks);save('desi_queue.csv',queue)
    save('gaia_queue.csv',[gaia[r['source_id']] for r in queue])

if __name__=='__main__':build(sys.argv[1])
