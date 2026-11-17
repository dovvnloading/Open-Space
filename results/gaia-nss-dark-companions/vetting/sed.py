"""Fit broadband composite SEDs with a scaled primary and empirical dwarf companions."""
import argparse
import csv
import hashlib
import json
import urllib.request
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

HERE=Path(__file__).parent
URL='https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt'
BANDS=['B','V','G','J','H','Ks','W1']
# Fixed effective-band A_lambda/A_V, approximate R_V=3.1. Assessed by free A_V.
EXT=np.array([1.32,1,.84,.282,.19,.114,.06])
TARGETS={'3378588057203660160':(1.77,7500,10400),'1830677037678583680':(4.05,8500,24500)}


def templates(path):
    if not path.exists(): path.write_bytes(urllib.request.urlopen(URL,timeout=30).read())
    lines=path.read_text().splitlines()
    head=next(l.split() for l in lines if l.startswith('#SpT'))
    vals=[]
    for l in lines:
        cols=l.split()
        if len(cols)!=len(head) or l.startswith('#'):continue
        d=dict(zip(head,cols))
        def num(k):
            try:return float(d[k])
            except ValueError:return np.nan
        mg=num('M_G')
        if not np.isfinite(mg):mg=num('Mv')+num('G-V')
        vals.append([num('Msun'),num('Teff'),num('Mv')+num('B-V'),num('Mv'),mg,
                     num('M_J'),num('M_Ks')+num('H-Ks'),num('M_Ks'),num('M_Ks')-num('Ks-W1')])
    return np.asarray(vals)


def mags(tab,x,column):
    out=[]
    for k in range(len(BANDS)):
        d=tab[np.isfinite(tab[:,column])&np.isfinite(tab[:,2+k])]
        d=d[np.argsort(d[:,column])]
        out.append(np.interp(x,d[:,column],d[:,2+k],left=np.nan,right=np.nan))
    return np.array(out)


def read_one(sid,key):
    with (HERE/'photometry'/(sid+'_'+key+'.csv')).open() as f: rows=list(csv.DictReader(f))
    if len(rows)!=1:raise ValueError('Ambiguous or missing match')
    return rows[0]


def run(cache):
    tmp=cache/'sequence.txt'
    tab=templates(tmp); records=[]
    for sid,(mass,tmin,tmax) in TARGETS.items():
        phot=read_one(sid,'optical'); tm=read_one(sid,'2mass'); wise=read_one(sid,'allwise')
        g=next(r for r in csv.DictReader((HERE.parent/'gaia/source_diagnostics.csv').open()) if r['source_id']==sid)
        o=next(r for r in csv.DictReader((HERE.parent/'gaia/mass_screen.csv').open()) if r['source_id']==sid)
        plx=float(o['parallax'] if o['parallax']!='null' else g['parallax'])
        dm=10-5*np.log10(plx)
        observed=np.array([float(a) for a in [phot['B'],phot['V'],g['phot_g_mean_mag'],tm['j_m'],tm['h_m'],tm['k_m'],wise['w1mpro']]])
        err=np.array([.06,.06,.035,*[np.hypot(float(tm[k]),.035) for k in ('j_cmsig','h_cmsig','k_cmsig')],np.hypot(float(wise['w1sigmpro']),.035)])
        keep=np.ones(len(BANDS),bool)
        if wise['cc_flags'][0]!='0' or wise['ph_qual'][0]!='A':keep[-1]=False
        for model,n in [('single',0),('one_dwarf',1),('two_equal_dwarfs',2)]:
            second=10**(-.4*mags(tab,mass/n,0))*n if n else np.zeros(len(BANDS))
            def pred(theta):
                teff,scale,av=theta
                fprimary=10**(-.4*(mags(tab,teff,1)+scale))
                return dm+EXT*av-2.5*np.log10(fprimary+second)
            def residual(theta): return ((observed-pred(theta))/err)[keep]
            # Try several temperatures to avoid local optima at empirical sequence knots.
            starts=np.linspace(tmin,tmax,9)[1:-1]
            out=min((least_squares(residual,[t,-.3,.2],bounds=([tmin,-6,0],[tmax,4,2.5]))
                     for t in starts),key=lambda r:r.cost)
            frac=second/(second+10**(-.4*(mags(tab,out.x[0],1)+out.x[1])))
            records.append(dict(source_id=sid,model=model,m2_total=mass,plx_mas=plx,chi2=2*out.cost,
                 dof=int(keep.sum()-3),teff_K=out.x[0],primary_mag_offset=out.x[1],av=out.x[2],
                 bands=np.array(BANDS)[keep].tolist(),residual_mags=(observed-pred(out.x))[keep].tolist(),
                 fractions_by_band=dict(zip(BANDS,frac))))
    return dict(dwarf_sequence_url=URL,sha256=hashlib.sha256(tmp.read_bytes()).hexdigest(),fits=records)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    a=p.parse_args();a.cache.mkdir(parents=True,exist_ok=True)
    (HERE/'sed_fits.json').write_text(json.dumps(run(a.cache),indent=2)+'\n')
