"""Apply a second spectrophotometric estimator and retain diagnostics on every screen hit."""
import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from astropy.io import fits
from screen_desi import cartesian

SHA='ef05275301814b6a1065cff99217ed3a2628472774a81750ecde5db5cb4526d5'

def run(path):
    with open(path,'rb') as f:
        assert hashlib.file_digest(f,'sha256').hexdigest()==SHA
    screen=list(csv.DictReader(open('desi_screen.csv')))
    ids={int(r['SOURCE_ID']) for r in screen}
    by_id=defaultdict(list)
    with fits.open(path, memmap=True) as f:
        data=f[1].data
        for off in range(0,len(data),100000):
            a=data[off:off+100000]; portion=a[np.isin(a['SOURCE_ID'],list(ids))]
            for r in portion:
                by_id[int(r['SOURCE_ID'])].append({k:r[k] for k in ('TARGETID','FLAG_GOOD','QFLAG_INPUT','G0','MG_16','MG_50','MG_84','LOGG_COR','TEFF_COR')})
    coadds=defaultdict(list)
    for r in csv.DictReader(open('desi/queries/coadds.csv')):
        coadds[r['targetid']].append(r)
    rows=[]
    for s in screen:
        sid=int(s['SOURCE_ID']);rv=float(s['VRAD']); p=float(s['PARALLAX'])-float(s['PARALLAX_ZPC']);ep=float(s['PARALLAX_ERR'])
        d=float(s['DIST']); ed=float(s['DISTERR']); tid=int(s['TARGETID']);hits=[r for r in by_id[sid] if r['TARGETID']==tid]
        # A repeated SPDist row can represent another spectrum. Use the median eligible distance.
        params=[x for x in hits if x['FLAG_GOOD'] and x['QFLAG_INPUT'] and np.isfinite(x['MG_50']+x['G0']+x['MG_16']+x['MG_84'])]
        ds=np.nan; dsig=np.nan; vs=np.nan
        if params:
            ds=np.median([10**(.2*(x['G0']-x['MG_50'])-2) for x in params]);dsig=ds*np.log(10)/5*np.median([np.hypot(.5*(x['MG_84']-x['MG_16']),.167) for x in params]); vs=np.linalg.norm(cartesian(float(s['RA']),float(s['DEC']),ds,float(s['PMRA']),float(s['PMDEC']),rv)[1])
        cands=coadds[str(tid)]
        match=min(cands,key=lambda x:abs(float(x['vrad'])-rv)) if cands else None
        r=match or {}; ok=bool(match and r['source_id']==str(sid) and abs(float(r['vrad'])-rv)<1 and int(r['rvs_warn'])==0 and r['success'] in ('1','t','True') and r['rr_spectype']=='STAR' and 0<float(r['vrad_err'])<20 and abs(rv-299792.458*float(r['rr_z']))<60)
        plx_agree=abs(p-1/d)/np.hypot(ep,ed/d**2)
        rows.append(dict(source_id=sid,targetid=tid,nearby=s['nearby'],fast_spec=s['fast_spec'],spec_v=s['v_spec'],inverse_v=s['v_plx'],spdist_d=ds,spdist_sigma=dsig,spdist_v=vs,spec_d=d,spec_sigma=ed,plx_discrepancy_sigma=plx_agree,rv=rv,rv_err=r.get('vrad_err',''),rvs_warn=r.get('rvs_warn',''),rr_type=r.get('rr_spectype',''),rr_rv=299792.458*float(r['rr_z']) if match else '',rv_ok=int(ok),mg50_good_count=len(params),spdist_join_rows=len(hits)))
    with open('desi_checks.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

if __name__=='__main__':run(sys.argv[1])
