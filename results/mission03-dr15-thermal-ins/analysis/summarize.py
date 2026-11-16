"""Summarize catalogue neighbors; keep survey coverage and flux-limit analysis separate."""
import argparse
import csv
import json
import math
from pathlib import Path
from counterparts import separation

def f(x):
    try:
        y=float(x);return y if math.isfinite(y) else None
    except (ValueError,TypeError): return None

def get_nearest(survey,rows,ra,dec):
    def measure(row):
        if survey=='ps1':
            if row.get('primaryDetection')!='1':return
            a,b=f(row.get('raStack')),f(row.get('decStack'))
        else:a,b=f(row.get('ra')),f(row.get('dec'))
        if a is None or b is None: return
        sep=separation(ra,dec,a,b)
        if survey=='gaia':
            # Closest point of linear proper-motion track during the release interval.
            # This is deliberately conservative until an observation epoch is known.
            east=((a-ra+180)%360-180)*3600*math.cos(math.radians(dec))
            north=(b-dec)*3600
            ve,vn=(f(row.get('pmra')) or 0)/1000,(f(row.get('pmdec')) or 0)/1000
            if ve*ve+vn*vn:
                dt=max(-16,min(8.80,-(east*ve+north*vn)/(ve*ve+vn*vn)))
                sep=math.hypot(east+ve*dt,north+vn*dt)
            obj,mag=row['source_id'],f(row.get('phot_g_mean_mag'))
        elif survey=='legacy':
            flux=f(row.get('flux_r'))
            mag=(22.5-2.5*math.log10(flux)) if flux and flux>0 else None
            obj=row['ls_id']
        elif survey=='ps1':
            mag=f(row.get('rPSFMag'));mag=mag if mag and mag>0 else None;obj=row['objID']
        else:obj,mag=row['designation'],f(row.get('w1mpro'))
        return sep,obj,mag
    candidates=[v for row in rows if (v:=measure(row)) is not None]
    return min(candidates,key=lambda item: item[0],default=None)

def summarize(targets,cache):
    for item in csv.DictReader(open(targets)):
        r={k.lower():v for k,v in item.items()};sid=r['srcid'];ra,dec=map(float,(r['ra'],r['dec']))
        sigma=math.hypot(float(r.get('radec_err') or r.get('error_radius'))/math.sqrt(2),0.88)
        pos95=sigma*math.sqrt(-2*math.log(.05));pos99=sigma*math.sqrt(-2*math.log(.01))
        out=dict(srcid=sid,name=(r.get('iauname') or r.get('name')),ra=ra,dec=dec,r95_arcsec=round(pos95,3),r99_arcsec=round(pos99,3),fsoft=r.get('f_soft') or r.get('fsoft'))
        for cat in ('gaia','allwise','legacy','ps1'):
            path=cache/f'{sid}-{cat}.json'
            data=json.loads(path.read_text()) if path.exists() else {}
            out[cat+'_status']=data.get('status','not_queried')
            near=get_nearest(cat,data.get('rows',[]),ra,dec)
            out[cat+'_sep']=round(near[0],3) if near else ''
            out[cat+'_id']=near[1] if near else ''
            out[cat+'_mag']=round(near[2],3) if near and near[2] is not None else ''
            out[cat+'_ncone']=len(data['rows']) if data.get('status')=='ok' else ''
        yield out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('targets');p.add_argument('cache',type=Path);p.add_argument('output');a=p.parse_args()
    rows=list(summarize(a.targets,a.cache))
    with open(a.output,'w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
