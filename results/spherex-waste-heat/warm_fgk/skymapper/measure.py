"""SkyMapper DR4 g/r detection photometry at the two red warm-FGK Gaia IDs."""
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import requests

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
TARGET={'FGK18':'3017430205815280384','FGK31':'6359523376149191936'}
SERVICES=[('https://api.skymapper.nci.org.au/public/tap/sync','dr4'),
          ('https://datalab.noirlab.edu/tap/sync','skymapper_dr4')]
MASTER=('object_id,raj2000,dej2000,mean_epoch,g_psf,e_g_psf,r_psf,e_r_psf,'
        'class_star,flags_psf,gaia_dr3_id1,gaia_dr3_dist1,gaia_dr3_dist2,self_dist1')
receipts={}
bad={}


def query(name,adql,required):
    path=HERE/'raw'/(name+'.csv')
    if path.exists():
        raw=path.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=receipts[name]['sha256']: raise ValueError('Cached data checksum mismatch')
        return list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))))
    errors=[]
    for url,schema in SERVICES:
        q=adql.replace('DRTABLE.',schema+'.')
        params=dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC='10000',QUERY=q)
        try:
            res=requests.get(url,params=params,timeout=(12,100))
            res.raise_for_status()
            reader=csv.DictReader(io.StringIO(res.content.decode('utf-8-sig')))
            if not set(required)<=set(reader.fieldnames or []):
                raise ValueError('Unexpected schema '+str(reader.fieldnames)[:200])
            data=list(reader)
            if len(data)>=10000: raise ValueError('Possible truncation')
            receipts[name]=dict(url=res.url,query=q,utc=datetime.now(timezone.utc).isoformat(),
                                http=res.status_code,bytes=len(res.content),sha256=hashlib.sha256(res.content).hexdigest())
            path.write_bytes(res.content)
            print(name,len(data),'rows',url,flush=True)
            return data
        except Exception as exc:
            errors.append(url+' '+str(exc)[:300])
    raise RuntimeError(' | '.join(errors))


def num(row,key):
    return float(row.get(key) or np.nan)


def write(name,rows,keys):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def dist(a,b,ra,de):
    return 3600*math.hypot(((a-ra+180)%360-180)*math.cos(math.radians(de)),b-de)


def measure(label,star):
    sid=star['source_id'];ra=num(star,'ra');de=num(star,'dec')
    def cone(rad):
        return f"1=CONTAINS(POINT('ICRS',raj2000,dej2000),CIRCLE('ICRS',{ra:.9f},{de:.9f},{rad:.8f}))"
    target=query(label+'_master','SELECT '+MASTER+' FROM DRTABLE.master WHERE '+cone(6/3600),['object_id','gaia_dr3_id1'])
    assoc=[s for s in target if str(s['gaia_dr3_id1'])==sid and num(s,'gaia_dr3_dist1')<1.5]
    if len(assoc)!=1:
        return [],[],[dict(label=label,source_id=sid,status='no unique exact Gaia association',all_sky_matches=len(target))]
    anchor=assoc[0];oid=anchor['object_id']
    qr=('SELECT '+MASTER+' FROM DRTABLE.master WHERE '+cone(.12)+
        f" AND r_psf BETWEEN {num(anchor,'r_psf')-1.5:.4f} AND {num(anchor,'r_psf')+1.5:.4f}")
    allpeers=query(label+'_peers',qr,['object_id','g_psf','r_psf'])
    color=num(anchor,'g_psf')-num(anchor,'r_psf')
    pool=[]
    for s in allpeers:
        col=num(s,'g_psf')-num(s,'r_psf')
        if (s['object_id']==oid or not np.isfinite(col) or abs(col-color)>.6 or num(s,'class_star')<.8
            or not 0<num(s,'e_r_psf')<.09 or not 0<num(s,'e_g_psf')<.09 or int(s['flags_psf'] or '0')&12):
            continue
        score=((num(s,'r_psf')-num(anchor,'r_psf'))/1.5)**2+((col-color)/.6)**2+(
            dist(num(s,'raj2000'),num(s,'dej2000'),ra,de)/432)**2
        pool.append((score,s['object_id']))
    peers=[s[1] for s in sorted(pool)[:16]]
    ids=','.join([oid]+peers)
    q=('SELECT p.object_id,p.image_id,p.ccd,p.filter,p.ra_img,p.decl_img,p.flags,p.nimaflags,'
       'p.img_qual,p.mag_psf,p.e_mag_psf,p.use_in_clipped,i.night_mjd,i.date AS exposure_mjd '
       'FROM DRTABLE.photometry AS p JOIN DRTABLE.images AS i ON i.image_id=p.image_id '
       f"WHERE p.object_id IN ({ids}) AND p.filter IN ('g','r') ORDER BY p.image_id,p.object_id")
    phot=query(label+'_exposures',q,['object_id','image_id','filter','exposure_mjd','night_mjd','mag_psf'])
    valid=[]
    for row in phot:
        tt=num(row,'exposure_mjd')
        if not (num(row,'flags')==0 and num(row,'nimaflags')==0 and 1<=num(row,'img_qual')<=2
                and 0<num(row,'e_mag_psf')<=.15 and 8<num(row,'mag_psf')<30 and np.isfinite(tt)):continue
        if row['object_id']==oid:
            years=(tt-57388.5)/365.25
            a=ra+num(star,'pmra')*years/3.6e6/math.cos(math.radians(de));b=de+num(star,'pmdec')*years/3.6e6
            if dist(num(row,'ra_img'),num(row,'decl_img'),a,b)>1.2:continue
        valid.append(row)
    summaries=[];epochs=[];nights=[]
    for filt in ['g','r']:
        allp=[r for r in phot if r['filter'].strip()==filt and r['object_id']==oid]
        pp=[r for r in valid if r['filter'].strip()==filt]
        tars=[r for r in pp if r['object_id']==oid]
        grouped=defaultdict(list)
        for r in pp:
            if r['object_id'] in peers:grouped[r['object_id']].append(num(r,'mag_psf'))
        med={k:np.median(v) for k,v in grouped.items() if len(v)>=3}
        matched=defaultdict(list)
        for row in pp:
            if row['object_id'] in med:
                matched[row['image_id'],row['ccd']].append(num(row,'mag_psf')-med[row['object_id']])
        epoch=[]
        for r in tars:
            off=matched[r['image_id'],r['ccd']];n=len(off)
            delta=float(np.median(off)) if n>=3 else 0.
            off_err=1.253*1.4826*np.median(abs(np.array(off)-delta))/math.sqrt(n) if n>=3 else 0.
            epoch.append(dict(label=label,source_id=sid,object_id=oid,band=filt,image_id=r['image_id'],
                      ccd=r['ccd'],night_mjd=r['night_mjd'],mjd=num(r,'exposure_mjd'),
                      raw_mag=num(r,'mag_psf'),raw_error=num(r,'e_mag_psf'),comparison_stars=n,
                      image_offset=delta,corrected_mag=num(r,'mag_psf')-delta,
                      error=float(np.sqrt(num(r,'e_mag_psf')**2+off_err**2+.02**2)),
                      master_kept=r['use_in_clipped']))
        epochs+=epoch
        groups=defaultdict(list)
        for ep in epoch:groups[ep['night_mjd']].append(ep)
        local=[]
        for night,points in sorted(groups.items()):
            local.append(dict(label=label,source_id=sid,band=filt,night_mjd=int(night),
                         mean_mjd=float(np.mean([p['mjd'] for p in points])),detections=len(points),
                         differential_exposures=sum(p['comparison_stars']>=3 for p in points),
                         differential_or_calibrated_mag=float(np.median([p['corrected_mag'] for p in points])),
                         assigned_error=float(np.max([p['error'] for p in points]))))
        report=dict(label=label,source_id=sid,object_id=oid,band=filt,raw_detections=len(allp),
                    quality_detections=len(epoch),nights=len(local),master_neighbor_bias_bit=bool(int(anchor['flags_psf'] or '0') & {'g':8,'r':4}[filt]),
                    reference_stars=peers,gaia_separation_arcsec=num(anchor,'gaia_dr3_dist1'),
                    secondary_gaia_separation_arcsec=anchor.get('gaia_dr3_dist2',''))
        if len(local)<3:
            report['status']='fewer than 3 observed nights'
            for row in local: row.update(residual_mag='',sigma='',upper_3sigma_mag='')
        else:
            mags=np.array([p['differential_or_calibrated_mag'] for p in local]);base=np.median(mags)
            scatter=1.4826*np.median(abs(mags-base))
            for p in local:
                z=p['differential_or_calibrated_mag']-base
                e=max(scatter,p['assigned_error'])*math.sqrt(1+1.253**2/len(mags))
                p.update(residual_mag=float(z),sigma=float(e),upper_3sigma_mag=float(z+3*e))
            report.update(status='sampled-night test',median_mag=float(base),
                          robust_scatter_mag=float(scatter),first_mjd=local[0]['mean_mjd'],last_mjd=local[-1]['mean_mjd'],
                          max_sampled_gap_days=float(np.max(np.diff([p['mean_mjd'] for p in local]))),
                          n_differential_exposures=sum(p['comparison_stars']>=3 for p in epoch),
                          faintest_residual_mag=max(p['residual_mag'] for p in local),
                          deepest_pointwise_3sigma_upper_mag=max(p['upper_3sigma_mag'] for p in local),
                          nights_excluding_0p2mag_3sigma=sum(p['upper_3sigma_mag']<.2 for p in local),
                          nights_at_least_0p2mag_and_3sigma=sum(p['residual_mag']>=.2 and p['residual_mag']>=3*p['sigma'] for p in local))
        nights+=local;summaries.append(report)
    return epochs,nights,summaries


def main():
    (HERE/'raw').mkdir(exist_ok=True)
    if (HERE/'receipts.json').exists():receipts.update(json.loads((HERE/'receipts.json').read_text())['requests'])
    astrom=ROOT/'literature/fgk_astrometry.csv'
    stars={s['source_id']:s for s in csv.DictReader(astrom.open())}
    er=[];nr=[];reports=[]
    for label,sid in TARGET.items():
        try:
            e,n,s=measure(label,stars[sid]);er+=e;nr+=n;reports+=s
        except Exception as exc:
            bad[label]=str(exc);print('ERROR',label,exc,flush=True)
    if er:write('epochs.csv',er,list(er[0]))
    if nr:write('nights.csv',nr,list(nr[0]))
    (HERE/'summary.json').write_text(json.dumps(dict(complete=not bad,errors=bad,series=reports),indent=2)+'\n')
    (HERE/'receipts.json').write_text(json.dumps(dict(astrometry_sha256=hashlib.sha256(astrom.read_bytes()).hexdigest(),
                                       requests=receipts,errors=bad),indent=2)+'\n')
    print(json.dumps(reports,indent=2))


if __name__=='__main__':
    main()
