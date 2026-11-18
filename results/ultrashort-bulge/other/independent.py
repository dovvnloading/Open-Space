"""Retrieve independent KMT DIAPL field stars at two MOA spike positions and check SkyBoT.

Restrict the joint 0211 shape comparison to the nearest star in each overlapping field.
Two windows in the same field are not independent measurements of two stars.
"""
import csv
import io
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from xml.etree import ElementTree as ET
import numpy as np
from scipy.optimize import minimize

from fetch import get, OUT, write_csv, SHA

CACHE=Path('.survey-cache')
KMT='https://kmtnet.kasi.re.kr/search/?'
BOT='https://ssp.imcce.fr/webservices/skybot/api/conesearch.php?'
TARGETS={'0093':(273.4972841667,-27.1183508889,2460756.149), '0211':(268.7528345833,-30.5618525833,2460877.14),
         '111':(270.8985416667,-28.3891972222,2460496.29)}


def retrieve(path,url):
    if not path.exists():path.write_bytes(get(url))
    return path.read_bytes()


def fields(target):
    ra,dec,peak=TARGETS[target];result=[]
    url=KMT+urlencode(dict(search_type='radec',ra=ra,dec=dec,ra_type='decimal_degree',arcsec=2))
    raw=retrieve(CACHE/('k'+target),url)
    hits=json.loads(re.search(r'const kmtnetSearchMarkers = (.*?);',raw.decode())[1])
    for hit in hits:
        name=hit['lcid']
        curve_url=KMT+urlencode(dict(action='download_lc',lcid=name))
        curve=retrieve(CACHE/(name+'.csv'),curve_url)
        data=list(csv.DictReader(io.StringIO(curve.decode())))
        dates=np.array([float(r['date'])+2400000 for r in data])
        offsets=dates-peak
        near=abs(offsets)<=.5
        result.append(dict(event=('MOA-2024-BLG-' if target=='111' else 'MOA-2025-BLG-')+target,field_star=name,separation_arcsec=hit['dist'],
                    source_ra=hit['ra'],source_dec=hit['dec'],query=url,query_sha256=SHA(raw),url=curve_url,
                    curve_sha256=SHA(curve),bytes=len(curve),rows=len(data),first_hjd=dates.min(),last_hjd=dates.max(),
                    near_half_day=int(near.sum()),near_0p18d=int((abs(offsets)<=.18).sum()),peak_jd_adopted=peak,
                    nearest_before_hjd=dates[offsets<0].max(),nearest_after_hjd=dates[offsets>0].min()))
    return result


def moving(target):
    ra,dec,t=TARGETS[target]
    url=BOT+urlencode({'-observer':'474','-from':'OpenSpace','RA':f'{ra:.6f}','DEC':f'{dec:.6f}',
                       'TIME':str(t),'SR':'0.0333333','VERB':3,'RESPONSEFORMAT':'votable'})
    path=CACHE/('bot0' if target=='0093' else 'bot'+target)
    try:
        raw=retrieve(path,url);root=ET.fromstring(raw)
        tables=[v for v in root.iter() if v.tag.endswith('TABLE') and v.get('nrows')]
        infos=[(x.get('value'),x.text or '') for x in root.iter() if x.get('name')=='QUERY_STATUS']
        if not tables or any(k=='ERROR' and 'No solar system object was found' not in s for k,s in infos):
            raise ValueError(str(infos))
        return dict(event=('MOA-2024-BLG-' if target=='111' else 'MOA-2025-BLG-')+target,url=url,sha256=SHA(raw),bytes=len(raw),known_objects=int(tables[0].get('nrows')),status='ok',reply=' '.join(s for _,s in infos).strip())
    except Exception as e:return dict(event=('MOA-2024-BLG-' if target=='111' else 'MOA-2025-BLG-')+target,url=url,status='error',error=str(e))


def joint():
    import VBMicrolensing as vb
    eng=vb.VBMicrolensing()
    eng.LoadESPLTable(str(Path(vb.__file__).parent/'data/ESPL.tbl'))
    eng.Tol=.0001
    name='MOA-2025-BLG-0211'
    t0=TARGETS['0211'][2]
    moa=np.load(CACHE/'curves'/f'{name}.npz')['data']
    chunks=[(moa[:,0],moa[:,1],moa[:,2],np.zeros(len(moa)),'MOA')]
    ids=['BLG01M0805.025254','BLG41M0706.026273']
    for i in ids:
        rows=list(csv.DictReader((CACHE/(i+'.csv')).open()))
        for site in ['SSO','SAO','BLG']:
            rr=[r for r in rows if r['site']==site]
            data=np.array([[float(r['date'])+2400000,float(r['flux']),float(r['err']),float(r['chi2'])] for r in rr])
            chunks.append((*data.T,i+'-'+site))
    results=[]
    for quality in ['all','0_lt_kmt_diapl_chi2_le_40','0_lt_kmt_diapl_chi2_le_20']:
      for correction in [-.0058,0.,.0058]:
        arr=[]
        for j,(t,f,e,chi,site) in enumerate(chunks):
            u=t-t0+(correction if j==0 else 0.)
            sel=(abs(u)<10)&(e>0)&np.isfinite(f+e)
            if quality!='all' and j:sel &= (chi>0)&(chi<=int(quality.rsplit('_',1)[1]))
            u,y,sig=u[sel],f[sel],e[sel]
            off=abs(u)>.7
            if len(y)<4 or off.sum()<3:continue
            weight=1/(np.maximum(sig,.2*np.median(sig))*max(1.,np.median(abs((y[off]-np.median(y[off]))/sig[off]))/.67448975))**2
            arr.append((u,y,weight,site))
        t=np.concatenate([v[0] for v in arr]); f=np.concatenate([v[1] for v in arr]);w=np.concatenate([v[2] for v in arr])
        g=np.repeat(np.arange(len(arr)),[len(a[0]) for a in arr])
        def center(v):return v-(np.bincount(g,weights=w*v)/np.bincount(g,weights=w))[g]
        y=center(f)
        def chi(p,kind):
            dt=t-p[0];scale,a=np.exp(p[1:3])
            if kind=='flare': model=np.exp(np.clip(np.where(dt<0,dt/scale,-dt/a),-744,0))
            else:
                sep=np.hypot(dt/scale,a);square=sep**2
                model=4/(sep*np.sqrt(square+4)*(square+2+sep*np.sqrt(square+4)))
                if kind=='fs':
                    rho=np.exp(p[3]);near=sep<max(6*rho,2.)
                    if near.any(): model[near]=np.array(eng.ESPLLightCurve([p[2],p[1],p[0],p[3]],t[near].tolist())[0])-1
            mod=center(model)
            amp=np.maximum(0,np.bincount(g,weights=w*mod*y)/np.maximum(1e-100,np.bincount(g,weights=w*mod**2)))
            return np.dot(w,(y-amp[g]*mod)**2)
        row=dict(event=name,quality=quality,moa_time_shift_d=correction,epochs=len(t),streams=len(arr),start_jd=t0-10,stop_jd=t0+10)
        for kind in ['ps','fs','flare']:
            lim=[(-.12,.12),(np.log(.001),np.log(100)),(np.log(.001),np.log(50))]
            if kind=='fs':lim.append((np.log(.01),np.log(25)))
            pars=[(.007,.15),(.06,.3),(2,.02),(10,.002)]
            if kind=='flare':pars=[(.003,.03),(.015,.025),(.03,.005)]
            allsol=[]
            for scale,u0 in pars:
              for rho in ([.3,3,10] if kind=='fs' else [1]):
                x=np.array([-.003,np.log(scale),np.log(u0)]+([np.log(rho)] if kind=='fs' else []))
                simplex=np.tile(x,(len(x)+1,1))
                for j in range(len(x)):simplex[j+1,j]+=.003 if j==0 else .2
                sol=minimize(chi,x,args=(kind,),bounds=lim,method='Nelder-Mead',options=dict(initial_simplex=simplex,maxiter=380,xatol=.0001,fatol=.001))
                allsol.append(sol)
            best=min(allsol,key=lambda x:x.fun)
            row[kind+'_chi2']=round(float(best.fun),3)
            row[kind+'_peak_jd']=round(t0+float(best.x[0]),6)
            for param,val in zip((['rise_d','fall_d'] if kind=='flare' else ['te_d','u0','rho']),np.exp(best.x[1:])):
                row[kind+'_'+param]=round(float(val),6)
        results.append(row)
    write_csv(OUT/'moa0211_joint_models.csv',results)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('--joint',action='store_true')
    a=p.parse_args()
    with ThreadPoolExecutor(2) as pool: data=list(pool.map(fields,TARGETS))
    write_csv(OUT/'kmt_field_queries.csv',[v for rr in data for v in rr])
    with ThreadPoolExecutor(2) as pool: write_csv(OUT/'skybot.csv',list(pool.map(moving,TARGETS)))
    if a.joint:joint()
