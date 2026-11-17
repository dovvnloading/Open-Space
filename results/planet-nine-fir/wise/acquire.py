"""Query the final two calendar years of NEOWISE for W2 images and select visits."""
import argparse, csv, hashlib, io, json, math, sys, time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import numpy as np
from astropy.io import fits
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'other'))
from candidates import act_rows,act_sky,epoch
from other_fields import fields,position

def load_url(url):
    for i in range(3):
        try:
            with urlopen(url,timeout=90) as f:return f.read()
        except Exception:
            if i==2:raise
            time.sleep(1+i)

def run(target,ra,dec,cache,down):
    home=cache/target;home.mkdir(exist_ok=True,parents=True)
    url='https://irsa.ipac.caltech.edu/SIA?'+urlencode(dict(COLLECTION='neowiser',TIME='59945 60532',POS=f'CIRCLE {ra:.5f} {dec:.5f} 0.04',BAND='4.6e-6',RESPONSEFORMAT='CSV',MAXREC=30000));meta=home/'query.csv'
    if not meta.exists():meta.write_bytes(load_url(url))
    raw=meta.read_bytes();rows=list(csv.DictReader(io.StringIO(raw.decode())));science={r['obs_id']:r for r in rows if r.get('energy_bandpassname')=='W2' and r['dataproduct_subtype']=='science' and '-w2-int-' in r['access_url']}
    if not rows or 'obs_id' not in rows[0]:raise RuntimeError(raw[:250])
    groups=[]
    for r in sorted(science.values(),key=lambda r:float(r['t_min'])):
        if not groups or float(r['t_min'])-float(groups[-1][-1]['t_min'])>30:groups.append([])
        groups[-1].append(r)
    receipt={'target':target,'query':url,'query_sha256':hashlib.sha256(raw).hexdigest(),'rows':len(rows),'science_count':len(science),'visits':[len(g) for g in groups],'center':[ra,dec],'products':[],'retrieved_utc':datetime.now(timezone.utc).isoformat()}
    def distance(r): return math.hypot((float(r['s_ra'])-ra)*math.cos(math.radians(dec)),float(r['s_dec'])-dec)
    for j,g in enumerate(groups):
        # Choose distinct exposures nearest common sky center, without looking at pixels.
        for r in sorted(g,key=distance)[:(8 if target=='Phan' else 3)]:
            mid=(float(r['t_min'])+float(r['t_max']))/2
            record={'visit':j,'id':r['obs_id'],'mjd':mid,'url':r['access_url'],'sep_deg':round(distance(r),6)};receipt['products'].append(record)
    prior_path=HERE/'receipts'/(target+'.json')
    prior={} if not prior_path.exists() else {p['id']:p for p in json.loads(prior_path.read_text())['products']}
    def image(p):
        # One common tile across all visits for stationary-background differencing.
        size=570 if target=='Phan' else 660
        for kind in ['int','msk']:
            path=home/(p['id']+'_'+kind+'.fits');base=p['url'].replace('-int-','-msk-')+'.gz' if kind=='msk' else p['url']
            u=base+'?'+urlencode(dict(center=f'{ra:.5f},{dec:.5f}',size=f'{size}pix',gzip='false'))
            p[kind+'_url']=u
            previous=prior.get(p['id'],{})
            if path.exists() and kind+'_original_url' in previous and hashlib.sha256(path.read_bytes()).hexdigest()==previous.get(kind+'_sha256'):
                p.update({k:v for k,v in previous.items() if k.startswith(kind+'_')});continue
            try:
                if not path.exists():
                    data=load_url(u)
                    if not data.startswith(b'SIMPLE'):raise ValueError(repr(data[:90]))
                    path.write_bytes(data)
                h=fits.getheader(path);p[kind+'_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();p[kind+'_shape']=[h['NAXIS2'],h['NAXIS1']]
            except Exception as e:p[kind+'_error']=repr(e)
    if down:
        with ThreadPoolExecutor(5 if target=='Phan' else 2) as ex:list(ex.map(image,receipt['products']))
    out=HERE/'receipts';out.mkdir(exist_ok=True);(out/(target+'.json')).write_text(json.dumps(receipt,indent=1)+'\n');print(target,receipt['visits'],len(receipt['products']),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--targets',nargs='+',default=['Phan']);p.add_argument('--download',action='store_true');args=p.parse_args()
    sites={'Phan':(35.05,-49.63)}
    for k,r in act_rows().items():
        a=r['ra']+r['vx']*(60200-57688)/365.25/60/np.cos(np.radians(r['dec']));d=r['dec']+r['vy']*(60200-57688)/365.25/60;sites[k]=(a%360,d)
    for key in fields():
        # Center midway between opposite seasons for this nominal circular solution.
        coords=[position(key,epoch(day)) for day in ['2023-03-01','2023-09-01','2024-03-01','2024-09-01']]
        sites[key]=tuple(np.mean(np.array(coords),axis=0).ravel())
    names=list(act_rows()) if args.targets==['ACT'] else (fields() if args.targets==['other'] else args.targets)
    def one(k):
        try:run(k,*sites[k],args.cache,args.download)
        except Exception as e:print(k,'FAILED',repr(e),flush=True);raise
    with ThreadPoolExecutor(4) as ex:list(ex.map(one,names))
if __name__=='__main__':main()
