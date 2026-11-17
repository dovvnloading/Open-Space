"""Add native frames to fields with no eligible initial tests, using geometry only."""
import argparse,csv,gzip,hashlib,json,math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from astropy.io import fits
from acquire import HERE, load_url

def extend(target,cache):
    home=cache/target;ledger=HERE/'receipts'/(target+'.json');r=json.loads(ledger.read_text());ra,dec=r['center'];base={p['id']:p for p in r['products']}
    sky=[x for x in csv.DictReader((home/'query.csv').open()) if x.get('energy_bandpassname')=='W2' and '-w2-int-' in x['access_url']]
    groups=[]
    for x in sorted(sky,key=lambda x:float(x['t_min'])):
        if not groups or float(x['t_min'])-float(groups[-1][-1]['t_min'])>30:groups.append([])
        groups[-1].append(x)
    new=[]
    for i,g in enumerate(groups):
        def dist(x):return math.hypot((float(x['s_ra'])-ra)*math.cos(math.radians(dec)),float(x['s_dec'])-dec)
        scans=set()
        for x in sorted(g,key=dist):
            scan=x['obs_id'].split('_')[1]
            if scan in scans:continue
            scans.add(scan)
            if x['obs_id'] not in base:
                new.append(dict(id=x['obs_id'],visit=i,mjd=(float(x['t_min'])+float(x['t_max']))/2,url=x['access_url'],sep_deg=round(dist(x),6)))
            if len(scans)==10:break
    def fetch(p):
        for k in ('int','msk'):
            part=p['url'].split('/p1bm_frm/')[1].replace('-int-','-'+k+'-')+('.gz' if k=='msk' else '')
            url='https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/neowiser/images/p1bm_frm/'+part;p[k+'_url']=url;p[k+'_original_url']=url
            try:
                raw=load_url(url);data=gzip.decompress(raw) if raw[:2]==bytes([31,139]) else raw
                if not data.startswith(b'SIMPLE'):raise ValueError(repr(data[:60]))
                dest=home/(p['id']+'_'+k+'.fits');dest.write_bytes(data);h=fits.getheader(dest);p[k+'_shape']=[h['NAXIS2'],h['NAXIS1']];p[k+'_sha256']=hashlib.sha256(data).hexdigest();p[k+'_original_sha256']=hashlib.sha256(raw).hexdigest()
            except Exception as e:p[k+'_fallback_error']=repr(e)
    with ThreadPoolExecutor(4) as ex:list(ex.map(fetch,new))
    r['products'].extend(new);r['extended_selection']='up to 10 distinct scan IDs per visit for blank eligibility failure';ledger.write_text(json.dumps(r,indent=1)+'\n')
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',type=Path,required=True);a.add_argument('target',nargs='+');p=a.parse_args()
    for target in p.target:extend(target,p.cache)
