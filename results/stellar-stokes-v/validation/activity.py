"""Screen the nearest single-sector 120 s TESS light curve for each dated catalogue host."""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from scipy.ndimage import median_filter

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
# NASA PSCompPars identifiers; record the exact query and raw hash in the receipt.
TIC_QUERY="select distinct hostname,tic_id from pscomppars where hostname in ('AU Mic','GJ 1151','GJ 4274','GJ 896 A','Proxima Cen','COCONUTS-2 A')"
def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cache',type=Path)
    args=p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    def cached(name,url):
        path=args.cache/name
        if not path.exists():
            with urlopen(url,timeout=120) as r:
                raw=r.read()
            path.write_bytes(raw)
        return path
    query_url='https://exoplanetarchive.ipac.caltech.edu/TAP/sync?'+urlencode(dict(query=TIC_QUERY,format='csv'))
    tic_file=cached('host_tics',query_url)
    targets=pd.read_csv(tic_file).set_index('hostname').tic_id.str.split().str[-1].to_dict()
    params={'service':'Mast.Caom.Filtered','params':{'columns':'target_name,obs_id,t_min,t_max,t_exptime,sequence_number,provenance_name,dataURL',
        'filters':[{'paramName':'obs_collection','values':['TESS']},{'paramName':'target_name','values':sorted(targets.values())},
                   {'paramName':'dataproduct_type','values':['timeseries']}]},'format':'json','pagesize':10000,'page':1}
    url='https://mast.stsci.edu/api/v0/invoke?'+urlencode({'request':json.dumps(params)})
    invfile=cached('tess_sorted_inventory.json',url)
    inv=json.loads(invfile.read_bytes())
    if inv['status']!='COMPLETE' or len(inv['data'])>=params['pagesize']:
        raise ValueError('Incomplete MAST response')
    inv=pd.DataFrame(inv['data'])
    inv=inv[(inv.t_exptime==120)&inv.dataURL.str.endswith('_lc.fits',na=False)&(inv.provenance_name=='SPOC')].copy()
    srsc=pd.read_csv(HERE/'srsc_hosts.csv')
    events=[(r.hostname,Time(r.start_utc,scale='utc').mjd,r.exposure_s) for r in srsc[srsc.has_v & srsc.start_utc.notna()].itertuples()]
    racs=pd.read_csv(ROOT/'askap/racs_matches.csv').query("sample=='hosts'")
    raw=pd.read_csv(ROOT/'askap/racs_low2_v.csv').set_index('id')
    for r in racs.itertuples():
        events.append((r.target,r.mjd,raw.loc[r.radio_id,'scan_length']))
    def measure(event):
        host,mjd,exposure=event
        tic=targets[host]
        choices=inv[inv.target_name==tic].copy()
        choices['days_to_sector']=np.maximum(0,np.maximum(choices.t_min-mjd,mjd-choices.t_max))
        ob=choices.sort_values(['days_to_sector','obs_id']).iloc[0]
        fits_url='https://mast.stsci.edu/api/v0.1/Download/file?'+urlencode({'uri':ob.dataURL})
        file=cached(ob.dataURL.rsplit('/',1)[-1],fits_url)
        with fits.open(file) as f:
            tab=f[1].data
            assert str(f[0].header['TICID'])==tic
            tref=f[1].header['BJDREFI']+f[1].header.get('BJDREFF',0)
            good=(tab['QUALITY']==0)&np.isfinite(tab['TIME'])&np.isfinite(tab['PDCSAP_FLUX'])
            t=tref+tab['TIME'][good]
            flux=tab['PDCSAP_FLUX'][good].astype(float)
        starts=np.r_[0,np.flatnonzero(np.diff(t)>3/1440)+1,len(t)]
        optical_events=[]
        good_exposure=0.
        for a,b in zip(starts[:-1],starts[1:]):
            if b-a<73: continue
            x=flux[a:b]
            residual=x/median_filter(x,size=73,mode='nearest')-1
            robust_sigma=1.482602*np.median(abs(residual-np.median(residual)))
            if robust_sigma<=0: continue
            flags=residual>4.5*robust_sigma
            flags[:36]=False
            flags[-36:]=False
            bounds=np.flatnonzero(np.diff(np.r_[False,flags,False]))
            good_exposure+=(b-a-72)*120/86400
            for lo,hi in zip(bounds[::2],bounds[1::2]):
                if hi-lo<2: continue
                i=lo+int(np.argmax(residual[lo:hi]))
                optical_events.append(dict(hostname=host,sector=ob.sequence_number,peak_bjd_tdb=t[a+i],
                     peak_fractional_excess=residual[i],cadences=hi-lo))
        # Observatory MJD and barycentric BJD can differ by several minutes.
        # All six sector windows are disjoint by >100 d; only claim that scoped result.
        overlaps=((choices.t_min<=mjd+float(exposure or 0)/86400)&(choices.t_max>=mjd)).sum() if np.isfinite(exposure) else ((choices.t_min<=mjd)&(choices.t_max>=mjd)).sum()
        row=dict(hostname=host,radio_mjd=mjd,tic=tic,spoc_single_sectors=len(choices),overlap_sectors=int(overlaps),nearest_sector=ob.sequence_number,
            sector_gap_days=ob.days_to_sector,analyzed_exposure_days=good_exposure,candidate_positive_excursions=len(optical_events),
            largest_fractional_excess=max((z['peak_fractional_excess'] for z in optical_events),default=np.nan),fits_url=fits_url,sha256=sha(file.read_bytes()))
        return row,optical_events
    with ThreadPoolExecutor(3) as ex:
        results=list(ex.map(measure,events))
    pd.DataFrame([r[0] for r in results]).sort_values('hostname').to_csv(HERE/'tess_activity.csv',index=False)
    pd.DataFrame([e for _,seq in results for e in seq]).to_csv(HERE/'tess_excursions.csv',index=False)
    (HERE/'tess_receipt.json').write_text(json.dumps(dict(nasa_url=query_url,nasa_sha256=sha(tic_file.read_bytes()),
        mast_url=url,mast_sha256=sha(invfile.read_bytes()),n_spoc_single_sector=len(inv)),indent=2)+'\n')

if __name__=='__main__':
    main()
