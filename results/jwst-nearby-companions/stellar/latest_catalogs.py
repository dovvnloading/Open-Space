"""Select host IDs, motions and F444W normalization for the latest public GO 8826 fields."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from urllib.parse import urlencode
from urllib.request import urlopen
import numpy as np
from astropy.time import Time

HERE=Path(__file__).resolve().parent
R=HERE/'latest_catalogs'
def read_csv(path):
    with path.open() as f:return list(csv.DictReader(f))

def sep(ra,dec,qra,qdec):
    return 3600*np.hypot((ra-qra)*np.cos(np.deg2rad(dec)),dec-qdec)

def rows(tag,name,url,q,refresh):
    dest=R/f"{tag}_{name}.csv"; receipt=dest.with_suffix('.json')
    if refresh:
        url+="?"+urlencode(dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=q))
        with urlopen(url,timeout=60) as resp:raw=resp.read()
        if b'<' in raw[:4] or not raw:raise ValueError(raw[:120])
        dest.write_bytes(raw)
        receipt.write_text(json.dumps(dict(url=url,query=q,sha256=hashlib.sha256(raw).hexdigest(),retrieved_utc=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
    b=dest.read_bytes()
    if hashlib.sha256(b).hexdigest()!=json.loads(receipt.read_text())['sha256']:raise ValueError(dest)
    return list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))))

def cone(r,d,radius):
    return f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r:.9f},{d:.9f},{radius/3600:.9f}))"

def one(job):
    t,refresh=job
    name=t['targprop'];ra=float(t['targ_ra']);dec=float(t['targ_dec'])
    catalogs=rows(name,'simbad','https://simbad.cfa.harvard.edu/simbad/sim-tap/sync',
      "SELECT main_id,ra,dec,pmra,pmdec,plx_value,sp_type,otype FROM basic WHERE "+cone(ra,dec,240),refresh)
    # Use the nearest object after proper motion from J2000.0 to the MAST target epoch.
    def pos(x,years):
        cd=np.cos(np.deg2rad(float(x['dec'])))
        return float(x['ra'])+years*float(x['pmra'] or 0)/(3.6e6*cd),float(x['dec'])+years*float(x['pmdec'] or 0)/3.6e6
    date=Time(t['date_obs']).jyear
    host=min((x for x in catalogs if x['otype'] not in ('Pl','Pl?')),key=lambda x:sep(*pos(x,date-2000),ra,dec))
    off=sep(*pos(host,date-2000),ra,dec)
    if off>2.:raise ValueError((name,"host separation",off))
    wr,wd=pos(host,10.5)
    wises=rows(name,'allwise','https://irsa.ipac.caltech.edu/TAP/sync',
      "SELECT designation,ra,dec,w2mpro,w2sigmpro,w2sat,w2snr,cc_flags,ph_qual FROM allwise_p3as_psd WHERE "+cone(wr,wd,20),refresh)
    wise=min(wises,key=lambda x:sep(float(x['ra']),float(x['dec']),wr,wd))
    ws=sep(float(wise['ra']),float(wise['dec']),wr,wd)
    if ws>2.:raise ValueError((name,"AllWISE separation",ws))
    kr,kd=pos(host,0.)
    ks_sources=rows(name,'twomass','https://irsa.ipac.caltech.edu/TAP/sync',
       "SELECT designation,ra,dec,k_m,k_cmsig,ph_qual,cc_flg FROM fp_psc WHERE "+cone(kr,kd,20),refresh)
    ks=min(ks_sources,key=lambda x:sep(float(x['ra']),float(x['dec']),kr,kd))
    gr,gd=pos(host,16.)
    gaia=rows(name,'gaia','https://gea.esac.esa.int/tap-server/tap/sync',
      "SELECT source_id,ra,dec,pmra,pmdec,parallax,phot_g_mean_mag,ruwe FROM gaiadr3.gaia_source WHERE "+cone(gr,gd,45),refresh)
    g=min(gaia,key=lambda x:sep(float(x['ra']),float(x['dec']),gr,gd))
    gsep=sep(float(g['ra']),float(g['dec']),gr,gd)
    planets=rows(name,'nasa_ps','https://exoplanetarchive.ipac.caltech.edu/TAP/sync',"SELECT DISTINCT pl_name,hostname,ra,dec FROM ps WHERE "+cone(float(host['ra']),float(host['dec']),90),refresh)
    grid=read_csv(HERE/'mcc124_sed.csv')
    # Reuse PHOENIX colour grid, calibrated to actual per-host W2 and Ks magnitudes.
    # Difference from MCC 124's reference magnitude is independent of the model.
    orig=float(read_csv(HERE/'catalogue_receipts/allwise.csv')[0]['w2mpro'])
    origk=float(read_csv(HERE/'catalogue_receipts/twomass.csv')[0]['k_m'])
    f=[float(z['w2_anchored_f444w_mJy'])/1000*10**(-.4*(float(wise['w2mpro'])-orig)) for z in grid]
    fk=[float(z['ks_anchored_f444w_mJy'])/1000*10**(-.4*(float(ks['k_m'])-origk)) for z in grid]
    return dict(target=name,host=host['main_id'],sp_type=host['sp_type'],simbad_host_sep_arcsec=off,
                pm_east_masyr=host['pmra'],pm_north_masyr=host['pmdec'],parallax_mas=host['plx_value'],
                gaia_id=g['source_id'],gaia_sep_at_2016_arcsec=gsep,gaia_rows_45arcsec=len(gaia),wise_id=wise['designation'],
                wise_sep_at_2010p5_arcsec=ws,W2=wise['w2mpro'],W2_err=wise['w2sigmpro'],W2_sat=wise['w2sat'],
                W2_cc_flags=wise['cc_flags'],Ks=ks['k_m'],Ks_sep_J2000_arcsec=sep(float(ks['ra']),float(ks['dec']),kr,kd),
                F444W_Jy=float(np.median(f)),W2grid_min_Jy=min(f),W2grid_max_Jy=max(f),Ksgrid_min_Jy=min(fk),Ksgrid_max_Jy=max(fk),
                nasa_ps_rows_90arcsec=len(planets))

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--refresh',action='store_true')
    a=parser.parse_args();R.mkdir(exist_ok=True)
    with ThreadPoolExecutor(4) as pool:
        out=list(pool.map(one,[(t,a.refresh) for t in read_csv(HERE/'newest_targets.csv')]))
    with (HERE/'latest_hosts.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(out[0]));w.writeheader();w.writerows(out)

if __name__=="__main__":main()
