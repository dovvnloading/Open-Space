"""Query catalogues around the 2025 JWST positions; retain cone answers and per-peak separations."""
from pathlib import Path
import sys,csv,json
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from astropy.time import Time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'stellar'))
from latest_catalogs import rows,cone,sep
import latest_catalogs
latest_catalogs.R=HERE/'fields'

def one(job):
    field,name,refresh=job
    targets=list(csv.DictReader((HERE.parent/'stellar/newest_targets.csv').open()))
    t=next(x for x in targets if x['targprop']==field);ra=float(t['targ_ra']);dec=float(t['targ_dec'])
    services={
        'GaiaDR3':('https://gea.esac.esa.int/tap-server/tap/sync','SELECT source_id,ra,dec,pmra,pmdec,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE ',2016.,'source_id'),
        'SIMBAD':('https://simbad.cfa.harvard.edu/simbad/sim-tap/sync','SELECT main_id,ra,dec,pmra,pmdec,otype FROM basic WHERE ',2000.,'main_id'),
    }
    url,q,epoch,label=services[name]
    cat=rows(field,name,url,q+cone(ra,dec,180),refresh)
    pts=[x for x in csv.DictReader((HERE/'positions.csv').open()) if x['field']==field]
    out=[]
    for p in pts:
        dist=[]
        for x in cat:
            r,d=float(x['ra']),float(x['dec'])
            dt=Time(t['date_obs']).jyear-epoch
            if x['pmra'] and x['pmdec']:
                r+=float(x['pmra'])*dt/(3.6e6*np.cos(np.deg2rad(d)));d+=float(x['pmdec'])*dt/3.6e6
            dist.append(sep(r,d,float(p['ra']),float(p['dec'])))
        if cat:
            i=np.argmin(dist);out.append(dict(field=field,peak=p['peak'],catalog=name,nearest=cat[i][label],separation_arcsec=dist[i]))
        else:out.append(dict(field=field,peak=p['peak'],catalog=name,nearest='',separation_arcsec=''))
    return out

if __name__=='__main__':
    tasks=[(field,name,'--refresh' in sys.argv) for field in ('L-788-34','Ross-619') for name in ('GaiaDR3','SIMBAD')]
    with ThreadPoolExecutor(4) as pool:res=[x for part in pool.map(one,tasks) for x in part]
    with (HERE/'catalogue_matches.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(res[0]));w.writeheader();w.writerows(res)
