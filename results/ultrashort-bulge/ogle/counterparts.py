"""Query four EWS sub-half-day coordinates. Each query records raw response SHA and status."""
import argparse
import concurrent.futures as cf
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tarfile
import urllib.parse
import urllib.request


HERE=Path(__file__).resolve().parent


def coords(ra,dec):
    h,m,s=map(float,ra.split(":")); sign=-1 if dec[0]=="-" else 1
    d,dm,ds=map(abs,map(float,dec.split(":")))
    return 15*(h+m/60+s/3600),sign*(d+dm/60+ds/3600)


def make_queries(row):
    ra,dec=coords(row["ra_hms"],row["dec_dms"])
    y,_,n=row["event"].split("-")
    def tap(base,q):
        return base+"?"+urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q))
    def viz(table,cols,radius):
        pars={'-source':table,'-c':f'{ra:.8f} {dec:+.8f}','-c.rs':radius,'-out.add':'_r','-sort':'_r','-out.max':20,'-out':cols}
        return 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+urllib.parse.urlencode(pars)
    gaia=(f"SELECT source_id,ra,dec,phot_g_mean_mag,bp_rp,phot_variable_flag,ruwe "
          f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.9f},{dec:.9f},0.000555556))")
    fp=f"1=CONTAINS(POINT('ICRS',{ra:.9f},{dec:.9f}),s_region)"
    fields="obs_collection,facility_name,instrument_name,obs_id,dataproduct_type,t_min,t_max"
    return {
        "gaiaDR3_2arcsec":tap('https://gea.esac.esa.int/tap-server/tap/sync',gaia),
        "vvvDR2_2arcsec":viz('II/348/vvv2','srcid,iauname,RAJ2000,DEJ2000,Jmag3,Hmag3,Ksmag3,Jperrbits,Hperrbits,Ksperrbits',2),
        "VSX_5arcsec":viz('B/vsx/vsx','Name,Type,RAJ2000,DEJ2000',5),
        "HST_JWST_CADC_footprint":tap("https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync",
            f"SELECT TOP 30 {fields} FROM ivoa.obscore WHERE {fp} AND (facility_name IN ('HST','JWST') OR obs_collection IN ('HST','JWST','HSTHLA')) AND dataproduct_type='image'"),
        "Euclid_EAS_footprint":tap("https://eas.esac.esa.int/tap-server/tap/sync",
            f"SELECT TOP 30 {fields} FROM ivoa.obscore WHERE {fp} AND dataproduct_type='image'"),
        "ogle_bundle":f"https://www.astrouw.edu.pl/ogle/ogle4/ews/{y}/blg-{n}.tar.gz",
    }


def run_one(task):
    cache,row,kind,url=task
    key=row['event']+'-'+kind
    path=cache/(key+"-"+hashlib.sha256(url.encode()).hexdigest()[:10])
    dest=dict(event=row['event'],product=kind,url=url,queried_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    try:
        if path.exists():
            data=path.read_bytes()
        else:
            with urllib.request.urlopen(url,timeout=100) as resp:
                data=resp.read()
            path.write_bytes(data)
        dest['queried_utc']=datetime.fromtimestamp(path.stat().st_mtime,timezone.utc).isoformat(timespec='seconds')
        dest['sha256']=hashlib.sha256(data).hexdigest()
        if kind=='ogle_bundle':
            with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as tf:
                maps=tf.extractfile(next(m for m in tf.getmembers() if m.name.endswith('/map.dat'))).read()
                rows=[line.split() for line in maps.decode().splitlines() if line.strip()]
                dest['map_star_count']=len(rows)
                dest['map_source_row']=next((v for v in rows if v[0]==row['star_no']),None)
                dest['map_sha256']=hashlib.sha256(maps).hexdigest()
        elif kind.startswith(('vvv','VSX')):
            lines=data.decode().strip().splitlines()
            if not any('VizieR Astronomical Server' in l for l in lines) or any('Error=' in l or 'fatal error' in l.lower() for l in lines):
                raise ValueError('invalid VizieR response')
            j=next((j for j,l in enumerate(lines) if re.match(r'^[\-\t ]{3,}$',l)),None)
            dest['rows']=[l.strip().split('\t') for l in lines[j+1:] if l.strip() and not l.startswith('#')] if j is not None else []
            dest['header']=lines[j-2] if j is not None else None
        else:
            lines=data.decode().strip().splitlines()
            if not lines or lines[0].lower().startswith(('<','error')): raise ValueError('TAP response not CSV')
            dest['rows']=list(csv.DictReader(lines))
            if kind=='gaiaDR3_2arcsec':
                a,d=coords(row['ra_hms'],row['dec_dms'])
                for v in dest['rows']:
                    v['sep_arcsec']=round(3600*math.hypot((a-float(v['ra']))*math.cos(math.radians(d)),d-float(v['dec'])),3)
        dest['status']='ok'
    except Exception as exc:
        dest['status']=str(exc)
    return dest


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('cache',type=Path)
    cache=ap.parse_args().cache
    cache.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((HERE/'catalog.csv').open()))
    jobs=[(cache,r,k,v) for r in rows if 0<float(r['tau_days'])<.5 for k,v in make_queries(r).items()]
    with cf.ThreadPoolExecutor(8) as ex:
        records=list(ex.map(run_one,jobs))
    (HERE/'counterparts.json').write_text(json.dumps(records,indent=2)+'\n')

if __name__=='__main__':
    main()
